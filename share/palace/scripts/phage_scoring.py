import numpy as np
import pandas as pd
import sys
import os
import types
from Bio import SeqIO
from multiprocessing import Pool
from functools import partial

# GPU settings: normally configured by SLURM, setting to GPU 0
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

import torch
import torch.nn as nn
import torch.nn.functional as F

# Import PyG modules
import torch_geometric.nn as pyg_nn
import torch_geometric.utils as pyg_utils
#from torch_geometric.loader import DataLoader
from torch_geometric.data import DataLoader
from torch_geometric.data import Data

import encode  # Ensure encode.pyx is recompiled for new environment

# ==========================================
# Patch 1: Fix model loading path compatibility (Pickle patch)
# ==========================================
import torch_geometric.utils

class MockInspector:
    def __init__(self, *args, **kwargs): pass
    def __getattr__(self, name): return lambda *args, **kwargs: None
    def inspect(self, *args, **kwargs): return None

# Inject old paths to prevent torch.load from failing to find classes
torch_geometric.utils.Inspector = MockInspector
mock_mod = types.ModuleType("inspector")
mock_mod.Inspector = MockInspector
sys.modules['torch_geometric.nn.conv.utils.inspector'] = mock_mod
sys.modules['torch_geometric.nn.conv.utils'] = types.ModuleType("utils")
sys.modules['torch_geometric.nn.conv.utils'].inspector = mock_mod

# ==========================================
# Model architecture definition
# ==========================================
HIDDEN_DIM = 3
PNODE_DIM = 3
PNODE_NUM = 4096
FNODE_NUM = 64
GCN_HIDDEN_DIM = 128
CNN_HIDDEN_DIM = 64
FC_HIDDEN_DIM = 100
GCN_LAYER_NUM = 2
DROP_RATE = 0.2

class GNN_Model(nn.Module):
    def __init__(self):
        super(GNN_Model, self).__init__()
        self.gcn_dim = GCN_HIDDEN_DIM
        self.num_layers = GCN_LAYER_NUM
        self.dropout = DROP_RATE

        # 1. Definitions corresponding to unexpected keys
        self.pnode_d = nn.Linear(PNODE_NUM * PNODE_DIM, PNODE_NUM * HIDDEN_DIM)
        self.fnode_d = nn.Linear(FNODE_NUM, FNODE_NUM * HIDDEN_DIM)

        # 2. Layer names must be convs_1 and convs_2
        self.convs_1 = nn.ModuleList([
            pyg_nn.SAGEConv((HIDDEN_DIM, HIDDEN_DIM), self.gcn_dim),
            pyg_nn.SAGEConv((self.gcn_dim, self.gcn_dim), self.gcn_dim)
        ])
        self.convs_2 = nn.ModuleList([
            pyg_nn.SAGEConv((self.gcn_dim, HIDDEN_DIM), self.gcn_dim),
            pyg_nn.SAGEConv((self.gcn_dim, self.gcn_dim), self.gcn_dim)
        ])

        # 3. Corresponds to lns.0
        self.lns = nn.ModuleList([nn.LayerNorm(self.gcn_dim)])

        # 4. Corresponds to conv1, conv2, conv3, d1, d2
        self.conv1 = nn.Conv1d(in_channels=self.gcn_dim, out_channels=64, kernel_size=8)
        self.conv2 = nn.Conv1d(in_channels=64, out_channels=64, kernel_size=8)
        self.conv3 = nn.Conv1d(in_channels=64, out_channels=64, kernel_size=8)
        self.d1 = nn.Linear(4075 * 64, 100)
        self.d2 = nn.Linear(100, 2)

    def forward(self, data):
        x_f, x_p = data.x_src, data.x_dst
        edge_index_forward = data.edge_index[:, ::2]
        edge_index_backward = data.edge_index[[1, 0], :][:, 1::2]

        x_p = torch.reshape(x_p, (-1, PNODE_NUM * PNODE_DIM))
        x_p = self.pnode_d(x_p)
        x_p = torch.reshape(x_p, (-1, HIDDEN_DIM))

        x_f = torch.reshape(x_f, (-1, FNODE_NUM))
        x_f = self.fnode_d(x_f)
        x_f = torch.reshape(x_f, (-1, HIDDEN_DIM))

        for i in range(self.num_layers):
            x_p = self.convs_1[i]((x_f, x_p), edge_index_forward)
            x_p = F.relu(x_p)
            x_p = F.dropout(x_p, p=self.dropout, training=self.training)
            x_f = self.convs_2[i]((x_p, x_f), edge_index_backward)
            x_f = F.relu(x_f)
            x_f = F.dropout(x_f, p=self.dropout, training=self.training)
            if i < self.num_layers - 1:
                x_p = self.lns[i](x_p)
                x_f = self.lns[i](x_f)

        x = torch.reshape(x_p, (-1, self.gcn_dim, PNODE_NUM))
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = F.relu(self.conv3(x))
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = x.flatten(start_dim=1)
        x = F.relu(self.d1(x))
        return F.softmax(self.d2(x), dim=1)

# Patch 2: Adapt to new PyG __inc__ interface
class BipartiteData(Data):
    def __inc__(self, key, value, *args, **kwargs):
        if key == 'edge_index':
            return torch.tensor([[self.x_src.size(0)], [self.x_dst.size(0)]])
        else:
            return super(BipartiteData, self).__inc__(key, value, *args, **kwargs)

def make_edge():
    edge = []
    for i in range(4096):
        edge.append([i // 64, i]); edge.append([i % 64, i])
    return np.array(edge).T

def generate_model_input(fasta_file, K, thread, batch_size=1000):
    pool = Pool(thread)
    partial_encode_seq = partial(encode.matrix_encoding, K=K)
    seq_batch, data_id_batch = [], []

    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq_batch.append(str(seq_record.seq))
        data_id_batch.append(seq_record.id)
        if len(seq_batch) == batch_size:
            feature = np.array(pool.map(partial_encode_seq, seq_batch))
            yield data_id_batch, feature
            seq_batch, data_id_batch = [], []

    if seq_batch:
        feature = np.array(pool.map(partial_encode_seq, seq_batch))
        yield data_id_batch, feature
    pool.close(); pool.join()

# ==========================================
# Main program
# ==========================================
if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python script.py <fasta> <out> <reverse> <threads> <model_path> [batch_size]")
        sys.exit(1)

    input_fasta, output_file, _, thread_num, model_path = sys.argv[1:6]
    batch_size = int(sys.argv[6]) if len(sys.argv) > 6 else 1000
    thread = int(thread_num)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    try:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
        model = GNN_Model()
        
        checkpoint = torch.load(model_path, map_location=device)
        
        if hasattr(checkpoint, 'state_dict'):
            state_dict = checkpoint.state_dict()
        else:
            state_dict = checkpoint
            
        model.load_state_dict(state_dict, strict=False)
        model.to(device)
        model.eval()
        print(f"Model loaded on {device} and patched.")
    except Exception as e:
        print(f"Load Error: {e}"); sys.exit(1)

    f_out = open(output_file, "w")
    first_line = True
    edge_template = make_edge()

    for data_id_batch, input_data in generate_model_input(input_fasta, 3, thread, batch_size):
        pnode_feature = input_data.reshape(-1, 3, 4096)
        pnode_feature = np.moveaxis(pnode_feature, 1, 2)
        zero_layer = input_data.reshape(-1, 3, 64, 64)[:, 0, :, :]
        fnode_feature = np.sum(zero_layer, axis=2).reshape(-1, 64, 1)

        data_list = []
        for i in range(pnode_feature.shape[0]):
            x_p = torch.tensor(pnode_feature[i], dtype=torch.float)
            x_f = torch.tensor(fnode_feature[i], dtype=torch.float)
            data = BipartiteData(x_src=x_f, x_dst=x_p, 
                                 edge_index=torch.tensor(edge_template, dtype=torch.long),
                                 num_nodes=x_f.size(0) + x_p.size(0))
            data_list.append(data)

        loader = DataLoader(data_list, batch_size=64, shuffle=False, follow_batch=['x_src', 'x_dst'])
        
        curr_idx = 0
        for batch in loader:
            with torch.no_grad():
                batch = batch.to(device)
                out = model(batch)
                probs = out[:, 1].cpu().numpy()
                for p in probs:
                    if not first_line: f_out.write("\n")
                    f_out.write(f"{data_id_batch[curr_idx]}\t{p}")
                    first_line = False
                    curr_idx += 1
        f_out.flush()

    f_out.close()
    print("Done.")
