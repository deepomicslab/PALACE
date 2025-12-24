import sys
import re
import subprocess
from collections import defaultdict
from operator import itemgetter


class SegmentProcessor:
    """Process sequence segments and graph data"""
    
    def __init__(self, samtools_path):
        self.samtools = samtools_path
        self.seg_graph = defaultdict(list)
        self.segs = []
        self.ref_name_records = {}
        self.ref_segs = {}
        
    def parse_blast(self, blast_in, blast_ratio, refs):
        """Parse BLAST output to identify reference segments"""
        ref_blast_segs = {ref: set() for ref in refs}
        
        prev_seg = ""
        prev_len = 0
        prev_ref = ""
        prev_seg_len = 0
        
        for line in blast_in:
            parts = line.strip().split("\t")
            if len(parts) < 6 or parts[1] not in refs:
                continue
                
            curr_seg = parts[0]
            curr_ref = parts[1]
            curr_identity = float(parts[2])
            curr_align_len = int(parts[5])
            curr_seg_len = int(parts[3])
            
            # Check if we need to process previous segment
            if (prev_seg != curr_seg and prev_seg) or (prev_ref != curr_ref and prev_ref):
                if prev_seg_len > 0 and prev_len / prev_seg_len > blast_ratio:
                    ref_blast_segs[prev_ref].add(prev_seg)
                prev_len = curr_align_len
                prev_seg_len = curr_seg_len
            else:
                if curr_identity > blast_ratio * 100:
                    prev_len += curr_align_len
                    
            prev_seg = curr_seg
            prev_ref = curr_ref
            prev_seg_len = curr_seg_len
        
        # Process last segment
        if prev_seg and prev_seg_len > 0 and prev_len / prev_seg_len > blast_ratio:
            ref_blast_segs[prev_ref].add(prev_seg)
            
        return ref_blast_segs
    
    def run_samtools_depth(self, input_bam, region):
        """Run samtools depth command for a specific region"""
        cmd = f'{self.samtools} depth -r {region} {input_bam}'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error running samtools depth: {result.stderr}")
            return []
            
        depths = []
        for line in result.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) >= 3:
                depths.append(int(parts[2]))
        return depths
    
    def calculate_average_depth(self, depths):
        """Calculate average depth from a list of depths"""
        return sum(depths) / len(depths) if depths else 0
    
    def get_depth_segments(self, segs, bam):
        """Calculate depth for segments and determine copy numbers"""
        total_depths = []
        average_contig_depth = {}
        
        # Calculate depth for each segment
        for seg in segs:
            depths = self.run_samtools_depth(bam, seg)
            if depths:
                avg_depth = self.calculate_average_depth(depths)
                average_contig_depth[seg] = avg_depth
                total_depths.extend(depths)
        
        # Calculate total average depth
        total_avg_depth = self.calculate_average_depth(total_depths)
        if total_avg_depth == 0:
            total_avg_depth = 1  # Avoid division by zero
        
        # Generate final segments with copy numbers
        final_segs = []
        for seg, avg_depth in average_contig_depth.items():
            copy_num = max(1, round(avg_depth / total_avg_depth))
            cov_value = seg.split('_')[-1]
            final_segs.append(f"SEG {seg} {cov_value} {copy_num}")
            
        return final_segs
    
    def load_segments(self, seg_file):
        """Load segments from file"""
        self.seg_graph.clear()
        self.segs.clear()
        self.ref_name_records.clear()
        self.ref_segs.clear()
        
        for idx, line in enumerate(seg_file):
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
                
            # Extract segments (remove direction indicators)
            segments = re.split(r'[+-]', parts[0])[:-1]
            ref_name = parts[1]
            
            self.ref_name_records[idx] = ref_name
            self.ref_segs[ref_name] = set(segments)
            
            # Build segment graph
            line_segs = []
            for seg in segments:
                if seg:
                    line_segs.append(seg)
                    self.seg_graph[seg].append(idx)
                    
            self.segs.append(line_segs)
    
    def find_segment_indices(self, seg1, seg2=None):
        """Find indices where segments appear"""
        if seg2 is None:
            # Check if seg1 exists in any segment list
            for seg_list in self.segs:
                if seg1 in seg_list:
                    return 1
            return 0
        else:
            # Find indices where both seg1 and seg2 exist
            indices = []
            for i, seg_list in enumerate(self.segs):
                if seg1 in seg_list and seg2 in seg_list:
                    indices.append(i + 1)
            return indices
    
    def process_graph(self, graph_file, min_support, bam, output_prefix):
        """Process graph file and generate outputs"""
        # Initialize output structures
        num_refs = len(self.ref_name_records)
        out_segs = [[] for _ in range(num_refs)]
        out_juncs = [[] for _ in range(num_refs)]
        
        # Process graph file
        for line in graph_file:
            parts = line.rstrip().split(" ")
            if not parts:
                continue
                
            if parts[0] == "SEG":
                # Process segment line
                if len(parts) > 1 and self.find_segment_indices(parts[1]):
                    for idx in self.seg_graph.get(parts[1], []):
                        out_segs[idx].append(line)
            else:
                # Process junction line
                if len(parts) >= 6:
                    seg1, seg3 = parts[1], parts[3]
                    indices = self.find_segment_indices(seg1, seg3)
                    for idx in indices:
                        out_juncs[idx - 1].append(parts)
        
        # Generate output files
        for idx, ref_name in self.ref_name_records.items():
            self._write_output_file(idx, ref_name, out_juncs[idx], 
                                  min_support, bam, output_prefix)
    
    def _write_output_file(self, idx, ref_name, juncs, min_support, bam, output_prefix):
        """Write output file for a specific reference"""
        output_file = f"{output_prefix}_{idx}ref{ref_name}ref.second"
        
        # Expand reference segments based on junctions
        ref_seg_set = self.ref_segs[ref_name].copy()
        prev_count = 0
        curr_count = len(ref_seg_set)
        
        while prev_count != curr_count:
            prev_count = curr_count
            for junc in juncs:
                if len(junc) >= 6 and int(junc[-1]) >= min_support:
                    seg1, seg3 = junc[1], junc[3]
                    if seg1 in ref_seg_set or seg3 in ref_seg_set:
                        ref_seg_set.add(seg1)
                        ref_seg_set.add(seg3)
            curr_count = len(ref_seg_set)
        
        # If no junctions, use original segments
        if not juncs:
            ref_seg_set = set(self.segs[idx])
        
        # Write output
        with open(output_file, 'w') as f:
            # Write segments with depth correction
            depth_corrected_segs = self.get_depth_segments(ref_seg_set, bam)
            for seg in depth_corrected_segs:
                f.write(f"{seg}\n")
            
            # Write valid junctions
            for junc in sorted(juncs):
                if len(junc) >= 6 and int(junc[-1]) >= min_support:
                    seg1, seg3 = junc[1], junc[3]
                    if seg1 in ref_seg_set or seg3 in ref_seg_set:
                        f.write(f"{' '.join(junc)}\n")


def main():
    """Main function to process graph and segment data"""
    # Parse command line arguments
    if len(sys.argv) < 5:
        print("Usage: script.py graph_file output_prefix seg_file samtools_path "
              "[min_support] [bam_file] [blast_ratio]")
        sys.exit(1)
    
    graph_file = sys.argv[1]
    output_prefix = sys.argv[2]
    seg_file = sys.argv[3]
    samtools_path = sys.argv[4]
    min_support = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    bam_file = sys.argv[6] if len(sys.argv) > 6 else None
    blast_ratio = float(sys.argv[7]) if len(sys.argv) > 7 else 0.7
    
    # Initialize processor
    processor = SegmentProcessor(samtools_path)
    
    # Load segments
    with open(seg_file, 'r') as f:
        processor.load_segments(f)
    
    # Process graph
    with open(graph_file, 'r') as f:
        processor.process_graph(f, min_support, bam_file, output_prefix)


if __name__ == "__main__":
    main()
