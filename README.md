# PALACE
PALACE is a graph-based algorithm to reconstruct the phage genome. PALACE currently supports normal pair-end reads, Oxford Nanopore(ONT) and PacBio SMRT(PB) reads. Please check xxx for detail.
## Installation
### Approach 1, install with mamba/conda.
1. Clone the repository and enter the directory:

```
git clone https://github.com/deepomicslab/PALACE
cd ./PALACE
```
2. Create a conda environment with all dependencies and enter the environment:
```
mamba env create --prefix=./PALACE -f environment.yml # conda env create --prefix=./PALACE -f environment.yml
mamba activate ./PALACE # conda activate ./PALACE
```
3. Create a build directory and compile PALACE under it (use **sudo**, if required):

```
cd seqgraph_phage/
mkdir build
cd build
cmake ..
make && make install
```
### Approach 2, from scratch
### Prerequisites
### Python packages
* pysam==0.17.0
* numpy==1.20.2
* sklearn==1.1.1
* biopython==1.78
* pysam==0.17.0
* matplotlib==3.4.2
### Torch packages, (cpu or gpu)
Please check https://pytorch.org/get-started/previous-versions/ for installation
* torch==1.7.1
* torch-cluster==1.5.9
* torch-geometric==1.7.0
* torch-scatter==2.0.6
* torch-sparse==0.6.9
* torch-spline-conv==1.2.1
* torch-summary==1.4.5
* torchvision==0.8.0a0
### Other packages
* [bwa](https://github.com/lh3/bwa) BWA is a software package for reads mapping.
* [samtools](http://www.htslib.org/download/) Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format.
* [fastp](https://github.com/OpenGene/fastp) Provide fast all-in-one preprocessing for FastQ files.
* [spades](https://github.com/ablab/spades) Genome assembler
* [ncbi-blast](https://www.ncbi.nlm.nih.gov/books/NBK569861/) Sequence alignment tool.
* [htslib](http://www.htslib.org/download/)

Install the prerequisites first, clone the repository and enter the directory:
```
git clone https://github.com/deepomicslab/PALACE
#create a new mamba(conda) env
mamba create -n palace ## or conda create -n palace
mamba activate palace ## or conda activate palace
cd ./PALACE/seqGraph_phage/
mkdir build
cd build
cmake ..
make
```


## Using PALACE
1. Config the config.txt file, [here](https://github.com/deepomicslab/PALACE/blob/main/config.txt) is a demo file.
* ```fastq1```, Read1 paired fastq file.
* ```fastq2```, Read2 paired fastq file.
* ```phagedb```, Phage reference database, this is a phage reference database, the latest database can be download from [here](https://portal.nersc.gov/CheckV/checkv-db-v1.5/genome_db/checkv_reps.fna).
* ```protein_db```, Phage protein database, this is a phage protein database, the latest database can be download from [here](https://portland-my.sharepoint.com/:f:/g/personal/gzpan2-c_my_cityu_edu_hk/EpVA0ISAp4FBrclyldwpjEwBBHujF4zG2Gu3Vxa5AZICJw?e=5z2qUe).
* ```threads```, Threads to be used.
* ```out_dir```, Output directory.
* ```prefix```, Intermediate file prefix, can be sample name.
* ```PYTHON```, Python path.
* ```BWA```, bwa path.
* ```SAMTOOLS```, samtools path.
* ```FASTP```, fastp path.
* ```SPADES```, spades.py path.
* ```NCBI_BIN```, ncbi-blast bin path, must contains makeblastdb, blastn and tblastn 
* ```PALACE```, PALACE path.
2. Runing PALACE.

* ``` bash PALACE_PATH/pipe.sh config.txt```
## Output
* ```01-qc/```, fastp output.
* ```02-assembly/```, Raw assembly result with spades with --meta.
* ```03-search/```, This directory contains three main intermediate files: ```hit_seqs.out``` contains contigs with phage protein. ```node_scores.out```, the second column is the score predicted by deeplearning network. ```{prefix}_ref_names.txt```, contains phage references identified by kmer alignment. 
* ```04-match/```, This directory contains the graph structure of the conjugate graph(```{prefix}_filtered_graph.txt```), the results of the graph decompose(```{prefix}_all_result.txt```).
* ```05-furth```, This directory contains the local matching result based on the phage reference.
* ```final_result```, This directory contains the final result, final contig paths for phages(```{prefix}_final.txt```), cycle paths for phages(```{prefix}_cycle.txt```), phages fasta(```{prefix}_final.fasta)
## Author
PALACE is developed by DeepOmics lab under the supervision of Dr. Li Shuaicheng, City University of Hong Kong, Hong Kong, China. Should you have any queries, please feel free to contact us by gzpan2-c@my.cityu.edu.hk or ruhawang2@my.c.
## License
This project is licensed under the MIT License - see the [LICENSE.txt](https://github.com/deepomicslab/PALACE/blob/main/LICENSE.txt) file for details.
