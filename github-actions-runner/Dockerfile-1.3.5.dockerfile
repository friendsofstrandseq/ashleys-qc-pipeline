Config file config/config.yaml is extended by additional config specified via the command line.
['/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/config/fastqc_output_touch.txt', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/cell_selection/labels.tsv', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/counts/CountComplete.classic.pdf', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/counts/CountComplete.GC_corrected.pdf', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/config/alfred_output_touch.txt', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/MERGE/merged_bam_gc_dist.merge.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/MERGE/merged_bam_gc_devi.merge.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/A_gc_dist.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/A_gc_devi.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/B_gc_dist.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/B_gc_devi.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/C_gc_dist.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/C_gc_devi.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/D_gc_dist.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/D_gc_devi.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/E_gc_dist.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/E_gc_devi.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/F_gc_dist.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/F_gc_devi.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/G_gc_dist.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/G_gc_devi.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/H_gc_dist.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/alfred/PLATE_ROW/H_gc_devi.row.png', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/plate/ashleys_plate_predictions.pdf', '/scratch/tweber/DATA/MC_DATA/BUSRA_23_10_22/KM1086/plots/plate/ashleys_plate_probabilities.pdf']
Building DAG of jobs...
Error: Directory cannot be locked. This usually means that another Snakemake instance is running on this directory. Another possibility is that a previous run exited unexpectedly.
Hashing conda environment https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/bwa/index/environment.yaml.
Hashing conda environment https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/fastqc/environment.yaml.
Hashing conda environment workflow/envs/ashleys_base.yaml.
Hashing conda environment workflow/envs/ashleys_rtools.yaml.
FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="cfe305df3c66d562f9326be7e9b42fc3d49c37623a7323f93ceb6de1e2c25ba8"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/bwa/index/environment.yaml
#   prefix: /conda-envs/5681728a49bd83ceed09ba194330c858
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bwa ==0.7.17
RUN mkdir -p /conda-envs/5681728a49bd83ceed09ba194330c858
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/bwa/index/environment.yaml /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/08d4368302a4bdf7eda6b536495efe7d
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - fastqc ==0.11.9
RUN mkdir -p /conda-envs/08d4368302a4bdf7eda6b536495efe7d
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/fastqc/environment.yaml /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml

# Conda environment:
#   source: workflow/envs/ashleys_base.yaml
#   prefix: /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc
#   name: ashleys_base
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - defaults
#   dependencies:
#     - samtools
#     - tabix
#     - bwa
#     - sambamba
#     - mosaicatcher
#     - alfred
#     - ashleys-qc
#     - pandas
#     # - pysam
RUN mkdir -p /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc
COPY workflow/envs/ashleys_base.yaml /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc/environment.yaml

# Conda environment:
#   source: workflow/envs/ashleys_rtools.yaml
#   prefix: /conda-envs/4cda6d03454db08ca24e6d039a2ce789
#   name: rtools
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#     - anaconda
#   dependencies:
#     # - bioconductor-biocparallel
#     # - bioconductor-bsgenome
#     # - bioconductor-bsgenome.hsapiens.ucsc.hg19
#     # - bioconductor-bsgenome.hsapiens.ucsc.hg38
#     # - bioconductor-fastseg
#     # - bioconductor-genomicalignments
#     - bioconductor-genomicranges
#     # - bioconductor-rsamtools
#     # - bioconductor-s4vectors
#     - r-assertthat
#     - r-base
#     # - r-biocmanager
#     - r-cowplot
#     - r-data.table
#     # - r-devtools
#     # - r-doparallel
#     # - r-foreach
#     - r-ggplot2
#     # - r-gtools
#     - r-reshape2
#     # - r-zoo
#     # - r-dplyr
#     # - r-mc2d
#     # - r-pheatmap
#     # - bioconductor-complexheatmap
#     # - r-gplots
#     - r-scales
#     - r-rcolorbrewer
#     # - r-stringr
#     - r-cairo
#     - fonts-anaconda
#     # NEW
#     - bioconductor-edger
#     - r-r.utils
#     # PLATE PLOT
#     - r-dplyr
#     - r-platetools
#     - r-viridis
RUN mkdir -p /conda-envs/4cda6d03454db08ca24e6d039a2ce789
COPY workflow/envs/ashleys_rtools.yaml /conda-envs/4cda6d03454db08ca24e6d039a2ce789/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc --file /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc/environment.yaml && \
    mamba env create --prefix /conda-envs/4cda6d03454db08ca24e6d039a2ce789 --file /conda-envs/4cda6d03454db08ca24e6d039a2ce789/environment.yaml && \
    mamba clean --all -y
