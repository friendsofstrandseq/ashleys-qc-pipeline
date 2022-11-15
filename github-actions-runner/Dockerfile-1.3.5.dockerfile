FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="8cf785feb91c64169420444ebf4d9987b9ba85206ebca0ebd2dd212f6e480c10"

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
#   prefix: /conda-envs/257f364c2fc06b9bef48012f8e00427c
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
RUN mkdir -p /conda-envs/257f364c2fc06b9bef48012f8e00427c
COPY workflow/envs/ashleys_rtools.yaml /conda-envs/257f364c2fc06b9bef48012f8e00427c/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc --file /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc/environment.yaml && \
    mamba env create --prefix /conda-envs/257f364c2fc06b9bef48012f8e00427c --file /conda-envs/257f364c2fc06b9bef48012f8e00427c/environment.yaml && \
    mamba clean --all -y
