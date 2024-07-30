FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="93e821ac2944264216e1176964afe6c8a6ee61b9daa831c249f6912f9ebb267f"

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
#   prefix: /conda-envs/32c736a65a401b33605acfa7a0241299
#   name: ashleys_base
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - ashleys-qc
#     - bwa
#     - mosaicatcher
#     - multiqc
#     - pandas
#     - python=3.10
#     - pysam
#     - rsync
#     - sambamba
#     - samtools
#     - scikit-learn=1.2.2
#     - tabix
RUN mkdir -p /conda-envs/32c736a65a401b33605acfa7a0241299
COPY workflow/envs/ashleys_base.yaml /conda-envs/32c736a65a401b33605acfa7a0241299/environment.yaml

# Conda environment:
#   source: workflow/envs/ashleys_rtools.yaml
#   prefix: /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0
#   name: rtools
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bioconductor-edger
#     - bioconductor-genomicranges
#     - fonts-conda-forge
#     - r-assertthat
#     - r-base
#     - r-cairo
#     - r-cowplot
#     - r-data.table
#     - r-dplyr
#     - r-ggplot2
#     - r-ggpubr
#     - r-platetools
#     - r-r.utils
#     - r-rcolorbrewer
#     - r-reshape2
#     - r-scales
#     - r-stringi=1.7.12
#     - r-tidyr
#     - r-viridis
RUN mkdir -p /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0
COPY workflow/envs/ashleys_rtools.yaml /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/32c736a65a401b33605acfa7a0241299 --file /conda-envs/32c736a65a401b33605acfa7a0241299/environment.yaml && \
    mamba env create --prefix /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0 --file /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0/environment.yaml && \
    mamba clean --all -y
