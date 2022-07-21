FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="1a8e4ba1bc686a7e19b8885f826469bd28a11cfce58a623d88a3808dee686fd1"

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
#   source: workflow/envs/ashleys.yaml
#   prefix: /conda-envs/40a7bb20d75edfcabb61d7db444cc9ea
#   name: ashleys
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - ashleys-qc
RUN mkdir -p /conda-envs/40a7bb20d75edfcabb61d7db444cc9ea
COPY workflow/envs/ashleys.yaml /conda-envs/40a7bb20d75edfcabb61d7db444cc9ea/environment.yaml

# Conda environment:
#   source: workflow/envs/mc_bioinfo_tools.yaml
#   prefix: /conda-envs/638f610ec9ecb52e489f031fa8ac523b
#   name: mc-bioinfo-tools
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools
#     - tabix
#     - bwa
#     - sambamba
#     - mosaicatcher
RUN mkdir -p /conda-envs/638f610ec9ecb52e489f031fa8ac523b
COPY workflow/envs/mc_bioinfo_tools.yaml /conda-envs/638f610ec9ecb52e489f031fa8ac523b/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/40a7bb20d75edfcabb61d7db444cc9ea --file /conda-envs/40a7bb20d75edfcabb61d7db444cc9ea/environment.yaml && \
    mamba env create --prefix /conda-envs/638f610ec9ecb52e489f031fa8ac523b --file /conda-envs/638f610ec9ecb52e489f031fa8ac523b/environment.yaml && \
    mamba clean --all -y
