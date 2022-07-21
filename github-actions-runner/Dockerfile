FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="240ec4738f90572c25b4456fce057e148ac0c5238337781b1a9d4b39ced9aecb"

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

# Conda environment:
#   source: workflow/envs/notebook.yaml
#   prefix: /conda-envs/9df293fe70bdfb960cde51bd8c039072
#   name: notebook
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - ipysheet
#     - jupyter
#     - ipywidgets
#     - pandas
#     - nbconvert==5.4.1
RUN mkdir -p /conda-envs/9df293fe70bdfb960cde51bd8c039072
COPY workflow/envs/notebook.yaml /conda-envs/9df293fe70bdfb960cde51bd8c039072/environment.yaml

# Conda environment:
#   source: workflow/envs/rtools.yaml
#   prefix: /conda-envs/fb2e0991658d68d49ab06d2bd2375acf
#   name: rtools
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#     - anaconda
#   dependencies:
#     - bioconductor-biocparallel=1.16.6
#     - bioconductor-bsgenome
#     - bioconductor-bsgenome.hsapiens.ucsc.hg19
#     - bioconductor-bsgenome.hsapiens.ucsc.hg38=1.4.1
#     - bioconductor-fastseg=1.28.0
#     - bioconductor-genomicalignments=1.18.1
#     - bioconductor-genomicranges=1.34.0
#     - bioconductor-rsamtools=1.34.0
#     - bioconductor-s4vectors=0.20.1
#     - r-assertthat=0.2.1
#     - r-base=3.5.1
#     - r-biocmanager
#     - r-cowplot=1.0.0
#     - r-data.table=1.12.6
#     - r-devtools=2.2.2
#     - r-doparallel
#     - r-foreach
#     - r-ggplot2=3.3.0
#     - r-gtools=3.8.1
#     - r-reshape2=1.4.3
#     - r-scales=1.1.0
#     - r-zoo=1.8_3
#     - r-dplyr=0.8.5
#     - r-mc2d=0.1_18
#     - r-pheatmap=1.0.12
#     - bioconductor-complexheatmap=2.0.0
#     - r-gplots=3.0.3
#     - r-scales=1.1.0
#     - r-rcolorbrewer=1.1_2
#     - r-stringr=1.4.0
#     - r-cairo
#     - fonts-anaconda
RUN mkdir -p /conda-envs/fb2e0991658d68d49ab06d2bd2375acf
COPY workflow/envs/rtools.yaml /conda-envs/fb2e0991658d68d49ab06d2bd2375acf/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/40a7bb20d75edfcabb61d7db444cc9ea --file /conda-envs/40a7bb20d75edfcabb61d7db444cc9ea/environment.yaml && \
    mamba env create --prefix /conda-envs/638f610ec9ecb52e489f031fa8ac523b --file /conda-envs/638f610ec9ecb52e489f031fa8ac523b/environment.yaml && \
    mamba env create --prefix /conda-envs/9df293fe70bdfb960cde51bd8c039072 --file /conda-envs/9df293fe70bdfb960cde51bd8c039072/environment.yaml && \
    mamba env create --prefix /conda-envs/fb2e0991658d68d49ab06d2bd2375acf --file /conda-envs/fb2e0991658d68d49ab06d2bd2375acf/environment.yaml && \
    mamba clean --all -y
