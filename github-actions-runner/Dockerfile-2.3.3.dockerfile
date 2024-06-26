FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="b4bf346d2b78ca78f3eab94bacc2a487b682a40690bd07db57f6c50bb52d1ce0"

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
#   prefix: /conda-envs/eba7cf011bffb12bcf25d066f6955913
#   name: ashleys_base
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools
#     - tabix
#     - bwa
#     - sambamba
#     - mosaicatcher
#     # - alfred
#     - ashleys-qc
#     - pandas
#     # PUBLISHDIR
#     - rsync
#     # MULTIQC
#     - multiqc
#     # Fix sklearn update
#     - scikit-learn=1.2.2
#     # Fix ashleys GH issue
#     - python=3.10
RUN mkdir -p /conda-envs/eba7cf011bffb12bcf25d066f6955913
COPY workflow/envs/ashleys_base.yaml /conda-envs/eba7cf011bffb12bcf25d066f6955913/environment.yaml

# Conda environment:
#   source: workflow/envs/ashleys_rtools.yaml
#   prefix: /conda-envs/9b847fc31baae8e01dfb7ce438a56b71
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
#     # GC_correction
#     - r-tidyr
#     - r-ggpubr
#     # SOLVE R lib issue
#     - r-stringi=1.7.12
RUN mkdir -p /conda-envs/9b847fc31baae8e01dfb7ce438a56b71
COPY workflow/envs/ashleys_rtools.yaml /conda-envs/9b847fc31baae8e01dfb7ce438a56b71/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/eba7cf011bffb12bcf25d066f6955913 --file /conda-envs/eba7cf011bffb12bcf25d066f6955913/environment.yaml && \
    mamba env create --prefix /conda-envs/9b847fc31baae8e01dfb7ce438a56b71 --file /conda-envs/9b847fc31baae8e01dfb7ce438a56b71/environment.yaml && \
    mamba clean --all -y
