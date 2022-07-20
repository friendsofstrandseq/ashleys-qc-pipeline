# Pipeline for ashleys-qc
Automated installation of ASHLEYS
Feature generation for new data and prediction of the quality of each cell

## Setup
Clone the repository via
``` python
git clone https://github.com/friendsofstrandseq/ashleys-qc-pipeline.git
cd ashleys-qc-pipeline
```
Then create and activate the conda environment:
``` python
conda env create -f environment/pipeline_env.yml
conda activate ashleys-pipeline
```
Run the pipeline with
 ``` python
snakemake --use-conda --cores n
```

Specify file directories in config_snakemake.yaml

# TODO 
- [ ] link predictions/selection to mosaicatcher
- [ ] handle samples & config in a cleaner way
