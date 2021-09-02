# Pipeline for the automated quality control of Strand-seq libraries (ASHLEYS QC)
This pipeline automates (i) data (pre-) processing, i.e. alignment of Strand-seq libraries to a
genome reference, (ii) training an ASHLEYS classification model for labeling Strand-seq
libraries, (iii) using an ASHLEYS classification model to label new Strand-seq libraries.

Data (pre-) processing tasks are realized in a generic manner and should be applicable
to any type of Strand-seq data set (see below for hard input requirements).

Metadata preprocessing such as deriving library labels on the basis of manually curated
annotation tables is only possible for selected projects where the necessary transformations
have been implemented (e.g., as for the HGSVC project).

## Hard input requirements
- Strand-seq libraries come as gzipped FASTQ files
  - file extension can be configured, will be normalized
- file names are composed using only the characters A-Z, a-z, 0-9, "-", "_", "."
- each (library) file name contains the sample name
- only paired-end libraries are supported
  - the mate number (1 or 2) is the last component of the file name before the file extension
  - file names are identical per library except for the mate number
- the data input folder is organized "by sample" as exemplified below:

```
STRANDSEQ_DATA_ROOT
  |
  -- sampleA
  |     |
  |     -- sampleA_library1_1.fastq.gz
  |     -- sampleA_library1_2.fastq.gz
  |     -- sampleA_library2_1.fastq.gz
  |     -- sampleA_library2_2.fastq.gz
  |     -- ...
  |
  -- sampleB
  |     |
  |     -- sampleB_library1_1.fastq.gz
        -- sampleB_library1_2.fastq.gz
        -- ...
```

## Using the pipeline

### Environment setup
Clone the pipeline repository via

```bash
git clone https://github.com/friendsofstrandseq/ashleys-qc-pipeline.git
```

The repository includes a [Conda](https://docs.conda.io/en/latest/miniconda.html)
environment file containing all necessary packages to run the (Snakemake) pipeline.
Build the Conda environment (in this example: in the folder "run_ashqc") as follows:

```bash
conda env create -f ashleys-qc-pipeline/environment/conda_smk.yml -p ./run_ashqc
conda activate ./run_ashqc
```

### Configuring your pipeline run
Before running the pipeline, you need to configure...

**(i)** the location of the Strand-seq input and reference data on your system:  
Please refer to the file `smk_config/cfg_data_hhu_hgsvc.yml` to see an example
for a DATA CONFIG file.

**(ii)** the number of available CPU cores in your environment:  
Please refer to the files `smk_config/cfg_run_hhu.yml` and `smk_config/cfg_run_laptop.yml`
to see examples for RUN CONFIG files for HPC cluster and
testing (laptop) environments.

**(iii)** Snakemake's behavior on your system by creating a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles):  
Please refer to the files `smk_profile/hhu_hpc/config.yaml` and `smk_profile/laptop/config.yaml`
to see examples for HPC cluster and testing (laptop) environments.

### Executing your pipeline run

#### Preliminaries: execution environment
Snakemake pipelines can be executed on a local machine (e.g., a laptop or desktop), or
on compute clusters distributing the work over several servers. In the latter case,
it is helpful to use a cluster status script that enables Snakemake to communicate
with the job scheduling system for more precise monitoring of the job status. This
pipeline repository does not include such a job status script because these scripts
are inherently system-specific. If in doubt, contact your local cluster IT support
and ask them if they happen to provide such a script. You can find more info in the
[official Snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html?highlight=cluster%20status#using-cluster-status).

#### Stage 1: collect input files
Before a complete pipeline run is triggered, it is necessary to check if all
input samples and corresponding library files have been properly identified by
the pipeline. Collecting and linking the input files takes only a few seconds.
Run the pipeline as follows:

```bash
snakemake -d ./YOUR_PIPELINE_RUN_FOLDER \
    --configfiles smk_config/cfg_params.yml YOUR_DATA_CONFIG.yml YOUR_RUN_CONFIG.yml \
    --profile PATH_TO_YOUR_PROFILE \
    link_sseq_input
```

You can now check that all Strand-seq samples have been properly picked up
by the pipeline. Note that the pipeline automatically checks that (i) each sample
has a number N of single-cell libraries where N has to be an integer multiple
of the config parameter `libs_per_plate` (default: 96); (ii) the names of the FASTQ
files enable an unambiguous pairing of 1st and 2nd mate per library, i.e. what is
commonly indicated via the FASTQ suffix `_1` and `_2`. You can glance at the linked
input files under this path:

```
YOUR_PIPELINE_RUN_FOLDER/input/fastq/
```

#### Stage 2: trigger a complete pipeline run

**TODO** Incomplete - performs only checksum computation and alignments

```bash
snakemake -d ./YOUR_PIPELINE_RUN_FOLDER \
    --configfiles smk_config/cfg_params.yml YOUR_DATA_CONFIG.yml YOUR_RUN_CONFIG.yml \
    --profile PATH_TO_YOUR_PROFILE \
    [--cluster-status CLUSTER_STATUS_SCRIPT] \
    run_sseq_checksums \
    run_sseq_alignments
```

## Citation
If you are using this pipeline (or parts of it) in you own work, please cite the following paper:

```
@article{2021ashleys,
  title    = "{ASHLEYS}: automated quality control for single-cell Strand-seq
              data",
  author   = "Gros, Christina and Sanders, Ashley D and Korbel, Jan O and
              Marschall, Tobias and Ebert, Peter",
  abstract = "SUMMARY: Single-cell DNA template strand sequencing (Strand-seq)
              enables chromosome length haplotype phasing, construction of
              phased assemblies, mapping sister-chromatid exchange events and
              structural variant discovery. The initial quality control of
              potentially thousands of single-cell libraries is still done
              manually by domain experts. ASHLEYS automates this tedious task,
              delivers near-expert performance and labels even large data sets
              in seconds. AVAILABILITY AND IMPLEMENTATION:
              github.com/friendsofstrandseq/ashleys-qc, MIT license.
              SUPPLEMENTARY INFORMATION: Supplementary data are available at
              Bioinformatics online.",
  journal  = "Bioinformatics",
  month    =  apr,
  year     =  2021,
  language = "en",
  doi      = "10.1093/bioinformatics/btab221"
}
```
