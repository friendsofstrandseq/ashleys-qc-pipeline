[![Workflow checks](https://github.com/friendsofstrandseq/ashleys-qc-pipeline/actions/workflows/main.yaml/badge.svg)](https://github.com/friendsofstrandseq/ashleys-qc-pipeline/actions/workflows/main.yaml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.4.0-brightgreen.svg)](https://snakemake.github.io)


_____________
Strand-Seq Quality Control Pipeline based on ashleys-qc
_____________
#  Overview of this workflow

This workflow uses [Snakemake](https://github.com/snakemake/snakemake) perform Quality Control analysis on Strand-Seq single-cell
sequencing data. The starting point are single-cell FASTQ files from Strand-seq experiments and the final output produced is a folder with clean selected BAM files. The pipeline can identify automatically high-quality libraries through ML-based analysis tool [ashleys-qc](https://github.com/friendsofstrandseq/ashleys-qc). Thus, the workflow goes through the following steps:

  1. FASTQ sequencing Quality Control through [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  2. Mapping FASTQ against a reference genome throught [BWA](http://bio-bwa.sourceforge.net/)
  3. Sorting, Deduplicating and Indexing of BAM files through [Samtools](http://www.htslib.org/) & [sambaba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
  4. Generating features and use [ashleys-qc](https://github.com/friendsofstrandseq/ashleys-qc) model to identify high-quality cells

# Quick Start

1. Install [Singularity](https://www.sylabs.io/guides/3.0/user-guide/) 
2. To prevent conda channel errors
```
conda config --set channel_priority strict
```
3. Install snakemake through conda
```
conda create -n snakemake -c conda-forge -c bioconda "snakemake>=7.4.1" && conda activate snakemake
```
4. Clone the repository 
``` 
git clone https://github.com/friendsofstrandseq/ashleys-qc-pipeline.git && cd ashleys-qc-pipeline
```
5. Download reference data for running your own analysis
```s
snakemake --cores 1 --config dl_external_files=True
```
6. Run on example data on only one small chromosome (`<disk>` must be replaced by your disk letter/name, `/g` or `/scratch` at EMBL for example)
```
snakemake --cores 6 --config input_bam_location=<PATH>/ashleys-qc-pipeline/.tests/data_CHR21 --use-conda --use-singularity --singularity-args "-B /<disk>:/<disk>" --latency-wait 60 
```
7. Run your own analysis (`<disk>` must be replaced by your disk letter/name, `/g` or `/scratch` at EMBL for example)
```
snakemake --cores 6 --config input_bam_location=<PATH> --use-conda --use-singularity --singularity-args "-B /<disk>:/<disk>" --latency-wait 60 
```
## Parameters


# Parameters

## MosaiCatcher arguments
________

**ℹ️ Note**
  
All these arguments can be specified in two ways:
1. In the config/config.yaml file, by replacing existing values
2. Using the `--config` snakemake argument (`--config` must be called only one time with all the arguments behind it, e.g: `--config input_bam_location=<INPUT>`)

________

### Input/output options

| Parameter            | Comment                                  | Parameter type | Default            |
| -------------------- | ---------------------------------------- | -------------- | ------------------ |
| `input_bam_location` | Path to parent folder containing samples | String         | .tests/data_CHR21/ |
| `email`              | Email address for completion summary     | String         | None               |

### Boolean parameters

| Parameter           | Comment                                                                                                                     | Default | Experimental |
| ------------------- | --------------------------------------------------------------------------------------------------------------------------- | ------- | ------------ |
| `dl_external_files` | Allow to retrieve automatically external files (GRCh38 reference genome + 1000G SNV VCF file) required to run the pipeline. | False   |              |
| `hand_selection`    | Allow to identify manually high-quality strand-seq libraries.                                                               | False   | X            |


### External files

| Parameter   | Comment          | Default                                                                                 |
| ----------- | ---------------- | --------------------------------------------------------------------------------------- |
| `reference` | Reference genome | sandbox.zenodo.org/record/1074721/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna |


## Snakemake arguments

Here are presented some essential snakemake options that could help you. 

```
--cores, -c
```
Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores. In case of cluster/cloud execution, this argument sets the number of total cores used over all jobs (made available to rules via workflow.cores).

```
--printshellcmds, -p
```
Recommended to print out the shell commands that will be executed.

```
--use-conda
```
If defined in the rule, run job in a conda environment. If this flag is not set, the conda directive is ignored and use the current environment (and path system) to execute the command.

```
--conda-frontend [mamba|conda] 
```
Choose the conda frontend for installing environments. Mamba is much faster and highly recommended but could not be installed by default on your system. Default: “conda”


```
--use-singularity 
```
If defined in the rule, run job within a singularity container. If this flag is not set, the singularity directive is ignored.

```
--singularity-args "-B /mounting_point:/mounting_point"
```
Pass additional args to singularity. `-B` stands for binding point between the host and the container.

```
--dryrun, -n 
```
Do not execute anything, and display what would be done. If you have a very large workflow, use –dry-run –quiet to just print a summary of the DAG of jobs.

```
--rerun-incomplete, --ri
```
Re-run all jobs the output of which is recognized as incomplete.


```
--keep-going, -k
```
Go on with independent jobs if a job fails.

```
-T, --retries, --restart-times
```
Number of times to restart failing jobs (defaults to 0).

```
--forceall, -F
```
Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.

---
**ℹ️ Note**

Currently, the binding command needs to correspond to the mounting point of your system (i.e: "/tmp:/tmp").
On seneca for example (EMBL), use `"/g:/g"` if you are working on `/g/korbel[2]` or `"/scratch:/scratch"` if you plan to work on `scratch`.

---

Obviously, all other [snakemake CLI options](https://snakemake.readthedocs.io/en/stable/executing/cli.html) can also be used. 


## Output

### TODO


## Roadmap

- [ ] HTML report
- [ ] Jupyter Notebook hand selection of cells
- [ ] Zenodo FASTA + index files

### Experimental feature: hand-selection of cells (Jupyter notebook)

If you wish to identify yourself the cells that seem uncorrect according to your expertise, you can use the experimental interactive Jupyter Notebook by passing to the config argument `hand_selection=True`.
By enabling this feature, the pipeline will run first [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher) count binning-based function, plot Strand-Seq karyotype figures, and then create a Jupyter Notebook for analysis.

---
**⚠️ Warning**

If you are running the pipeline remotely and not on your local computer, you need first to open a [SSH tunnel (with Local Forwarding)](https://www.ssh.com/academy/ssh/tunneling/example#local-forwarding) in order to access the Jupyter Notebook webpage. 

---

The following command (that comprise snakemake `--notebook-listen` and `--edit-notebook` arguments, need to be passed to test this feature:

```
snakemake --cores 12 --use-conda --config hand_selection=True input_bam_location=<INPUT> \
  --notebook-listen localhost:5500 --edit-notebook <INPUT>/<SAMPLE>/predictions/predictions_raw.tsv
```

However, as the previous command point to a specific output file related to the Jupyter Notebook snakemake rule, the snakemake command need to be runned again as the following to complete its execution (current snakemake limitation (7.9.0)):
```
snakemake --cores 12 --use-conda --config hand_selection=True input_bam_location=<INPUT>
```

---
**ℹ️ Note**

Singularity execution is not available for this mode of execution at the moment.

---