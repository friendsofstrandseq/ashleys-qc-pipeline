![logo](docs/images/logo.png)

[![Workflow checks](https://github.com/friendsofstrandseq/ashleys-qc-pipeline/actions/workflows/main.yaml/badge.svg)](https://github.com/friendsofstrandseq/ashleys-qc-pipeline/actions/workflows/main.yaml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.4.0-brightgreen.svg)](https://snakemake.github.io)

---

Strand-Seq Quality Control Pipeline based on ashleys-qc

---

# Overview of this workflow

This workflow uses [Snakemake](https://github.com/snakemake/snakemake) perform Quality Control analysis on Strand-Seq single-cell
sequencing data. The starting point are single-cell FASTQ files from Strand-seq experiments and the final output produced is a folder with clean selected BAM files. The pipeline can identify automatically high-quality libraries through ML-based analysis tool [ashleys-qc](https://github.com/friendsofstrandseq/ashleys-qc). Thus, the workflow goes through the following steps:

1. FASTQ sequencing Quality Control through [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Mapping FASTQ against a reference genome throught [BWA](http://bio-bwa.sourceforge.net/)
3. Sorting, Deduplicating and Indexing of BAM files through [Samtools](http://www.htslib.org/) & [sambaba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
4. Generating features and use [ashleys-qc](https://github.com/friendsofstrandseq/ashleys-qc) model to identify high-quality cells

# Quick Start

0. [Optional] Install [Singularity](https://www.sylabs.io/guides/3.0/user-guide/)
1. To prevent conda channel errors

```bash
conda config --set channel_priority strict
```

2. Install snakemake through conda

```bash
conda create -n snakemake -c defaults -c anaconda -c conda-forge -c bioconda snakemake && conda activate snakemake
```

3. Clone the repository

```bash
git clone https://github.com/friendsofstrandseq/ashleys-qc-pipeline.git && cd ashleys-qc-pipeline
```

4. Run on example data on only one small chromosome (`<disk>` must be replaced by your disk letter/name, `/g` or `/scratch` at EMBL for example)

```bash
snakemake --cores 6 --configfile .tests/config/simple_config.yaml --profile workflow/snakemake_profiles/local/conda_singularity --singularity-args "-B /<disk>:/<disk>"
```

5. Run your own analysis **locally** (`<disk>` must be replaced by your disk letter/name, `/g` or `/scratch` at EMBL for example)

```bash
snakemake --cores 6 --config data_location=<PATH> --profile workflow/snakemake_profiles/local/conda_singularity --singularity-args "-B /<disk>:/<disk>" --latency-wait 60
```

---

**ℹ️ Note**

- Steps 1 - 3 are required only during first execution
- After the first execution, do not forget to go in the git repository and to activate the snakemake environment

---

# Parameters

## MosaiCatcher arguments

---

**ℹ️ Note**

All these arguments can be specified in two ways:

1. In the config/config.yaml file, by replacing existing values
2. Using the `--config` snakemake argument (`--config` must be called only one time with all the arguments behind it, e.g: `--config data_location=<INPUT>`)

---

### Input/output options

| Parameter       | Comment                                  | Parameter type | Default            |
| --------------- | ---------------------------------------- | -------------- | ------------------ |
| `data_location` | Path to parent folder containing samples | String         | .tests/data_CHR17/ |
| `email`         | Email address for completion summary     | String         | None               |

### Boolean parameters

| Parameter           | Comment                                                                                                                                                        | Default | Experimental |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------- | ------------ |
| `dl_external_files` | Allow to retrieve automatically external files (GRCh38 reference genome + 1000G SNV VCF file) required to run the pipeline.                                    | False   |              |
| `hand_selection`    | Allow to identify manually high-quality strand-seq libraries.                                                                                                  | False   | X            |
| `GC_analysis`       | Enable/Disable GC analysis and correction of Strand-Seq libraries libraries.                                                                                   | False   | X            |
| `plate_orientation` | If GC_analysis enabled and conditions tested by rows/columns, set the orientation (landscape/portrait) to perform a row/column-wise analysis of the libraries. | False   | X            |

### External files

| Parameter   | Comment          | Default | Other possibilities |
| ----------- | ---------------- | ------- | ------------------- |
| `reference` | Reference genome | hg38    | hg19, T2T           |

### Experimental: hand-selection related parameters

| Parameter     | Comment                                                                                             | Default       |
| ------------- | --------------------------------------------------------------------------------------------------- | ------------- |
| `window`      | Window size used for binning by mosaic count (Can be of high importance regarding library coverage) | 200000        |
| `chromosomes` | List of chromosomes to be processed in the pipeline                                                 | chr1..22,chrX |

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

### FastQC files

---

File path: `<FOLDER>/<SAMPLE>/fastqc/<CELL>.[1|2].[html|zip]`

---

You will be able to find all QC analysis of raw FastQ files in the path above.

![fastqc](docs/images/fastqc.png)

### ashleys-qc prediction

---

File path: `<FOLDER>/<SAMPLE>/predictions/predictions.tsv`

---

Ashleys predictions can be found at the path above.

Predictions table are composed of 3 columns: the cell name, the binary prediction (1: Selected, 0: Unselected) and the associated probability of the SVC model (normal model with a binary cutoff of 0.5).

| cell                         | prediction | probability |
| ---------------------------- | ---------- | ----------- |
| BM510x3PE20405.sort.mdup.bam | 0          | 0           |
| BM510x3PE20413.sort.mdup.bam | 0          | 0.43        |
| BM510x3PE20409.sort.mdup.bam | 0          | 0.32        |
| BM510x3PE20414.sort.mdup.bam | 1          | 0.85        |
| BM510x3PE20412.sort.mdup.bam | 0          | 0.29        |
| BM510x3PE20418.sort.mdup.bam | 1          | 0.84        |
| BM510x3PE20401.sort.mdup.bam | 1          | 0.87        |
| BM510x3PE20408.sort.mdup.bam | 1          | 0.88        |
| BM510x3PE20410.sort.mdup.bam | 1          | 0.93        |
| BM510x3PE20407.sort.mdup.bam | 1          | 0.89        |
| BM510x3PE20402.sort.mdup.bam | 1          | 0.95        |
| BM510x3PE20411.sort.mdup.bam | 1          | 0.91        |
| BM510x3PE20404.sort.mdup.bam | 1          | 0.89        |
| BM510x3PE20416.sort.mdup.bam | 1          | 0.84        |
| BM510x3PE20406.sort.mdup.bam | 1          | 0.91        |
| BM510x3PE20417.sort.mdup.bam | 1          | 0.93        |
| BM510x3PE20403.sort.mdup.bam | 1          | 0.91        |
| BM510x3PE20415.sort.mdup.bam | 1          | 0.9         |
| BM510x3PE20419.sort.mdup.bam | 1          | 0.91        |

### BAM selected folder

Selected libraries BAM files can be retrieved at the path above and can be used as an input for the [mosaicatcher-pipeline](https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git).

## Roadmap

### Major features

- [x] Jupyter Notebook hand selection of cells
- [x] HTML report
- [x] Multiple FASTA reference

### Minor features

- [x] replace `input_bam_location` by `data_location` (harmonization with [mosaicatcher-pipeline](https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git))

### Experimental feature: hand-selection of cells via Jupyter notebook

---

**ℹ️ Note**

Singularity execution (`--use-singularity`) is not available for this mode of execution at the moment due to path system issues.

---

If you wish to identify yourself the cells that seem uncorrect according to your expertise, you can use the experimental interactive Jupyter Notebook by passing to the config argument `hand_selection=True`.
By enabling this feature, the pipeline will run :

- [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher) count binning-based function
- plot Strand-Seq karyotype figures
- fire a Jupyter Notebook for analysis.

---

**⚠️ Warning**

If you are running the pipeline remotely and not on your local computer, you need first to open a [SSH tunnel (with Local Forwarding to port 5500)](https://www.ssh.com/academy/ssh/tunneling/example#local-forwarding) in order to access the Jupyter Notebook webpage.

---

The following command, including snakemake `--notebook-listen` (allow to chose the port oppened: arbitrary selected to 5500) and `--edit-notebook` (fire a jupyter server that allow graphical interaction instead of directly run the complete notebook) arguments, need to be passed to test this feature:

```bash
snakemake --cores 12 --profile workflow/snakemake_profiles/local/conda --config hand_selection=True data_location=<INPUT> \
  --notebook-listen localhost:5500 --edit-notebook <INPUT>/<SAMPLE>/cell_selection/labels_raw.tsv
```

Then, you can accessing Jupyter Notebook with your favorite web browser through the following URL:

```
http://localhost:5500
```

You will need to open the following untitled with the following pattern : `tmp[XXX].hand_selection.py.ipynb` and follow the instructions inside it.

Instructions listed in the notebook are also listed here:

---

#### _Jupyter Notebook instructions_

1. Run first the top cell (enable to use snakemake arguments)
2. Follow the different cells
   1. Symlink PDF plots to jupyter nb directory
   2. Enable jupyter widgets
   3. Display Mosaic Count PDF inside Jupyter Notebook
   4. Display graphical table to allow you to unselect low-quality libraries (that will be not processed in the remaining analysis)
      _Please <ins>do not</ins> RE-select cells that were automatically unselected cells (corresponding to low-coverage by [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher) count program and not possible to process by the pipeline)_
   5. Export to pandas dataframe
   6. Save & clean data
3. File > Close and Halt
4. Click on Quit (top right)
5. Close the webpage

---

![nb0](docs/images/nb_0.png)

![nb1](docs/images/nb_1.png)

![nb2](docs/images/nb_2.png)

However, as the previous command point to a specific output file related to the Jupyter Notebook snakemake rule, the snakemake command need to be runned again as the following to complete its execution (current snakemake limitation (7.9.0)).

Thus, after closing the notebook, run the following commands:

```bash
# Fix snakemake issue regarding metadata & incomplete jobs (to remove when solved)
rm .snakemake/incomplete/
```

```bash
# Snakemake > 7.8 changes its rerun behavior. Before, rerunning jobs relied purely on file modification times. https://github.com/snakemake/snakemake/issues/1694
# Run snakemake touch function to prevent timestamps errors
snakemake --cores 12 --use-conda --config hand_selection=True data_location=<INPUT> --touch
```

```bash
snakemake --cores 12 --use-conda --config hand_selection=True data_location=<INPUT>
```

# Authors (alphabetical order)

## Contributors

- Ebert Peter
- Grimes Karen
- Gros Christina
- Korbel Jan
- Marschall Tobias
- Sanders Ashley
- Weber Thomas (maintainer and current developer)

# References

> Gros, Christina, Ashley D Sanders, Jan O Korbel, Tobias Marschall, and Peter Ebert. “ASHLEYS: Automated Quality Control for Single-Cell Strand-Seq Data.” Bioinformatics 37, no. 19 (October 1, 2021): 3356–57. https://doi.org/10.1093/bioinformatics/btab221.
