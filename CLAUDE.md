# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is ashleys-qc-pipeline, a Snakemake-based bioinformatics pipeline for Quality Control analysis on Strand-Seq single-cell sequencing data. The pipeline processes single-cell FASTQ files through mapping, QC analysis, and automated cell selection using the ML-based ashleys-qc tool.

**Key Documentation**: https://friendsofstrandseq.github.io/mosaicatcher-docs/

## Common Commands

### Running the Pipeline
```bash
# Basic pipeline execution
snakemake --cores 1 --use-conda --configfile config/config.yaml

# With specific data location
snakemake --cores 1 --use-conda --config data_location=.tests/data_CHR17

# List available commands
snakemake --cores 1 --config list_commands=True --verbose --debug

# Dry run with detailed output
snakemake --cores 1 --use-conda --configfile config/config.yaml -n -p

# Using conda with mamba frontend
snakemake --cores 1 --use-conda --conda-frontend mamba --configfile config/config.yaml
```

### Development and Testing
```bash
# Linting (from GitHub Actions)
snakemake --config data_location=.tests/data_CHR17 --lint

# Run with test data
snakemake --cores 1 --use-conda --configfile .tests/config/simple_config_ashleys.yaml --conda-frontend mamba -p --verbose --debug

# Format code using snakefmt (as shown in CI)
snakefmt workflow/

# Pre-commit hooks (install first: pip install pre-commit)
pre-commit install
pre-commit run --all-files  # Run on all files
pre-commit run              # Run on staged files only
```

### HPC Execution
```bash
# SLURM execution example
snakemake --cores 1 --use-conda --profile workflow/snakemake_profiles/HPC/slurm_generic/
```

## Architecture and Structure

### Core Components
- **Snakefile**: Main workflow entry point (`workflow/Snakefile`)
- **Rules**: Modular workflow components in `workflow/rules/`
  - `common.smk`: Common functions and setup
  - `rules.smk`: Main processing rules
  - `count.smk`: Read counting and binning
  - `gc.smk`: GC correction
  - `multiqc.smk`: Quality control reporting
  - `external_data.smk`: External data handling
  - `aggregate_fct.smk`: Aggregation functions

### Key Directories
- `workflow/`: Contains the Snakemake workflow files
- `config/`: Configuration files (main config.yaml and metadata)
- `workflow/envs/`: Conda environment definitions
- `workflow/data/`: Reference data and normalization files
- `workflow/scripts/`: Python and R scripts for processing
- `.tests/`: Test data and configurations

### Configuration System
- Main config: `config/config.yaml`
- Metadata config: `config/config_metadata.yaml`
- Test configs: `.tests/config/`
- The pipeline auto-generates sample configurations based on input data structure

### Data Structure Requirements
The pipeline expects a specific directory structure:
```
data_location/
├── SAMPLE_NAME/
│   └── fastq/
│       ├── CELL_ID.1.fastq.gz
│       └── CELL_ID.2.fastq.gz (for paired-end)
```

### Docker Integration
- Pipeline can run in Docker containers: `docker://weber8thomas/ashleys-qc-pipeline:VERSION`
- Controlled by `mosaicatcher_pipeline` config flag
- Container versions tracked in `config/config.yaml`

### Key Processing Steps
1. **FastQC**: Quality control of raw FASTQ files
2. **BWA Mapping**: Alignment against reference genome
3. **Samtools/Sambamba**: BAM processing (sort, deduplicate, index)
4. **Ashleys-QC**: ML-based quality assessment and cell selection
5. **Counting**: Generate read counts in genomic bins
6. **GC Correction**: Normalize for GC content bias
7. **Plotting**: Generate QC visualizations

### Environment Management
- Uses conda environments defined in `workflow/envs/`
- Main environment: `ashleys_base.yaml`
- Specialized environments for R tools, plotting, etc.
- Python 3.10 is the standard version

### Testing Infrastructure
- Test data in `.tests/data_CHR17/`
- GitHub Actions CI/CD in `.github/workflows/`
- Supports both local and SLURM execution testing
- Linting with Super-Linter and snakefmt

### Integration Points
- Can be used as submodule in mosaicatcher-pipeline
- Controlled by `mosaicatcher_pipeline` config flag
- Supports GENECORE data structure integration
- Hand selection through Jupyter notebooks

## Development Notes

### Code Style and Quality
- Python scripts use standard conventions (formatted with ruff)
- R scripts follow bioconductor patterns
- Snakemake files use consistent indentation and structure (formatted with snakefmt)
- Pre-commit hooks ensure code quality:
  - **Basic hooks**: trailing whitespace, end-of-file fixes, YAML validation
  - **Ruff**: Python formatting and linting (replaces black + flake8)
  - **Snakefmt**: Snakemake file formatting
  - **Snakemake --lint**: Workflow validation

### Key Configuration Parameters
- `reference`: Reference genome (hg38, hg19, T2T, mm10, mm39)
- `window`: Bin size for counting (default: 200000)
- `ashleys_threshold`: ML classification threshold (default: 0.5)
- `paired_end`: Whether data is paired-end (default: True)
- `bypass_ashleys`: Skip automatic quality control
- `multistep_normalisation`: Advanced normalization options

### HPC Profile Usage
Multiple HPC profiles available in `workflow/snakemake_profiles/HPC/`:
- `slurm_generic/`: Generic SLURM configuration
- `slurm_EMBL/`: EMBL-specific settings
- `slurm_BIH/`: BIH-specific settings
- `lsf_generic/`: LSF support

### Troubleshooting
- Check log files in output directories
- Verify conda environment creation
- Ensure reference genome files are accessible
- Test with provided `.tests/` data first