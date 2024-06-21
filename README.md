![logo](docs/images/logo.png)

[![Workflow checks](https://github.com/friendsofstrandseq/ashleys-qc-pipeline/actions/workflows/main_test.yaml/badge.svg)](https://github.com/friendsofstrandseq/ashleys-qc-pipeline/actions/workflows/main_test.yaml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.12.0-brightgreen.svg)](https://snakemake.github.io)

# ashleys-qc-pipeline

ashleys-qc-pipeline performs Quality Control analysis on Strand-Seq single-cell
sequencing data. The starting point are single-cell FASTQ files from Strand-seq experiments and the final output produced is a folder with clean selected BAM files. The pipeline can identify automatically high-quality libraries through ML-based analysis tool [ashleys-qc](https://github.com/friendsofstrandseq/ashleys-qc). Thus, the workflow goes through the following steps:

1. FASTQ sequencing Quality Control through [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Mapping FASTQ against a reference genome throught [BWA](http://bio-bwa.sourceforge.net/)
3. Sorting, Deduplicating and Indexing of BAM files through [Samtools](http://www.htslib.org/) & [sambaba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
4. Generating features and use [ashleys-qc](https://github.com/friendsofstrandseq/ashleys-qc) model to identify high-quality cells

## Documentation

**📚 Homepage:** [https://friendsofstrandseq.github.io/mosaicatcher-docs/](https://friendsofstrandseq.github.io/mosaicatcher-docs/)

## Authors (alphabetical order)

- Ebert Peter
- Grimes Karen
- Gros Christina
- Korbel Jan
- Marschall Tobias
- Sanders Ashley
- Weber Thomas (maintainer and current developer)

## References

> MosaiCatcher v2 publication: Weber Thomas, Marco Raffaele Cosenza, and Jan Korbel. 2023. ‘MosaiCatcher v2: A Single-Cell Structural Variations Detection and Analysis Reference Framework Based on Strand-Seq’. Bioinformatics 39 (11): btad633. https://doi.org/10.1093/bioinformatics/btad633.

> Gros, Christina, Ashley D Sanders, Jan O Korbel, Tobias Marschall, and Peter Ebert. “ASHLEYS: Automated Quality Control for Single-Cell Strand-Seq Data.” Bioinformatics 37, no. 19 (October 1, 2021): 3356–57. https://doi.org/10.1093/bioinformatics/btab221.
