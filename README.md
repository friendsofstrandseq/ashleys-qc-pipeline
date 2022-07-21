Strand-Seq Quality Control Pipeline based on ashleys-qc

#  Overview of this workflow

This workflow uses [Snakemake](https://github.com/snakemake/snakemake) perform Quality Control analysis on Strand-Seq single-cell
sequencing data. The starting point are single-cell FASTQ files from Strand-seq experiments and the final output produced is a folder with 
clean selected BAM files. Based on the user wish, the pipeline can (A) identify automatically high-quality libraries through ML-based analysis tool [ashleys-qc](https://github.com/friendsofstrandseq/ashleys-qc) or (B) wait for the user to realise a hand-selection of correct libraries through a Jupyter Notebook. Thus, the workflow goes through the following steps:

  1. FASTQ sequencing Quality Control through [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  2. Mapping FASTQ against a reference genome throught [BWA](http://bio-bwa.sourceforge.net/)
  3. Sorting, Deduplicating and Indexing of BAM files through [Samtools](http://www.htslib.org/) & [sambaba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)

At this point the pipeline can automatically continue through either:

  - 4A. Generating features and use ashleys-qc model to identify high-quality cells

  - 4B. Run [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher), plot figures and launch jupyter notebook for analysis

## Setup
```
Run the pipeline with
``` python
snakemake --use-conda --cores n
```

## Parameters

### TODO 

## Output

### TODO


## Roadmap

- [ ] HTML report
- [ ] Zenodo FASTA + index files