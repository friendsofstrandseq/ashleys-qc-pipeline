## Rules dedicated to download index Fasta files
## ---------------------------------------------------------------

import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


rule download_hg19_reference:
    input:
        HTTP.remote(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/hg19.fa",
    log:
        "workflow/data/ref_genomes/log/hg19.ok",
    run:
        directory = "workflow/data/ref_genomes/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        shell("mv {input} workflow/data/ref_genomes/hg19.fa.gz")
        shell("gunzip workflow/data/ref_genomes/hg19.fa.gz")


rule download_hg38_reference:
    input:
        HTTP.remote(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/hg38.fa",
    log:
        "workflow/data/ref_genomes/log/hg38.ok",
    run:
        directory = "workflow/data/ref_genomes/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        shell("mv {input} workflow/data/ref_genomes/hg38.fa.gz")
        shell("gunzip workflow/data/ref_genomes/hg38.fa.gz")


rule download_T2T_reference:
    input:
        HTTP.remote(
            "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/T2T.fa",
    log:
        "workflow/data/ref_genomes/log/T2T.ok",
    run:
        directory = "workflow/data/ref_genomes/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        shell("mv {input} workflow/data/ref_genomes/T2T.fa.gz")
        shell("gunzip workflow/data/ref_genomes/T2T.fa.gz")


rule samtools_faindex:
    input:
        ancient("{file}.fa"),
    output:
        "{file}.fa.fai",
    log:
        "{file}.log",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        "samtools faidx {input}"
