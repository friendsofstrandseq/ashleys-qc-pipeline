## Rules dedicated to download remote files
## ---------------------------------------------------------------

import os
# from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

storage:
    provider="http"

HTTP = HTTPRemoteProvider()


rule download_hg19_reference:
    input:
        # HTTP.remote(
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/hg19.fa",
    log:
        "workflow/data/ref_genomes/log/hg19.ok",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/hg19.fa.gz
        gunzip workflow/data/ref_genomes/hg19.fa.gz
        """


rule download_hg38_reference:
    input:
        # HTTP.remote(
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/hg38.fa",
    log:
        "workflow/data/ref_genomes/log/hg38.ok",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/hg38.fa.gz
        gunzip workflow/data/ref_genomes/hg38.fa.gz
        """


rule download_T2T_reference:
    input:
        # HTTP.remote(
        storage.http(
            "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/T2T.fa",
    log:
        "workflow/data/ref_genomes/log/T2T.ok",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/T2T.fa.gz
        gunzip workflow/data/ref_genomes/T2T.fa.gz
        """


rule download_mm10_reference:
    input:
        # HTTP.remote(
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/mm10.fa",
    log:
        "workflow/data/ref_genomes/log/mm10.ok",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/mm10.fa.gz
        gunzip workflow/data/ref_genomes/mm10.fa.gz
        """


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
