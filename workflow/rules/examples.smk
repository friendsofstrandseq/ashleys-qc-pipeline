import os
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()


# TODO: Adapt according reference
rule dl_external_data:
    """
    rule fct: Download External files 
    input: files stored on Zenodo
    output: touch file to check if everything was running correctly
    """
    input:
        HTTP.remote(
            "https://sandbox.zenodo.org/record/1074721/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
            keep_local=True,
        ),
    output:
        touch("sandbox.zenodo.org/config/dl_external_data.ok"),
    log:
        touch("sandbox.zenodo.org/log/config_output/dl_external_data.log"),


# TODO: Adapt according reference
rule dl_external_data_index:
    """
    rule fct: Download External files 
    input: files stored on Zenodo
    output: touch file to check if everything was running correctly
    """
    input:
        HTTP.remote(
            "https://sandbox.zenodo.org/record/1074721/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
            keep_local=True,
        ),
    output:
        touch("sandbox.zenodo.org/config/dl_external_data_index.ok"),
    log:
        touch("sandbox.zenodo.org/log/config_output/dl_external_data_index.log"),

# rule download_hg19_reference:
#     input:
#         HTTP.remote(
#             "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz",
#             keep_local=True,
#         ),
#     output:
#         "workflow/data/ref_genomes/hg19.fa",
#     log:
#         "workflow/data/ref_genomes/log/hg19.ok",
#     run:
#         directory = "workflow/data/ref_genomes/"
#         if not os.path.exists(directory):
#             os.makedirs(directory)
#         shell("mv {input} workflow/data/ref_genomes/hg19.fa.gz")
#         shell("gunzip workflow/data/ref_genomes/hg19.fa.gz")


# rule download_hg38_reference:
#     input:
#         HTTP.remote(
#             "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
#             keep_local=True,
#         ),
#     output:
#         "workflow/data/ref_genomes/hg38.fa",
#     log:
#         "workflow/data/ref_genomes/log/hg38.ok",
#     run:
#         directory = "workflow/data/ref_genomes/"
#         if not os.path.exists(directory):
#             os.makedirs(directory)
#         shell("mv {input} workflow/data/ref_genomes/hg38.fa.gz")
#         shell("gunzip workflow/data/ref_genomes/hg38.fa.gz")


# rule download_T2T_reference:
#     input:
#         HTTP.remote(
#             "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
#             keep_local=True,
#         ),
#     output:
#         "workflow/data/ref_genomes/T2T.fa",
#     log:
#         "workflow/data/ref_genomes/log/T2T.ok",
#     run:
#         directory = "workflow/data/ref_genomes/"
#         if not os.path.exists(directory):
#             os.makedirs(directory)
#         shell("mv {input} workflow/data/ref_genomes/T2T.fa.gz")
#         shell("gunzip workflow/data/ref_genomes/T2T.fa.gz")
