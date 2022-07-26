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
