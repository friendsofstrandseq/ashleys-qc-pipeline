
rule save_conda_versions_ashleys_base:
    output:
        "{folder}/{sample}/config/conda_export/ashleys_base.yaml",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        "conda env export > {output}"

rule save_conda_versions_ashleys_rtools:
    output:
        "{folder}/{sample}/config/conda_export/ashleys_rtools.yaml",
    conda:
        "../envs/ashleys_rtools.yaml"
    shell:
        "conda env export > {output}"

