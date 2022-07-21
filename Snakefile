from snakemake.utils import min_version

min_version("7.4.0")

configfile: "config/config.yaml"

containerized: "docker://weber8thomas/ashleys-qc-pipeline:1.1"


include: "workflow/rules/common.smk"
include: "workflow/rules/rules.smk"
include: "workflow/rules/examples.smk"

if config["mosaicatcher_pipeline"] is False:
    if config["dl_external_files"] is True:
        rule all_ashleys:
            input:
                rules.dl_external_data.output,
                rules.dl_external_data_index.output
    else:
        rule all_ashleys:
            input:
                get_final_output()