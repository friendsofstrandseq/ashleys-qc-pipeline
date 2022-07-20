from snakemake.utils import min_version

min_version("7.4.1")

configfile: "config/config.yaml"


include: "workflow/rules/common.smk"
include: "workflow/rules/rules.smk"

if config["mosaicatcher_pipeline"] is False:
    rule all:
        input:
            get_final_output()