
rule generate_exclude_file_for_mosaic_count:
    input:
        bam=lambda wc: expand(
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[str(wc.sample)],
        ),
    output:
        excl="{folder}/{sample}/config/chroms_to_exclude.txt",
    log:
        "{folder}/log/config/{sample}/exclude_file.log",
    conda:
        "../envs/mc_base.yaml"
    params:
        chroms=config["chromosomes"],
    script:
        "../scripts/utils/generate_exclude_file.py"


rule mosaic_count:
    input:
        bam=lambda wc: expand(
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[str(wc.sample)],
        ),
        bai=lambda wc: expand(
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[str(wc.sample)],
        ),
        excl="{folder}/{sample}/config/chroms_to_exclude.txt",
    output:
        counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
        info="{folder}/{sample}/counts/{sample}.info_raw",
    log:
        "{folder}/log/counts/{sample}/mosaic_count.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    params:
        window=config["window"],
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        mosaicatcher count \
            --verbose \
            --do-not-blacklist-hmm \
            -o {output.counts} \
            -i {output.info} \
            -x {input.excl} \
            -w {params.window} \
            {input.bam} \
        > {log} 2>&1
        """

if (
    (config["window"] in [50000, 100000, 200000])
    and (config["reference"] == "hg38")
    and (config["normalized_counts"] is True)
):

    rule merge_blacklist_bins:
        input:
            norm="workflow/data/normalization/HGSVC.{window}.txt",
            whitelist="workflow/data/normalization/inversion-whitelist.tsv",
        output:
            merged="{folder}/{sample}/normalizations/HGSVC.{window}.merged.tsv",
        log:
            "{folder}/log/normalizations/{sample}/HGSVC.{window}.merged.tsv"
        conda:
            "../envs/mc_base.yaml"
        shell:  
            """
            workflow/scripts/normalization/merge-blacklist.py --merge_distance 500000 {input.norm} --whitelist {input.whitelist} --min_whitelist_interval_size 100000 > {output.merged} 2>> {log}
            """

    rule normalize_counts:
        input:
            counts="{folder}/{sample}/counts/{sample}.txt.filter.gz",
            norm=expand(
                "{folder}/{sample}/normalizations/HGSVC.{window}.merged.tsv",
                folder=config["data_location"],
                sample=samples,
                window=config["window"],
            ),
        output:
            # "{folder}/{sample}/counts/{sample}.txt.norm.gz",
            "{folder}/{sample}/counts/{window}.txt.gz",
        log:
            "{folder}/log/normalize_counts/{sample}_{window}.log",
        conda:
            "../envs/rtools.yaml"
        shell:
            """
            Rscript workflow/scripts/normalization/normalize.R {input.counts} {input.norm} {output} 2>&1 > {log}
            """

else:

    rule cp_mosaic_count:
        input:
            "{folder}/{sample}/counts/{sample}.txt.filter.gz",
        output:
            "{folder}/{sample}/counts/{sample}.txt.gz",
        log:
            "{folder}/log/counts/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        shell:
            "cp {input} {output}"

rule order_mosaic_count_output:
    input:
        "{folder}/{sample}/counts/{sample}.all.txt.fixme.gz",
    output:
        "{folder}/{sample}/counts/{sample}.all.txt.gz",
    log:
        "{folder}/log/counts/{sample}/{sample}.log",
    run:
        df = pd.read_csv(input[0], compression="gzip", sep="\t")
        df = df.sort_values(by=["sample", "cell", "chrom", "start"])
        df.to_csv(output[0], index=False, compression="gzip", sep="\t")


rule plot_mosaic_counts:
    input:
        counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
        info="{folder}/{sample}/counts/{sample}.info_raw",
    output:
        "{folder}/{sample}/plots/counts/CountComplete.classic.pdf",
    log:
        "{folder}/log/plot_mosaic_counts/{sample}.log",
    conda:
        "../envs/rtools.yaml"
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
        """