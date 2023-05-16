
rule fastqc:
    input:
        "{folder}/{sample}/fastq/{cell}.{pair}.fastq.gz",
    output:
        html=report(
            "{folder}/{sample}/fastqc/{cell}_{pair}_fastqc.html",
            category="FastQC",
            subcategory="{sample}",
            labels={"Sample": "{sample}", "Cell": "{cell}", "Pair": "{pair}"},
        ),
        zip="{folder}/{sample}/fastqc/{cell}_{pair}_fastqc.zip",
    params:
        "--quiet",
    log:
        "{folder}/log/fastqc/{sample}/{cell}_{pair}.log",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
    wrapper:
        "v1.7.0/bio/fastqc"


rule fastqc_aggregate:
    input:
        lambda wc: expand(
            "{folder}/{sample}/fastqc/{cell}_{pair}_fastqc.html",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
            pair=[1, 2],
        ),
    output:
        touch("{folder}/{sample}/fastqc/config/fastqc_output_touch.ok"),


rule samtools_idxstats:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/samtools_idxstats/{cell}.txt",
    log:
        "{folder}/{sample}/log/samtools_idxstats/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        "samtools idxstats {input} > {output}"


rule samtools_idxstats_aggr:
    input:
        lambda wc: expand(
            "{folder}/{sample}/samtools_idxstats/{cell}.txt",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
        ),
    output:
        touch("{folder}/{sample}/samtools_idxstats/config/samtools_idxstats_aggr_touch.ok"),



rule samtools_flagstats:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/samtools_flagstats/{cell}.txt",
    log:
        "{folder}/{sample}/log/samtools_flagstats/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        "samtools flagstats {input} > {output}"


rule samtools_flagstats_aggr:
    input:
        lambda wc: expand(
            "{folder}/{sample}/samtools_flagstats/{cell}.txt",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
        ),
    output:
        touch("{folder}/{sample}/samtools_flagstats/config/samtools_flagstats_aggr_touch.ok"),

rule samtools_stats:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/samtools_stats/{cell}.txt",
    log:
        "{folder}/{sample}/log/samtools_stats/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        "samtools stats {input} > {output}"


rule samtools_stats_aggr:
    input:
        lambda wc: expand(
            "{folder}/{sample}/samtools_stats/{cell}.txt",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
        ),
    output:
        touch("{folder}/{sample}/samtools_stats/config/samtools_stats_aggr_touch.ok"),

rule multiqc:
    input:
        fastqc="{folder}/{sample}/fastqc/config/fastqc_output_touch.ok",
        samtools_idxstats="{folder}/{sample}/samtools_idxstats/config/samtools_idxstats_aggr_touch.ok",
        samtools_stats="{folder}/{sample}/samtools_stats/config/samtools_stats_aggr_touch.ok",
        samtools_flagstats="{folder}/{sample}/samtools_flagstats/config/samtools_flagstats_aggr_touch.ok",
    output:
        report="{folder}/{sample}/multiqc/multiqc_report.html",
        outdir=directory("{folder}/{sample}/multiqc")
    params: 
        dirs = lambda wc, input: ["/".join(e.split("/")[:-2]) for e in list(input)]
    conda:
        "multiqc_megaqc"
    shell:
        "multiqc {params.dirs} --outdir {output.outdir}"