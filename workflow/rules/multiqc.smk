
rule fastqc:
    input:
        "{folder}/{sample}/fastq/{cell}.{pair}.fastq.gz",
    output:
        html=report(
            "{folder}/{sample}/multiqc/fastqc/{cell}_{pair}_fastqc.html",
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
            "{folder}/{sample}/multiqc/fastqc/{cell}_{pair}_fastqc.html",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
            pair=[1, 2],
        ),
    output:
        touch("{folder}/{sample}/multiqc/fastqc/config/fastqc_output_touch.ok"),


rule samtools_idxstats:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/multiqc/samtools_idxstats/{cell}.txt",
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
            "{folder}/{sample}/multiqc/samtools_idxstats/{cell}.txt",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
        ),
    output:
        touch(
            "{folder}/{sample}/multiqc/samtools_idxstats/config/samtools_idxstats_aggr_touch.ok"
        ),


rule samtools_flagstats:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/multiqc/samtools_flagstats/{cell}.txt",
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
            "{folder}/{sample}/multiqc/samtools_flagstats/{cell}.txt",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
        ),
    output:
        touch(
            "{folder}/{sample}/multiqc/samtools_flagstats/config/samtools_flagstats_aggr_touch.ok"
        ),


rule samtools_stats:
    input:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    output:
        "{folder}/{sample}/multiqc/samtools_stats/{cell}.txt",
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
            "{folder}/{sample}/multiqc/samtools_stats/{cell}.txt",
            folder=config["data_location"],
            sample=wc.sample,
            cell=cell_per_sample[wc.sample],
        ),
    output:
        touch(
            "{folder}/{sample}/multiqc/samtools_stats/config/samtools_stats_aggr_touch.ok"
        ),


rule multiqc:
    input:
        fastqc="{folder}/{sample}/multiqc/fastqc/config/fastqc_output_touch.ok",
        samtools_idxstats="{folder}/{sample}/multiqc/samtools_idxstats/config/samtools_idxstats_aggr_touch.ok",
        samtools_stats="{folder}/{sample}/multiqc/samtools_stats/config/samtools_stats_aggr_touch.ok",
        samtools_flagstats="{folder}/{sample}/multiqc/samtools_flagstats/config/samtools_flagstats_aggr_touch.ok",
        multiqc_input="{folder}/{sample}/multiqc/",
    output:
        report="{folder}/{sample}/multiqc/multiqc_report.html",
        outdir=report(
            directory("{folder}/{sample}/multiqc"),
            htmlindex="multiqc_report.html",
            category="MultiQC",
            subcategory="{sample}",
        ),
    log:
        "{folder}/{sample}/log/multiqc/{sample}.log",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        "multiqc {input.multiqc_input} --outdir {output.outdir}"
