
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


rule bwa_index:
    input:
        ancient(config["references_data"][config["reference"]]["reference_fasta"]),
    output:
        idx=multiext(
            config["references_data"][config["reference"]]["reference_fasta"],
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "{}.log".format(
            config["references_data"][config["reference"]]["reference_fasta"]
        ),
    params:
        algorithm="bwtsw",
    threads: 16
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    wrapper:
        "v1.7.0/bio/bwa/index"


rule bwa_strandseq_to_reference_alignment:
    input:
        mate1="{folder}/{sample}/fastq/{cell}.1.fastq.gz",
        mate2="{folder}/{sample}/fastq/{cell}.2.fastq.gz",
        ref="{ref}".format(
            ref=config["references_data"][config["reference"]]["reference_fasta"]
        ),
        ref_index="{ref}.ann".format(
            ref=config["references_data"][config["reference"]]["reference_fasta"]
        ),
    output:
        bam="{folder}/{sample}/all/{cell}.bam",
    log:
        bwa="{folder}/{sample}/log/{cell}.bwa.log",
        samtools="{folder}/{sample}/log/{cell}.samtools.log",
    threads: 6
    params:
        idx_prefix=lambda wildcards, input: input.ref_index.rsplit(".", 1)[0],
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "bwa mem -t {threads}"
        ' -R "@RG\\tID:{wildcards.cell}\\tPL:Illumina\\tSM:{wildcards.sample}"'
        " -v 2 {input.ref} {input.mate1} {input.mate2} 2> {log.bwa} | "
        " samtools view -b /dev/stdin > {output.bam} 2> {log.samtools}"


rule samtools_sort_bam:
    input:
        "{folder}/{sample}/all/{cell}.bam",
    output:
        temp("{folder}/{sample}/all/{cell}.sort.bam"),
    log:
        "{folder}/{sample}/log/samtools_sort/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
        time="10:00:00",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools sort -O BAM -o {output} {input} 2>&1 > {log}"


rule mark_duplicates:
    input:
        bam="{folder}/{sample}/all/{cell}.sort.bam",
    output:
        "{folder}/{sample}/all/{cell}.sort.mdup.bam",
    log:
        "{folder}/{sample}/log/markdup/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
        time="10:00:00",
    shell:
        "sambamba markdup {input.bam} {output} 2>&1 > {log}"


rule generate_exclude_file_for_mosaic_count:
    input:
        bam=lambda wc: expand(
            "{folder}/{sample}/all/{cell}.sort.mdup.bam",
            folder=config["input_bam_location"],
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
            "{folder}/{sample}/all/{cell}.sort.mdup.bam",
            folder=config["input_bam_location"],
            sample=wc.sample,
            cell=cell_per_sample[str(wc.sample)],
        ),
        bai=lambda wc: expand(
            "{folder}/{sample}/all/{cell}.sort.mdup.bam.bai",
            folder=config["input_bam_location"],
            sample=wc.sample,
            cell=cell_per_sample[str(wc.sample)],
        ),
        excl="{folder}/{sample}/config/chroms_to_exclude.txt",
    output:
        counts="{folder}/{sample}/ashleys_counts/{sample}.all.txt.fixme.gz",
        info="{folder}/{sample}/ashleys_counts/{sample}.all.info",
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


rule order_mosaic_count_output:
    input:
        "{folder}/{sample}/ashleys_counts/{sample}.all.txt.fixme.gz",
    output:
        "{folder}/{sample}/ashleys_counts/{sample}.all.txt.gz",
    log:
        "{folder}/log/ashleys_counts/{sample}/{sample}.log",
    run:
        df = pd.read_csv(input[0], compression="gzip", sep="\t")
        df = df.sort_values(by=["sample", "cell", "chrom", "start"])
        df.to_csv(output[0], index=False, compression="gzip", sep="\t")


rule plot_mosaic_counts:
    input:
        counts="{folder}/{sample}/ashleys_counts/{sample}.all.txt.gz",
        info="{folder}/{sample}/ashleys_counts/{sample}.all.info",
    output:
        "{folder}/{sample}/plots/ashleys_counts/CountComplete.classic.pdf",
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


if config["mosaicatcher_pipeline"] is False:

    rule samtools_index:
        input:
            "{folder}/{sample}/all/{cell}.sort.mdup.bam",
        output:
            "{folder}/{sample}/all/{cell}.sort.mdup.bam.bai",
        log:
            "{folder}/{sample}/log/samtools_index/{cell}.log",
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools index {input} 2>&1 > {log}"


if config["hand_selection"] is False:

    rule generate_features:
        input:
            bam=lambda wc: expand(
                "{folder}/{sample}/all/{cell}.sort.mdup.bam",
                folder=config["input_bam_location"],
                sample=wc.sample,
                cell=cell_per_sample[str(wc.sample)],
            ),
            bai=lambda wc: expand(
                "{folder}/{sample}/all/{cell}.sort.mdup.bam.bai",
                folder=config["input_bam_location"],
                sample=wc.sample,
                cell=cell_per_sample[str(wc.sample)],
            ),
            plot=expand(
                "{{folder}}/{{sample}}/plots/ashleys_counts/CountComplete.{plottype}.pdf",
                plottype=plottype_counts,
            ),
        output:
            "{folder}/{sample}/predictions/ashleys_features.tsv",
        log:
            "{folder}/log/ashleys/{sample}/features.log",
        conda:
            "../envs/ashleys.yaml"
        threads: 64
        params:
            windows="5000000 2000000 1000000 800000 600000 400000 200000",
            jobs=64,
            extension=".sort.mdup.bam",
            folder=lambda wildcards, input: "{}all".format(input.bam[0].split("all")[0]),
        resources:
            mem_mb=get_mem_mb_heavy,
            time="10:00:00",
        shell:
            "ashleys -j {params.jobs} features -f {params.folder} -w {params.windows} -o {output} --recursive_collect -e {params.extension}"

    rule predict:
        input:
            folder="{folder}/{sample}/predictions/ashleys_features.tsv",
        output:
            "{folder}/{sample}/cell_selection/labels_raw.tsv",
        log:
            "{folder}/log/ashleys/{sample}/prediction_ashleys.log",
        conda:
            "../envs/ashleys.yaml"
        params:
            model_default="./workflow/ashleys_models/svc_default.pkl",
            model_stringent="./workflow/ashleys_models/svc_stringent.pkl",
        resources:
            mem_mb=get_mem_mb,
            time="10:00:00",
        shell:
            "ashleys predict -p {input.folder} -o {output} -m {params.model_default}"


elif config["hand_selection"] is True:

    localrules:
        notebook_hand_selection,

    rule notebook_hand_selection:
        input:
            pdf_raw=expand(
                "{{folder}}/{{sample}}/plots/ashleys_counts/CountComplete.{plottype}.pdf",
                plottype=plottype_counts,
            ),
            info="{folder}/{sample}/ashleys_counts/{sample}.all.info",
        output:
            folder="{folder}/{sample}/cell_selection/labels_raw.tsv",
        log:
            "{folder}/log/hand_selection/{sample}/prediction_probabilities.log",
        params:
            cell_per_sample=cell_per_sample,
        conda:
            "../envs/notebook.yaml"
        notebook:
            "../notebooks/hand_selection.py.ipynb"


if config["use_light_data"] is False:

    rule cp_predictions:
        input:
            folder="{folder}/{sample}/cell_selection/labels_raw.tsv",
        output:
            folder="{folder}/{sample}/cell_selection/labels.tsv",
        log:
            "{folder}/log/cp_predictions/{sample}.log",
        conda:
            "../envs/ashleys.yaml"
        shell:
            "cp {input.folder} {output.folder} > {log} 2>&1"


elif config["use_light_data"] is True:

    rule dev_all_cells_correct:
        input:
            folder="{folder}/{sample}/cell_selection/labels_raw.tsv",
        output:
            folder="{folder}/{sample}/cell_selection/labels.tsv",
        log:
            "{folder}/log/dev_all_cells_correct/{sample}.log",
        conda:
            "../envs/mc_base.yaml"
        script:
            "../scripts/utils/dev_all_cells_correct.py"
