rule fastqc:
    input:
        "{path}/{sample}/fastq/{cell}.{pair}.fastq.gz",
    output:
        html=report(
            "{path}/{sample}/fastqc/{cell}_{pair}_fastqc.html",
            category="FastQC",
            subcategory="{sample}",
            labels={"Sample": "{sample}", "Cell": "{cell}", "Pair": "{pair}"},
        ),
        zip="{path}/{sample}/fastqc/{cell}_{pair}_fastqc.zip",
    params:
        "--quiet",
    log:
        "{path}/log/fastqc/{sample}/{cell}_{pair}.log",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
        time="10:00:00",
    wrapper:
        "v1.7.0/bio/fastqc"



rule bwa_index:
    input:
        config["reference"],
    output:
        idx=multiext(config["reference"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "{}.log".format(config["reference"]),
    params:
        algorithm="bwtsw",
    resources:
        mem_mb=get_mem_mb_heavy,
        time="10:00:00",
    wrapper:
        "v1.7.0/bio/bwa/index"


rule bwa_strandseq_to_reference_alignment:
    input:
        mate1="{path}/{sample}/fastq/{cell}.1.fastq.gz",
        mate2="{path}/{sample}/fastq/{cell}.2.fastq.gz",
        ref="{ref}".format(ref=config["reference"]),
        ref_index="{ref}.ann".format(ref=config["reference"]),
    output:
        bam="{path}/{sample}/all/{cell}.bam",
    log:
        bwa="{path}/{sample}/log/{cell}.bwa.log",
        samtools="{path}/{sample}/log/{cell}.samtools.log",
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
        "{path}/{sample}/all/{cell}.bam",
    output:
        temp("{path}/{sample}/all/{cell}.sort.bam"),
    log:
        "{path}/{sample}/log/samtools_sort/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
        time="10:00:00",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    shell:
        "samtools sort -O BAM -o {output} {input} 2>&1 > {log}"


rule mark_duplicates:
    input:
        bam="{path}/{sample}/all/{cell}.sort.bam",
    output:
        "{path}/{sample}/all/{cell}.sort.mdup.bam",
    log:
        "{path}/{sample}/log/markdup/{cell}.log",
    conda:
        "../envs/mc_bioinfo_tools.yaml"
    resources:
        mem_mb=get_mem_mb,
        time="10:00:00",
    shell:
        "sambamba markdup {input.bam} {output} 2>&1 > {log}"


if config["mosaicatcher_pipeline"] is False:

    rule samtools_index:
        input:
            "{path}/{sample}/all/{cell}.sort.mdup.bam",
        output:
            "{path}/{sample}/all/{cell}.sort.mdup.bam.bai",
        log:
            "{path}/{sample}/log/samtools_index/{cell}.log",
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools index {input} 2>&1 > {log}"
### ASHLEYS AUTOMATED ANALYSIS



if config["hand_selection"] is False:

    rule generate_features:
        input:
            # ashleys="{path}/config/ashleys_install_success.txt",
            bam=expand(
                "{path}/{sample}/all/{cell}.sort.mdup.bam",
                zip,
                path=input_bam_location_expand,
                sample=samples_expand,
                cell=cell_expand,
            ),
        output:
            "{path}/{sample}/predictions/ashleys_features.tsv",
        log:
            "{path}/log/ashleys/{sample}/features.log",
        conda:
            "../envs/ashleys.yaml"
        params:
            windows="5000000 2000000 1000000 800000 600000 400000 200000",
            jobs=23,
            extension=".sort.mdup.bam",
            path=lambda wildcards, input: "{}all".format(input.bam[0].split("all")[0]),
        resources:
            mem_mb=get_mem_mb_heavy,
            time="10:00:00",
        shell:
            "ashleys -j {params.jobs} features -f {params.path} -w {params.windows} -o {output} --recursive_collect -e {params.extension}"

    rule predict:
        input:
            path="{path}/{sample}/predictions/ashleys_features.tsv",
        output:
            "{path}/{sample}/cell_selection/labels_raw.tsv",
        log:
            "{path}/log/ashleys/{sample}/prediction_ashleys.log",
        conda:
            "../envs/ashleys.yaml"
        params:
            model_default="./workflow/ashleys_models/svc_default.pkl",
            model_stringent="./workflow/ashleys_models/svc_stringent.pkl",
        resources:
            mem_mb=get_mem_mb,
            time="10:00:00",
        shell:
            "ashleys predict -p {input.path} -o {output} -m {params.model_default}"
########################################################



#                      DEV PART
########################################################
### HAND SELECTION VIA JUPYTER NB
elif config["hand_selection"] is True:

    rule generate_exclude_file_for_mosaic_count:
        input:
            ancient(
                "{path}/config/config_df_ashleys.tsv".format(
                    path=config["input_bam_location"]
                )
            ),
            # ancient("config/samples.tsv"),
            bam=config["input_bam_location"],
        output:
            "{path}/config/exclude_file",
        log:
            "{path}/log/config/exclude_file.log",
        conda:
            "../envs/mc_base.yaml"
        params:
            chroms=config["chromosomes"],
        script:
            "../scripts/utils/generate_exclude_file.py"

    rule mosaic_count:
        input:
            bam=lambda wc: expand(
                "{path}/{sample}/all/{cell}.sort.mdup.bam",
                path=config["input_bam_location"],
                sample=samples,
                cell=cell_per_sample[str(wc.sample)]
                if wc.sample in cell_per_sample
                else "FOOBAR",
            ),
            bai=lambda wc: expand(
                "{path}/{sample}/all/{cell}.sort.mdup.bam.bai",
                path=config["input_bam_location"],
                sample=samples,
                cell=cell_per_sample[str(wc.sample)],
            )
            if wc.sample in cell_per_sample
            else "FOOBAR",
            excl="{path}/config/exclude_file",
        output:
            counts="{path}/{sample}/ashleys_counts/{sample}.all.txt.fixme.gz",
            info="{path}/{sample}/ashleys_counts/{sample}.all.info",
        log:
            "{path}/log/counts/{sample}/mosaic_count.log",
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
            "{path}/{sample}/ashleys_counts/{sample}.all.txt.fixme.gz",
        output:
            "{path}/{sample}/ashleys_counts/{sample}.all.txt.gz",
        log:
            "{path}/log/ashleys_counts/{sample}/{sample}.log",
        run:
            df = pd.read_csv(input[0], compression="gzip", sep="\t")
            df = df.sort_values(by=["sample", "cell", "chrom", "start"])
            df.to_csv(output[0], index=False, compression="gzip", sep="\t")

    rule plot_mosaic_counts:
        input:
            counts="{path}/{sample}/ashleys_counts/{sample}.all.txt.gz",
            info="{path}/{sample}/ashleys_counts/{sample}.all.info",
        output:
            "{path}/{sample}/plots/ashleys_counts/CountComplete_{sample}.pdf",
        log:
            "{path}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            """
            LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """

    rule notebook_hand_selection:
        input:
            pdf_raw="{path}/{sample}/plots/ashleys_counts/CountComplete_{sample}.pdf",
        output:
            path="{path}/{sample}/cell_selection/labels_raw.tsv",
        log:
            "{path}/log/hand_selection/{sample}/prediction_probabilities.log",
        params:
            cell_per_sample=cell_per_sample,
        conda:
            "../envs/notebook.yaml"
        notebook:
            "../notebooks/hand_selection.py.ipynb"


if config["use_light_data"] is False:

    rule cp_predictions:
        input:
            path="{path}/{sample}/cell_selection/labels_raw.tsv",
        output:
            path="{path}/{sample}/cell_selection/labels.tsv",
        log:
            "{path}/log/cp_predictions/{sample}.log",
        conda:
            "../envs/ashleys.yaml"
        shell:
            "cp {input.path} {output.path} > {log} 2>&1"
# BM cells 05 & 12 I EXAMPLE DATA WERE IDENTIFIED AS NOT POSSIBLE TO BE PROCESSED BY MOSAIC COUNT



elif config["use_light_data"] is True:

    rule dev_all_cells_correct:
        input:
            path="{path}/{sample}/cell_selection/labels_raw.tsv",
        output:
            path="{path}/{sample}/cell_selection/labels.tsv",
        log:
            "{path}/log/dev_all_cells_correct/{sample}.log",
        run:
            df = pd.read_csv(input.path, sep="\t")
            df["prediction"] = 1
            df["probability"] = 1
            df.loc[df["cell"].str.contains("05"), "prediction"] = 0
            df.loc[df["cell"].str.contains("12"), "prediction"] = 0
            df.to_csv(output.path, sep="\t", index=False)
