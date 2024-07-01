## Rules to perform GC analysis & correction on Strand-Seq libraries
## ---------------------------------------------------------------
## fastqc: QC analysis of FASTQ files
## bwa_index: index fasta file using the BW transformation to perform next reads mapping
## bwa_strandseq_to_reference_alignment: mapping against FASTA reference (based on reference selected (hg19/hg38/T2T))
## samtools_sort_bam: sorting bam files
## mark_duplicates: mark duplicates in bam files
## *samtools index: index bam files (only if not loaded as a module into mosaicatcher_pipeline)
## generate_features/predict: features creation & prediction using ashleys-qc ML method to detect high/low quality libraries
## notebook_hand_selection: fire a jupyter notebook that allow hand selection of low quality cells based on QC plots


if config["genecore"] is True and config["genecore_date_folder"]:
    if config["mosaicatcher_pipeline"] is False:

        localrules:
            genecore_symlink,
            symlink_bam_ashleys,

    rule genecore_symlink:
        input:
            lambda wc: df_config_files.loc[
                (df_config_files["Sample"] == wc.sample)
                & (df_config_files["File"] == "{}.{}".format(wc.cell, wc.pair))
            ]["Genecore_path"]
            .unique()
            .tolist(),
        output:
            "{folder}/{sample}/fastq/{cell}.{pair}.fastq.gz",
        log:
            "{folder}/log/genecore_symlink/{sample}/{cell}_{pair}.log",
        shell:
            "ln -f -s {input} {output}"

    ruleorder: genecore_symlink > bwa_strandseq_to_reference_alignment


# if config["use_light_data"] is False:
localrules:
    symlink_bam_ashleys,


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


if config["mosaicatcher_pipeline"] is False:

    ruleorder: bwa_strandseq_to_reference_alignment > samtools_sort_bam > mark_duplicates > samtools_index

else:

    ruleorder: bwa_strandseq_to_reference_alignment > samtools_sort_bam > mark_duplicates


if config["paired_end"] is True:

    rule bwa_strandseq_to_reference_alignment:
        input:
            mate1="{folder}/{sample}/fastq/{cell}.1.fastq.gz",
            mate2="{folder}/{sample}/fastq/{cell}.2.fastq.gz",
            ref="{ref}".format(
                ref=config["references_data"][config["reference"]]["reference_fasta"]
            ),
            ref_index=multiext(
                config["references_data"][config["reference"]]["reference_fasta"],
                ".amb",
                ".ann",
                ".bwt",
                ".pac",
                ".sa",
            ),
        output:
            bam="{folder}/{sample}/bam/{cell}.bam.raw",
        log:
            bwa="{folder}/{sample}/log/{cell}.bwa.log",
            samtools="{folder}/{sample}/log/{cell}.samtools.log",
        threads: 6
        params:
            idx_prefix=lambda wildcards, input: input.ref_index[0].rsplit(".", 1)[0],
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        conda:
            "../envs/ashleys_base.yaml"
        shell:
            "bwa mem -t {threads}"
            ' -R "@RG\\tID:{wildcards.cell}\\tPL:Illumina\\tSM:{wildcards.sample}"'
            " -v 2 {input.ref} {input.mate1} {input.mate2} 2> {log.bwa} | "
            " samtools view -b /dev/stdin > {output.bam} 2> {log.samtools}"

else:

    rule bwa_strandseq_to_reference_alignment:
        input:
            mate1="{folder}/{sample}/fastq/{cell}.1.fastq.gz",
            ref="{ref}".format(
                ref=config["references_data"][config["reference"]]["reference_fasta"]
            ),
            ref_index=multiext(
                config["references_data"][config["reference"]]["reference_fasta"],
                ".amb",
                ".ann",
                ".bwt",
                ".pac",
                ".sa",
            ),
        output:
            bam="{folder}/{sample}/bam/{cell}.bam.raw",
        log:
            bwa="{folder}/{sample}/log/{cell}.bwa.log",
            samtools="{folder}/{sample}/log/{cell}.samtools.log",
        threads: 6
        params:
            idx_prefix=lambda wildcards, input: input.ref_index[0].rsplit(".", 1)[0],
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        conda:
            "../envs/ashleys_base.yaml"
        shell:
            "bwa mem -t {threads}"
            ' -R "@RG\\tID:{wildcards.cell}\\tPL:Illumina\\tSM:{wildcards.sample}"'
            " -v 2 {input.ref} {input.mate1} 2> {log.bwa} | "
            " samtools view -b /dev/stdin > {output.bam} 2> {log.samtools}"


rule samtools_sort_bam:
    input:
        "{folder}/{sample}/bam/{cell}.bam.raw",
    output:
        "{folder}/{sample}/bam/{cell}.bam.sort",
    log:
        "{folder}/{sample}/log/samtools_sort/{cell}.log",
    resources:
        mem_mb=get_mem_mb,
        time="01:00:00",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        "samtools sort -O BAM -o {output} {input} 2>&1 > {log}"


rule mark_duplicates:
    input:
        bam="{folder}/{sample}/bam/{cell}.bam.sort",
    output:
        "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
    log:
        "{folder}/{sample}/log/markdup/{cell}.log",
    conda:
        "../envs/ashleys_base.yaml"
    resources:
        mem_mb=get_mem_mb_heavy,
        time="01:00:00",
    shell:
        "sambamba markdup {input.bam} {output} 2>&1 > {log}"


if config["mosaicatcher_pipeline"] is False:

    rule samtools_index:
        input:
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam",
        output:
            "{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
        log:
            "{folder}/{sample}/log/samtools_index/{cell}.log",
        conda:
            "../envs/ashleys_base.yaml"
        shell:
            "samtools index {input} 2>&1 > {log}"

    rule gunzip_fasta:
        input:
            ancient("{file}.fa.gz"),
        output:
            "{file}.fa",
        log:
            "{file}.log",
        conda:
            "../envs/ashleys_base.yaml"
        resources:
            mem_mb=get_mem_mb_heavy,
        shell:
            "gunzip -cd {input} > {output}"


rule symlink_bam_ashleys:
    input:
        bam="{folder}/{sample}/bam/{cell}.sort.mdup.bam",
        bai="{folder}/{sample}/bam/{cell}.sort.mdup.bam.bai",
    output:
        bam="{folder}/{sample}/bam_ashleys/{cell}.sort.mdup.bam",
        bai="{folder}/{sample}/bam_ashleys/{cell}.sort.mdup.bam.bai",
    log:
        "{folder}/log/symlink_selected_bam/{sample}/{cell}.log",
    conda:
        "../envs/ashleys_base.yaml"
    script:
        "../scripts/utils/symlink_selected_bam.py"


if config["bypass_ashleys"] is False:

    rule generate_features:
        input:
            bam=selected_input_bam,
        output:
            "{folder}/{sample}/predictions/ashleys_features.tsv",
        log:
            "{folder}/log/ashleys/{sample}/features.log",
        conda:
            "../envs/ashleys_base.yaml"
        threads: 64
        params:
            windows="5000000 2000000 1000000 800000 600000 400000 200000",
            extension=".sort.mdup.bam",
            folder=lambda wildcards, input: "{}bam_ashleys".format(
                input.bam[0].split("bam_ashleys")[0]
            ),
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        shell:
            "ashleys -j {threads} features -f {params.folder} -w {params.windows} -o {output} --recursive_collect -e {params.extension}"

    rule predict:
        input:
            folder="{folder}/{sample}/predictions/ashleys_features.tsv",
        output:
            "{folder}/{sample}/cell_selection/labels_raw.tsv",
        log:
            "{folder}/log/ashleys/{sample}/prediction_ashleys.log",
        conda:
            "../envs/ashleys_base.yaml"
        params:
            model_default="./workflow/ashleys_models/svc_default.pkl",
            model_stringent="./workflow/ashleys_models/svc_stringent.pkl",
        resources:
            mem_mb=get_mem_mb,
            time="01:00:00",
        shell:
            "ashleys predict -p {input.folder} -o {output} -m {params.model_default}"

else:

    rule generate_default_labels:
        input:
            bam=selected_input_bam,
        output:
            "{folder}/{sample}/cell_selection/labels_raw.tsv",
        log:
            "{folder}/log/generate_default_labels/{sample}.log",
        conda:
            "../envs/ashleys_base.yaml"
        shell:
            """
            echo "cell\tprediction\tprobability\tsample" > {output}
            for bam in {input.bam} ; do
                # remove path only
                cell=$(basename $bam)
                echo -e "$cell\t1\t1\t{wildcards.sample}" >> {output}
            done
            """


# localrules:
#     notebook_hand_selection,

# rule notebook_hand_selection:
#     input:
#         pdf=expand(
#             "{{folder}}/{{sample}}/plots/counts/CountComplete.{plottype}.pdf",
#             plottype=plottype_counts,
#         ),
#         info="{folder}/{sample}/counts/{sample}.info_raw",
#         ashleys_labels="{folder}/{sample}/cell_selection/labels_raw.tsv",
#     output:
#         folder="{folder}/{sample}/cell_selection/labels_notebook.tsv",
#     log:
#         "{folder}/log/hand_selection/{sample}/prediction_probabilities.log",
#     params:
#         cell_per_sample=cell_per_sample,
#     conda:
#         "../envs/ashleys_notebook.yaml"
#     container:
#         None
#     notebook:
#         "../notebooks/hand_selection.py.ipynb"

# else:


rule copy_labels:
    input:
        labels="{folder}/{sample}/cell_selection/labels_raw.tsv",
    output:
        folder="{folder}/{sample}/cell_selection/labels_notebook.tsv",
    log:
        "{folder}/log/positive_control_bypass/{sample}.log",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        "cp {input} {output}"


if config["use_light_data"] is False:

    rule positive_negative_control_bypass:
        input:
            labels="{folder}/{sample}/cell_selection/labels_notebook.tsv",
            info="{folder}/{sample}/counts/{sample}.info_raw",
        output:
            labels_corrected="{folder}/{sample}/cell_selection/labels_positive_control_corrected.tsv",
            bypass_cell="{folder}/{sample}/config/bypass_cell.txt",
        log:
            "{folder}/log/positive_control_bypass/{sample}.log",
        conda:
            "../envs/ashleys_base.yaml"
        script:
            "../scripts/utils/positive_negative_control_bypass.py"

    checkpoint tune_predictions_based_on_threshold:
        input:
            "{folder}/{sample}/cell_selection/labels_positive_control_corrected.tsv",
        output:
            "{folder}/{sample}/cell_selection/labels.tsv",
        log:
            "{folder}/log/cp_predictions/{sample}.log",
        conda:
            "../envs/ashleys_base.yaml"
        script:
            "../scripts/utils/tune_predictions_based_on_threshold.py"

    rule plot_plate:
        input:
            labels="{folder}/{sample}/cell_selection/labels.tsv",
        output:
            predictions=report(
                "{folder}/{sample}/plots/plate/ashleys_plate_predictions.pdf",
                category="Ashleys plate plots",
                subcategory="{sample}",
                labels={"Sample": "{sample}", "Plot Type": "Predictions"},
            ),
            probabilities=report(
                "{folder}/{sample}/plots/plate/ashleys_plate_probabilities.pdf",
                category="Ashleys plate plots",
                subcategory="{sample}",
                labels={"Sample": "{sample}", "Plot Type": "Probabilities"},
            ),
            well_table="{folder}/{sample}/plots/plate/ashleys_well_table.tsv",
        log:
            "{folder}/log/plot_plate/{sample}.log",
        conda:
            "../envs/ashleys_rtools.yaml"
        script:
            "../scripts/plotting/plot_plate.R"

elif config["use_light_data"] is True:

    rule dev_all_cells_correct:
        input:
            folder="{folder}/{sample}/cell_selection/labels_notebook.tsv",
        output:
            folder="{folder}/{sample}/cell_selection/labels.tsv",
        log:
            "{folder}/log/dev_all_cells_correct/{sample}.log",
        conda:
            "../envs/ashleys_base.yaml"
        script:
            "../scripts/utils/dev_all_cells_correct.py"


if config["publishdir"] != "":

    rule publishdir_outputs_ashleys:
        input:
            list_publishdir=publishdir_fct,
        output:
            touch("{folder}/{sample}/config/publishdir_outputs_ashleys.ok"),
        log:
            "{folder}/log/publishdir_outputs_ashleys/{sample}.log",
        conda:
            "../envs/ashleys_base.yaml"
        script:
            "../scripts/utils/publishdir.py"


rule save_config:
    input:
        "config/config.yaml",
    output:
        "{folder}/{sample}/config/config_ashleys.yaml",
    log:
        "{folder}/log/save_config/{sample}.log",
    conda:
        "../envs/ashleys_base.yaml"
    resources:
        mem_mb=get_mem_mb,
    script:
        "../scripts/utils/dump_config.py"


rule ashleys_final_results:
    input:
        get_final_output,
    output:
        "{folder}/{sample}/config/ashleys_final_results.ok",
    log:
        "{folder}/log/ashleys_final_results/{sample}.log",
    conda:
        "../envs/ashleys_base.yaml"
    shell:
        "touch {output}"
