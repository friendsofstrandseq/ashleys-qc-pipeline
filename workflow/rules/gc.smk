if config["GC_analysis"] is True:

    rule mergeBams:
        input:
            lambda wc: expand(
                "{input_folder}/{sample}/all/{bam}.sort.mdup.bam",
                input_folder=config["input_bam_location"],
                sample=wc.sample,
                bam=cell_per_sample[wc.sample],
            ),
        output:
            "{output_folder}/{sample}/ashleys_merged_bam/{sample}/merged.raw.bam",
        log:
            "{output_folder}/log/ashleys_merged_bam/{sample}.log",
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        threads: 10
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools merge -@ {threads} {output} {input} 2>&1 > {log}"

    rule mergeSortBams:
        input:
            "{output_folder}/{sample}/ashleys_merged_bam/{sample}/merged.raw.bam",
        output:
            "{output_folder}/{sample}/ashleys_merged_bam/{sample}/merged.bam",
        log:
            "{output_folder}/log/ashleys_merged_bam/{sample}.mergeSortBams.log",
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        threads: 10
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools sort -@ {threads} -o {output} {input} 2>&1 > {log}"

    rule index_merged_bam:
        input:
            "{output_folder}/{sample}/ashleys_merged_bam/{sample}/merged.bam",
        output:
            "{output_folder}/{sample}/ashleys_merged_bam/{sample}/merged.bam.bai",
        log:
            "{output_folder}/log/ashleys_merged_bam/{sample}/index_merged_bam.log",
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            "samtools index {input} > {log} 2>&1"

    rule alfred:
        input:
            merged_bam="{folder}/{sample}/ashleys_merged_bam/{sample}/merged.bam",
            merged_bam_bai="{folder}/{sample}/ashleys_merged_bam/{sample}/merged.bam.bai",
            fasta=config["references_data"][config["reference"]]["reference_fasta"],
            fasta_index="{fasta}.fai".format(
                fasta=config["references_data"][config["reference"]]["reference_fasta"]
            ),
        output:
            alfred_json="{folder}/{sample}/alfred/{sample}.json.gz",
            alfred_tsv="{folder}/{sample}/alfred/{sample}.tsv.gz",
        log:
            "{folder}/{sample}/log/alfred/{sample}.log",
        resources:
            mem_mb=get_mem_mb_heavy,
        conda:
            "../envs/dev/mc_bioinfo_tools.yaml"
        shell:
            """
            alfred qc -r {input.fasta} -j {output.alfred_json} -o {output.alfred_tsv} {input.merged_bam}
            """

    rule alfred_table:
        input:
            "{folder}/{sample}/alfred/{sample}.tsv.gz",
        output:
            "{folder}/{sample}/alfred/{sample}.table",
        log:
            "{folder}/{sample}/log/alfred_table/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            zcat {input} | grep "^GC" > {output}
            """

    rule VST_correction:
        input:
            counts="{folder}/{sample}/ashleys_counts/{sample}.all.txt.gz",
        output:
            counts_vst="{folder}/{sample}/counts/GC_correction/{sample}.all.VST.txt.gz",
        log:
            "{folder}/{sample}/log/VST_correction/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/dev/GC.yaml"
        script:
            "../scripts/GC/variance_stabilizing_transformation.R"

    rule GC_correction:
        input:
            counts_vst="{folder}/{sample}/counts/GC_correction/{sample}.all.VST.txt.gz",
        output:
            counts_vst_gc="{folder}/{sample}/counts/GC_correction/{sample}.all.VST.GC.txt.gz",
        log:
            "{folder}/{sample}/log/GC_correction/{sample}.log",
        params:
            gc_matrix="workflow/data/GC/GC_matrix_200000.txt",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/dev/GC.yaml"
        script:
            "../scripts/GC/GC_correction.R"

    rule counts_scaling:
        input:
            counts_vst_gc="{folder}/{sample}/counts/GC_correction/{sample}.all.VST.GC.txt.gz",
        output:
            counts_vst_gc_scaled="{folder}/{sample}/counts/GC_correction/{sample}.all.VST.GC.scaled.txt.gz",
        log:
            "{folder}/{sample}/log/counts_scaling/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/dev/GC.yaml"
        script:
            "../scripts/GC/counts_scaling.R"

    rule plot_mosaic_gc_norm_counts:
        input:
            counts="{folder}/{sample}/counts/GC_correction/{sample}.all.VST.GC.scaled.txt.gz",
            info="{folder}/{sample}/ashleys_counts/{sample}.all.info",
        output:
            "{folder}/{sample}/plots/ashleys_counts/CountComplete.GC_corrected.pdf",
        log:
            "{folder}/{sample}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../envs/rtools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            """
            LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """

    rule alfred_plot:
        input:
            table="{folder}/{sample}/alfred/{sample}.table",
        output:
            gcdist_plot=report(
                "{folder}/{sample}/plots/{sample}/alfred/gc_dist.png",
                category="GC analysis",
                labels={"Sample": "{sample}", "Type": "GC distribution"},
            ),
            gcdevi_plot=report(
                "{folder}/{sample}/plots/{sample}/alfred/gc_devi.png",
                category="GC analysis",
                labels={"Sample": "{sample}", "Type": "GC deviation"},
            ),
        log:
            "{folder}/{sample}/log/alfred_plot/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/rtools.yaml"
        script:
            "../scripts/GC/gc.R"
