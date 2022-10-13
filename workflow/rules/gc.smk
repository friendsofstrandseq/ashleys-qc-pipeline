if config["GC_analysis"] is True:

    rule mergeBams:
        input:
            lambda wc: expand(
                "{folder}/{sample}/bam/{bam}.sort.mdup.bam",
                folder=config["data_location"],
                sample=wc.sample,
                bam=cell_per_sample[wc.sample],
            ),
        output:
            "{folder}/{sample}/merged_bam/merged.raw.bam",
        log:
            "{folder}/log/merged_bam/{sample}.log",
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        threads: 32
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools merge -@ {threads} {output} {input} 2>&1 > {log}"

    rule mergeSortBams:
        input:
            "{folder}/{sample}/merged_bam/merged.raw.bam",
        output:
            "{folder}/{sample}/merged_bam/merged.bam",
        log:
            "{folder}/log/merged_bam/{sample}.mergeSortBams.log",
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        threads: 32 
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools sort -@ {threads} -o {output} {input} 2>&1 > {log}"

    rule index_merged_bam:
        input:
            "{folder}/{sample}/merged_bam/merged.bam",
        output:
            "{folder}/{sample}/merged_bam/merged.bam.bai",
        log:
            "{folder}/log/{sample}/merged_bam/index_merged_bam.log",
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            "samtools index {input} > {log} 2>&1"

    ruleorder: mergeBams > mergeSortBams > mergeBams_plate_row > mergeSortBams_plate_row > index_merged_bam > index_merged_bam_plate_row > alfred_merged > alfred_plate_row > alfred_sc > alfred_table > alfred_table_merged > alfred_table_plate_row
    
    rule alfred_merged:
        input:
            merged_bam="{folder}/{sample}/merged_bam/merged.bam",
            merged_bam_bai="{folder}/{sample}/merged_bam/merged.bam.bai",
            fasta=config["references_data"][config["reference"]]["reference_fasta"],
            fasta_index="{fasta}.fai".format(
                fasta=config["references_data"][config["reference"]]["reference_fasta"]
            ),
        output:
            alfred_json="{folder}/{sample}/alfred/MERGE/merged_bam.merge.json.gz",
            alfred_tsv="{folder}/{sample}/alfred/MERGE/merged_bam.merge.tsv.gz",
        log:
            "{folder}/{sample}/log/alfred/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/dev/mc_bioinfo_tools.yaml"
        shell:
            """
            alfred qc -r {input.fasta} -j {output.alfred_json} -o {output.alfred_tsv} {input.merged_bam}
            """

    rule alfred_sc:
        input:
            bam="{folder}/{sample}/bam/{bam}.sort.mdup.bam",
            bam_bai="{folder}/{sample}/bam/{bam}.sort.mdup.bam.bai",
            fasta=config["references_data"][config["reference"]]["reference_fasta"],
            fasta_index="{fasta}.fai".format(
                fasta=config["references_data"][config["reference"]]["reference_fasta"]
            ),
        output:
            alfred_json="{folder}/{sample}/alfred/{bam}.json.gz",
            alfred_tsv="{folder}/{sample}/alfred/{bam}.tsv.gz",
        log:
            "{folder}/{sample}/log/alfred/{bam}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/dev/mc_bioinfo_tools.yaml"
        shell:
            """
            alfred qc -r {input.fasta} -j {output.alfred_json} -o {output.alfred_tsv} {input.bam}
            """


    rule alfred_table:
        input:
            bam="{folder}/{sample}/alfred/{bam}.tsv.gz",
        output:
            "{folder}/{sample}/alfred/{bam}.table",
        log:
            "{folder}/{sample}/log/alfred_table/{bam}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            zcat {input} | grep "^GC" > {output}
            """
    rule alfred_table_plate_row:
        input:
            bam="{folder}/{sample}/alfred/PLATE_ROW/{row}.row.tsv.gz",
        output:
            "{folder}/{sample}/alfred/PLATE_ROW/{row}.row.table",
        log:
            "{folder}/{sample}/log/alfred_table/PLATE_ROW/{row}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/mc_base.yaml"
        shell:
            """
            zcat {input} | grep "^GC" > {output}
            """

    rule alfred_table_merged:
        input:
            bam="{folder}/{sample}/alfred/MERGE/merged_bam.merge.tsv.gz",
        output:
            "{folder}/{sample}/alfred/MERGE/merged_bam.merge.table",
        log:
            "{folder}/{sample}/log/alfred_table/merged_bam.log",
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
            counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
        output:
            counts_vst="{folder}/{sample}/counts/GC_correction/{sample}.txt.VST.gz",
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
            counts_vst="{folder}/{sample}/counts/GC_correction/{sample}.txt.VST.gz",
        output:
            counts_vst_gc="{folder}/{sample}/counts/GC_correction/{sample}.txt.VST.GC.gz",
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
            counts_vst_gc="{folder}/{sample}/counts/GC_correction/{sample}.txt.VST.GC.gz",
        output:
            counts_vst_gc_scaled="{folder}/{sample}/counts/GC_correction/{sample}.txt.VST.GC.scaled.gz",
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
            counts="{folder}/{sample}/counts/GC_correction/{sample}.txt.VST.GC.scaled.gz",
            info="{folder}/{sample}/counts/{sample}.info_raw",
        output:
            "{folder}/{sample}/plots/counts/CountComplete.GC_corrected.pdf",
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
            # table="{folder}/{sample}/alfred/{sample}.table",
            table="{folder}/{sample}/alfred/{bam}.table",
        output:
            gcdist_plot=report(
                "{folder}/{sample}/plots/alfred/{bam}_gc_dist.png",
                category="GC analysis",
                labels={"Sample": "{sample}", "Cell(s)" : "{bam}", "Type": "GC distribution"},
            ),
            gcdevi_plot=report(
                "{folder}/{sample}/plots/alfred/{bam}_gc_devi.png",
                category="GC analysis",
                labels={"Sample": "{sample}", "Cell(s)" : "{bam}", "Type": "GC deviation"},
            ),
        log:
            "{folder}/{sample}/log/alfred_plot/{bam}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/rtools.yaml"
        script:
            "../scripts/GC/gc.R"



    rule mergeBams_plate_row:
        input:
            # "{folder}/{sample}/all/{bam}.sort.mdup.bam",
            lambda wc: expand(
                    "{folder}/{sample}/bam/{bam}.sort.mdup.bam",
                    folder=config["data_location"],
                    sample=wc.sample,
                    bam=d[wc.sample][wc.row]
            ),
        output:
            "{folder}/{sample}/merged_bam/PLATE_ROW/{row}.platerow.raw.bam",
        log:
            "{folder}/log/merged_bam/{sample}.{row}.log",
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        threads: 32
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools merge -@ {threads} {output} {input} 2>&1 > {log}"

    rule mergeSortBams_plate_row:
        input:
            "{folder}/{sample}/merged_bam/PLATE_ROW/{row}.platerow.raw.bam",
        output:
            "{folder}/{sample}/merged_bam/PLATE_ROW/{row}.platerow.bam",
        log:
            "{folder}/log/merged_bam/{sample}.{row}.mergeSortBams.log",
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        threads: 32
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools sort -@ {threads} -o {output} {input} 2>&1 > {log}"

    rule index_merged_bam_plate_row:
        input:
            "{folder}/{sample}/merged_bam/PLATE_ROW/{row}.platerow.bam",
        output:
            "{folder}/{sample}/merged_bam/PLATE_ROW/{row}.platerow.bam.bai",
        log:
            "{folder}/log/merged_bam/{sample}/{row}.log",
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            "samtools index {input} > {log} 2>&1"
            
    rule alfred_plate_row:
        input:
            bam="{folder}/{sample}/merged_bam/PLATE_ROW/{row}.platerow.bam",
            bam_bai="{folder}/{sample}/merged_bam/PLATE_ROW/{row}.platerow.bam.bai",
            fasta=config["references_data"][config["reference"]]["reference_fasta"],
            fasta_index="{fasta}.fai".format(
                fasta=config["references_data"][config["reference"]]["reference_fasta"]
            ),
        output:
            alfred_json="{folder}/{sample}/alfred/PLATE_ROW/{row}.row.json.gz",
            alfred_tsv="{folder}/{sample}/alfred/PLATE_ROW/{row}.row.tsv.gz",
        log:
            "{folder}/{sample}/log/alfred/{row}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/dev/mc_bioinfo_tools.yaml"
        shell:
            """
            alfred qc -r {input.fasta} -j {output.alfred_json} -o {output.alfred_tsv} {input.bam}
            """

    rule alfred_plot_merge:
        input:
            # table="{folder}/{sample}/alfred/{sample}.table",
            table="{folder}/{sample}/alfred/MERGE/merged_bam.merge.table",
        output:
            gcdist_plot=report(
                "{folder}/{sample}/plots/alfred/MERGE/merged_bam_gc_dist.merge.png",
                category="GC analysis",
                labels={"Sample": "{sample}", "Type": "GC distribution"},
            ),
            gcdevi_plot=report(
                "{folder}/{sample}/plots/alfred/MERGE/merged_bam_gc_devi.merge.png",
                category="GC analysis",
                labels={"Sample": "{sample}",  "Type": "GC deviation"},
            ),
        log:
            "{folder}/{sample}/log/alfred_plot/merge_bam.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/rtools.yaml"
        script:
            "../scripts/GC/gc.R"

    rule alfred_plot_plate_row:
        input:
            # table="{folder}/{sample}/alfred/{sample}.table",
            table="{folder}/{sample}/alfred/PLATE_ROW/{row}.row.table",
        output:
            gcdist_plot=report(
                "{folder}/{sample}/plots/alfred/PLATE_ROW/{row}_gc_dist.row.png",
                category="GC analysis",
                labels={"Sample": "{sample}", "Cell(s)" : "{row}", "Type": "GC distribution"},
            ),
            gcdevi_plot=report(
                "{folder}/{sample}/plots/alfred/PLATE_ROW/{row}_gc_devi.row.png",
                category="GC analysis",
                labels={"Sample": "{sample}", "Cell(s)" : "{row}", "Type": "GC deviation"},
            ),
        log:
            "{folder}/{sample}/log/alfred_plot/{row}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/rtools.yaml"
        script:
            "../scripts/GC/gc.R"