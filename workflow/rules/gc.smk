## Rules to perform GC analysis & correction on Strand-Seq libraries
## ---------------------------------------------------------------
## mergeBams/mergeBams_plate_row: merge all/fraction of bams
## mergeSortBams/mergeSortBams_plate_row: sort merged bam file
## index_merged_bam/index_merged_bam_plate_row: index merged bam file
## alfred_plate_row/alfred_merged/alfred_sc: alfred QC to retrieve stats on bam files
## alfred_table_plate_row/alfred_table_merged/alfred_table_sc: extract only GC rows
## alfred_plot_plate_row/alfred_plot_merged/alfred_plot_sc: plot statistics using R script
## VST_correction/GC_correction/counts_scaling: variance stabilizing transformation, GC correction & counts scaling based @MarcoCosenza method
## plot_mosaic_gc_norm_counts: plots QC counts plot after correction


if config["GC_analysis"] is True:


    rule library_size_normalisation:
        input:
            counts="{folder}/{sample}/counts/{sample}.txt.raw.gz",
        output:
            counts_scaled="{folder}/{sample}/counts/GC_correction/{sample}.txt.scaled.gz",
        log:
            "{folder}/{sample}/log/counts_scaling/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/ashleys_rtools.yaml"
        script:
            "../scripts/GC/library_size_normalisation.R"


    rule GC_correction:
        input:
            counts_scaled="{folder}/{sample}/counts/GC_correction/{sample}.txt.scaled.gz",
        output:
            counts_scaled_gc="{folder}/{sample}/counts/GC_correction/{sample}.txt.scaled.GC.gz",
            plot=report("{folder}/{sample}/plots/GC_correction/GC_correction_result_distribution.png",
                category="GC analysis",
                subcategory="{sample}",
                labels={
                    "Sample": "{sample}",
                    "Plot Type": "GC distribution",
                    "Cell/Row/Plate": "Row {row}",
                },
            ),
        log:
            "{folder}/{sample}/log/GC_correction/{sample}.log",
        params:
            gc_matrix="workflow/data/GC/GC_matrix_200000.txt",
            gc_min_reads=5,
            gc_n_subsample=1000,
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/ashleys_rtools.yaml"
        script:
            "../scripts/GC/GC_correction.R"


    rule VST_correction:
        input:
            counts_scaled_gc="{folder}/{sample}/counts/GC_correction/{sample}.txt.scaled.GC.gz",
        output:
            counts_scaled_gc_vst="{folder}/{sample}/counts/GC_correction/{sample}.txt.scaled.GC.VST.gz",
        log:
            "{folder}/{sample}/log/VST_correction/{sample}.log",
        resources:
            mem_mb=get_mem_mb,
        conda:
            "../envs/ashleys_rtools.yaml"
        script:
            "../scripts/GC/variance_stabilizing_transformation.R"

        
    rule populate_counts_GC:
        input:
            bin_bed="workflow/data/bin_200kb_all.bed",
            counts="{folder}/{sample}/counts/GC_correction/{sample}.txt.scaled.GC.VST.gz",
        output:
            populated_counts="{folder}/{sample}/counts/GC_correction/{sample}.txt.scaled.GC.VST.populated.gz",
        log:
            "{folder}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../envs/ashleys_base.yaml"
        resources:
            mem_mb=get_mem_mb,
        script:
            "../scripts/utils/populated_counts_for_qc_plot.py"

    rule plot_mosaic_gc_norm_counts:
        input:
            counts="{folder}/{sample}/counts/GC_correction/{sample}.txt.scaled.GC.VST.populated.gz",
            info="{folder}/{sample}/counts/{sample}.info_raw",
        output:
            "{folder}/{sample}/plots/counts/CountComplete.GC_corrected.pdf",
        log:
            "{folder}/{sample}/log/plot_mosaic_counts/{sample}.log",
        conda:
            "../envs/ashleys_rtools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            """
            LC_CTYPE=C Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output} > {log} 2>&1
            """
