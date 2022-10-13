
    rule mergeBams_plate_row:
        input:
            "{folder}/{sample}/all/{bam}.sort.mdup.bam",
            lambda wc: [
                expand(
                    "{folder}/{sample}/all/{bam}.sort.mdup.bam",
                    folder=config["input_bam_location"],
                    sample=wc.sample,
                    row=row
                    bam=d[wc.sample][row],
             )
                for row in d[wc.sample]
            ],
        output:
            "{folder}/{sample}/ashleys_merged_bam/PLATE_ROW/{row}.raw.bam",
        log:
            "{folder}/log/ashleys_merged_bam/{sample}.{row}.log",
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        threads: 10
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools merge -@ {threads} {output} {input} 2>&1 > {log}"

    rule mergeSortBams_plate_row:
        input:
            "{folder}/{sample}/ashleys_merged_bam/PLATE_ROW/{row}.raw.bam",
        output:
            "{folder}/{sample}/ashleys_merged_bam/PLATE_ROW/{row}.bam",
        log:
            "{folder}/log/ashleys_merged_bam/{sample}.{row}.mergeSortBams.log",
        resources:
            mem_mb=get_mem_mb_heavy,
            time="01:00:00",
        threads: 10
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        shell:
            "samtools sort -@ {threads} -o {output} {input} 2>&1 > {log}"

    rule index_merged_bam_plate_row:
        input:
            "{folder}/{sample}/ashleys_merged_bam/PLATE_ROW/{row}.bam",
        output:
            "{folder}/{sample}/ashleys_merged_bam/PLATE_ROW/{row}.bam.bai",
        log:
            "{folder}/log/ashleys_merged_bam/{sample}/{row}.log",
        conda:
            "../envs/mc_bioinfo_tools.yaml"
        resources:
            mem_mb=get_mem_mb,
        shell:
            "samtools index {input} > {log} 2>&1"
            
    rule alfred_plate_row:
        input:
            bam="{folder}/{sample}/ashleys_merged_bam/PLATE_ROW/{row}.sort.mdup.bam",
            bam_bai="{folder}/{sample}/ashleys_merged_bam/PLATE_ROW/{row}.sort.mdup.bam.bai",
            fasta=config["references_data"][config["reference"]]["reference_fasta"],
            fasta_index="{fasta}.fai".format(
                fasta=config["references_data"][config["reference"]]["reference_fasta"]
            ),
        output:
            alfred_json="{folder}/{sample}/alfred/{row}.json.gz",
            alfred_tsv="{folder}/{sample}/alfred/{row}.tsv.gz",
        log:
            "{folder}/{sample}/log/alfred/{row}.log",
        resources:
            mem_mb=get_mem_mb_heavy,
        conda:
            "../envs/dev/mc_bioinfo_tools.yaml"
        shell:
            """
            alfred qc -r {input.fasta} -j {output.alfred_json} -o {output.alfred_tsv} {input.bam}
            """