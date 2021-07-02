
include: 'smk_include/preprocessing/reads.smk'
include: 'smk_include/preprocessing/labeling.smk'

localrules: create_cluster_log_folders,
            run_sseq_checksums,
            run_sseq_alignments


rule create_cluster_log_folders:
    """
    This always needs to be executed,
    otherwise cluster job logs may be lost
    """
    output:
        stdout = directory('log/cluster_jobs/stdout'),
        stderr = directory('log/cluster_jobs/stderr')
    shell:
        'mkdir -p {output.stdout} && mkdir -p {output.stderr}'


rule run_sseq_checksums:
    input:
        rules.create_cluster_log_folders.output.stderr,
        expand(
            'output/checksums/{sample}.libchk.tsv',
            sample=SSEQ_SAMPLE_IDS
        )


rule run_sseq_alignments:
    input:
        rules.create_cluster_log_folders.output.stderr,
        expand(
            'output/alignments/{reference}_{sample}.libraries.txt',
            reference=REFERENCE_GENOME_NAME,  # set in preprocessing::reads.smk
            sample=SSEQ_SAMPLE_IDS
        )
