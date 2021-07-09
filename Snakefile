
SKIP_DATA_PROCESSING = bool(config.get('skip_data_processing', False))

if not SKIP_DATA_PROCESSING:
    # For debugging and testing purposes,
    # setting the config value "skip_data_processing=true"
    # prevents including any modules that contain
    # data processing tasks (checksum computation,
    # alignment, model training etc.).
    # This can be used to locally (=laptop) develop metadata
    # preprocessing rules (e.g., deriving library quality labels)
    # before deploying the pipeline on a server or cluster.
    include: 'smk_include/preprocessing/reads.smk'

include: 'smk_include/preprocessing/labels_hgsvc.smk'

localrules: create_cluster_log_folders,
            prepare_hgsvc_labels,
            run_sseq_checksums,
            run_sseq_alignments

###############################################
### Aggregate rules for metadata processing ###
###############################################

rule prepare_hgsvc_labels:
    input:
        'log/cluster_jobs/stderr',
        rules.hgsvc_annotation_format_conversion.output.table


if not SKIP_DATA_PROCESSING:

##################################################
### Aggregate rules for data (pre-) processing ###
##################################################

    rule run_sseq_checksums:
        input:
            'log/cluster_jobs/stderr',
            expand(
                'output/checksums/{sample}.libchk.tsv',
                sample=SSEQ_SAMPLE_IDS
            )


    rule run_sseq_alignments:
        input:
            'log/cluster_jobs/stderr',
            expand(
                'output/alignments/{reference}_{sample}.libraries.txt',
                reference=REFERENCE_GENOME_NAME,  # set in preprocessing::reads.smk
                sample=SSEQ_SAMPLE_IDS
            )


#############################################
### Auxiliary rules for environment setup ###
#############################################

rule create_cluster_log_folders:
    """
    This always needs to be executed in a cluster environment,
    otherwise cluster job logs may not be saved.
    ! Not important in non-cluster settings !
    """
    output:
        stdout = directory('log/cluster_jobs/stdout'),
        stderr = directory('log/cluster_jobs/stderr')
    shell:
        'mkdir -p {output.stdout} && mkdir -p {output.stderr}'
