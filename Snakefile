
SKIP_DATA_PROCESSING = bool(config.get('skip_data_processing', False))

CLUSTER_LOG_FOLDER = config.get('cluster_log_folder', 'log/cluster_jobs')

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
    include: 'smk_include/processing/prediction.smk'

include: 'smk_include/preprocessing/labels_hgsvc.smk'

localrules: create_cluster_log_folders,
            install_source_ashleys_qc,
            prepare_hgsvc_labels,
            run_sseq_checksums,
            run_sseq_alignments,
            run_sample_predictions

###############################################
### Aggregate rules for metadata processing ###
###############################################

rule prepare_hgsvc_labels:
    input:
        os.path.join(CLUSTER_LOG_FOLDER, 'stderr'),
        rules.hgsvc_annotation_format_conversion.output.table


if not SKIP_DATA_PROCESSING:

##################################################
### Aggregate rules for data (pre-) processing ###
##################################################

    rule link_sseq_input:
        input:
            os.path.join(CLUSTER_LOG_FOLDER, 'stderr'),
            expand(
                'input/fastq/{sample}',
                sample=SSEQ_SAMPLE_IDS
            )


    rule run_sseq_checksums:
        input:
            os.path.join(CLUSTER_LOG_FOLDER, 'stderr'),
            expand(
                'output/checksums/{sample}.libchk.tsv',
                sample=SSEQ_SAMPLE_IDS
            )


    rule run_sseq_alignments:
        input:
            os.path.join(CLUSTER_LOG_FOLDER, 'stderr'),
            expand(
                'output/alignments/{sample}_{reference}.libraries.txt',
                reference=REFERENCE_GENOME_NAME,  # set in preprocessing::reads.smk
                sample=SSEQ_SAMPLE_IDS
            )

    rule run_sample_predictions:
        input:
            os.path.join(CLUSTER_LOG_FOLDER, 'stderr'),
            expand(
                'output/library_predictions/{sample}_{reference}.predictions.tsv',
                reference=REFERENCE_GENOME_NAME,  # set in preprocessing::reads.smk
                sample=SSEQ_SAMPLE_IDS
            )


#############################################
### Auxiliary rules for environment setup ###
#############################################

rule create_cluster_log_folders:
    """
    Log paths for stdout/stderr of cluster scripts
    are configured for HHU-HILBERT. You may need/want
    to change this to fit your local setup.
    """
    output:
        stdout = directory(os.path.join(CLUSTER_LOG_FOLDER, 'stdout'),),
        stderr = directory(os.path.join(CLUSTER_LOG_FOLDER, 'stderr'),)
    priority: 1000
    shell:
        'mkdir -p {output.stdout} && mkdir -p {output.stderr}'


rule install_source_ashleys_qc:
    """
    This rule only needs to exist until
    ASHLEYS is available via bioconda
    (or to force a source build)
    """
    output:
        chk = 'repositories/ashleys-qc.setup.ok'
    log:
       'log/repositories/ashleys-qc.setup.log'
    conda:
        'environment/conda_ashleys.yml'
    params:
        repo_folder = lambda wildcards, output: output.chk.split('/')[0]
    shell:
         '( cd {params.repo_folder} && '
         'git clone https://github.com/friendsofstrandseq/ashleys-qc.git && '
         'cd ashleys-qc && '
         'git checkout develop && '
         'python setup.py install && cd ../../ && touch {output} ; ) > {log} 2>&1'
