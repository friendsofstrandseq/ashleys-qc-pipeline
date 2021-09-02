
rule compute_strandseq_library_features:
    input:
        setup_ok = 'repositories/ashleys-qc.setup.ok',
        bam = lambda wildcards: expand(
            'output/alignments/{{reference}}/{{sample}}/{library_id}.psort.mdup.bam',
            library_id=SSEQ_SMP_LIB_MAP[wildcards.sample]
        ),
        bai = lambda wildcards: expand(
            'output/alignments/{{reference}}/{{sample}}/{library_id}.psort.mdup.bam.bai',
            library_id=SSEQ_SMP_LIB_MAP[wildcards.sample]
        ),
    output:
        table = 'output/feature_tables/{sample}_{reference}.features.tsv'
    log:
        'log/output/feature_tables/{sample}_{reference}.features.log'
    benchmark:
        'rsrc/output/feature_tables/{sample}_{reference}.features' + '.t{}.rsrc'.format(config['num_cpu_medium'])
    conda:
        '../../environment/conda_ashleys.yml'
    threads: config['num_cpu_medium']
    params:
        bam_folder = lambda wildcards, input: os.path.dirname(input.bam[0]),
        bam_ext = '.psort.mdup.bam',
        windows = '5000000 2000000 1000000 800000 600000 400000 200000',
    shell:
        'ashleys.py -j {threads} --verbose features -f {params.bam_folder} '
        ' -w {params.windows} -e {params.bam_ext} -o {output.table} &> {log}'


rule predict_strandseq_library_quality:
    input:
        setup_ok = 'repositories/ashleys-qc.setup.ok',
        features = 'output/feature_tables/{sample}_{reference}.features.tsv'
    output:
        table = 'output/library_predictions/{sample}_{reference}.{model}.predictions.tsv',
    conda:
        '../../environment/conda_ashleys.yml'
    params:
        model_file = lambda wildcards: os.path.abspath(f'repositories/ashleys-qc/models/{wildcards.model}.pkl'))
    shell:
        'ashleys.py predict -p {input.features} -o {output.table} -m {params.model_file}'
