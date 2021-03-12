import os as os

configfile: 'config_snakemake.yaml'
INPUT_PATH = config['INPUT_PATH']
SRC_PATH = config['SRC_PATH']
OUTPUT_PATH = config['OUTPUT_PATH']
threads = config['cores']


def collect_input_data(collect_path):
    strand_seq_files = set()
    for f in os.listdir(collect_path):
        if os.path.isfile(os.path.join(collect_path, f)):
            # remove _1.fastq.gz endings
            strand_seq_files.add(f)#.rsplit('_', 1)[0])
            continue

    return strand_seq_files

def collect_input_bam(collect_path):
    strand_seq_files = set()
    for f in os.listdir(INPUT_PATH):
        if os.path.isfile(os.path.join(INPUT_PATH, f)):
            # remove .bam endings
            strand_seq_files.add(f.split('.', 1)[0])
            continue

    return strand_seq_files

all_files = collect_input_bam(INPUT_PATH)# + 'lansdorp/') # SRC_PATH
print(all_files)
window_list = '50_20_10_8_6_4_2' #[5000000, 3000000, 2000000, 1000000, 800000, 600000, 400000, 200000]

rule all:
    input:
           #expand(INPUT_PATH + "{sample_name}.bam", sample_name=all_files),
           #expand(INPUT_PATH + "{sample_name}.sort.bam", sample_name=all_files),
           #expand(INPUT_PATH + "{sample_name}.sort.mdup.bam", sample_name=all_files),
	   #expand(INPUT_PATH + "{sample_name}.sort.mdup.bam.bai", sample_name=all_files),
           expand(OUTPUT_PATH + '{feature_folder}/k562_features.tsv', feature_folder=window_list),
	   #expand(OUTPUT_PATH + '{feature_folder}/{model}/hgsvc_{quality}.tsv', feature_folder=window_list, model=['gb', 'svc'], quality=['high-qual', 'ok-qual', 'clust-qual']),
           #expand(OUTPUT_PATH + '{feature_folder}/prediction_probabilities.tsv', feature_folder=window_list),
           #expand('ashleys-qc/ashleys_install_success.txt')


rule bwa_strandseq_to_reference_alignment:
    input:
        mate1 = SRC_PATH + '{sample_name}_1.fastq.gz',
        mate2 = SRC_PATH + '{sample_name}_2.fastq.gz',
        ref_index = config['reference']
    output:
        bam = INPUT_PATH + '{sample_name}.bam'
    log:
        bwa = INPUT_PATH + 'log/' + '{sample_name}.bwa.log',
        samtools = INPUT_PATH + 'log/' + '{sample_name}.samtools.log',
    benchmark:
        os.path.join(INPUT_PATH + 'log/' + '{sample_name}' + '.t{}.rsrc'.format(threads))
    threads: 6
    params:
        idx_prefix = lambda wildcards, input: input.ref_index.rsplit('.', 1)[0],
        sample = lambda wildcards, input: input.mate1.rsplit('_', 5)[1][5:],
        cell = lambda wildcards, input: input.mate1.rsplit('_', 5)[4]
    shell:
        'bwa mem -t {threads}' # http://bio-bwa.sourceforge.net/bwa.shtml
            ' -R "@RG\\tID:{params.cell}\\tPL:Illumina\\tSM:{params.sample}"'
            ' -v 2 {params.idx_prefix} {input.mate1} {input.mate2} 2> {log.bwa} | '
            ' samtools view -b /dev/stdin > {output.bam} 2> {log.samtools}'


rule samtools_sort_bam:
    input:
        INPUT_PATH + '{sample_name}.bam'
    output:
        INPUT_PATH + '{sample_name}.sort.bam'
    shell:
        'samtools sort -O BAM -o {output} {input}'


rule mark_duplicates:
    input:
        INPUT_PATH + '{sample_name}.sort.bam'
    output:
        INPUT_PATH + '{sample_name}.sort.mdup.bam'
    shell:
        'sambamba markdup {input} {output}'


rule create_bai_files:
    input:
        INPUT_PATH + '{sample_name}.sort.mdup.bam'
    output:
        INPUT_PATH + '{sample_name}.sort.mdup.bam.bai'
    shell:
        'samtools index {input} {output}'


rule train_model:
    input:
        features = OUTPUT_PATH + '{feature_folder}/features.tsv',
        annotation = INPUT_PATH + 'annotation/hgsvc/{quality}.txt'
    output:
        OUTPUT_PATH + '{feature_folder}/{model}/hgsvc_{quality}.tsv'
    conda:
        'environment/ashleys_env.yml'
    params:
        model = '--{model}',
        jobs = 20,
        json = 'ashleys-qc/models/dict_{model}.json'
    shell:
        './ashleys-qc/bin/ashleys.py -j {params.jobs} train -p {input.features} -a {input.annotation} -o {output} {params.model} -js {params.json}'


rule generate_features:
    input:
        ashleys = 'ashleys-qc/ashleys_install_success.txt',
        #bam =  [(INPUT_PATH + "{sample_name}".format(sample_name=sample_name)) for sample_name in all_files],
        path = INPUT_PATH
    output:
        OUTPUT_PATH + '{feature_folder}/k562_features.tsv'
    conda:
        'environment/ashleys_env.yml' # run snakemake --use-conda
    params:
        windows = '5000000 2000000 1000000 800000 600000 400000 200000',
        #windows = '{feature_folder}',
        #windows = '2000000 1000000 800000 600000 400000',
        jobs = 23,
        extension = '.sort.mdup.bam'
    shell:
        './ashleys-qc/bin/ashleys.py -j {params.jobs} features -f {input.path} -w {params.windows} -o {output} --recursive_collect -e {params.extension}'


rule predict:
    input:
        ashleys = 'ashleys-qc/ashleys_install_success.txt',
        path = OUTPUT_PATH + '{feature_folder}/features.tsv'
    output:
        OUTPUT_PATH + '{feature_folder}/prediction_probabilities.tsv'
    conda:
        'environment/ashleys_env.yml'
    params:
        model = 'ashleys-qc/models/svc_default.pkl'
    shell:
        './ashleys-qc/bin/ashleys.py predict -p {input.path} -o {output} -m {params.model}'


rule install_ashleys:
    input:
        'ashleys_install.txt'
    output:
        touch('ashleys-qc/ashleys_install_success.txt')
    log:
        'install_ashleys.log'
    shell:
        '( git clone https://github.com/friendsofstrandseq/ashleys-qc.git &&'
        'cd ashleys-qc &&'
        'python setup.py install ) &> {log}'
