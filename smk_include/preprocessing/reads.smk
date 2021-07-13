"""
This Snakemake module handles data preprocessing starting from
(unfiltered) Strand-seq FASTQ files. Because Snakemake inherently
works on the level of file names for determining an execution path,
this module hard-codes a minimal set of assumptions about FASTQ file names:

- a file name can be matched by: ^[0-9A-Za-z_\-\.]+$
- a file name must contain the sample name as component
- a file name has the mate number (1 or 2) as last component
--- because the EMBL Strand-seq data have the uninformative word
--- "sequence" as last component, exactly this one special case is
--- handled as part of this module.
- the two mates of a library are consecutive entries in the list of file names after lex-sorting

If you need to modify this module, DO NOT hard-code file name transformations
in the rules below to match the idiosyncratic naming conventions of
project data. If you do this, it will most likely break this module for
other project data. Move this type of preprocessing to a separate module
and fix file names before using the ASHLEYS QC pipeline. Besides, you would
also most likely violate three mantras of the Zen of Python:

- Sparse is better than dense.
- Readability counts.
- Special cases aren't special enough to break the rules.

"""


localrules: link_strandseq_libraries,
            collect_sample_alignments,
            collect_sample_library_checksums


def _keep_file(file_extensions, file_name):
    """
    This little functions exists b/c
    file extensions will always be a problem
    until the end of time...
    """
    return any(file_name.endswith(e) for e in file_extensions)


def find_file_extension(file_name, file_extensions):
    """
    NB: at this point, this cannot fail
    """
    for fext in file_extensions:
        if file_name.endswith(fext):
            return fext
    raise ValueError(f'Cannot find file extension: {file_name} / {file_extensions}')


def annotate_libraries(mate1, mate2, fastq_ext):

    ext1 = find_file_extension(mate1, fastq_ext)
    ext2 = find_file_extension(mate2, fastq_ext)
    if not ext1 == ext2:
        raise ValueError(f'File extensions of file pair do not match: {mate1} / {mate2} -> {ext1} / {ext2}')

    lib1 = mate1.split(ext1)[0]
    lib2 = mate2.split(ext2)[0]

    # this is an explicit "fix" for EMBL data b/c of the
    # unnecessary "sequence" in the name as last static
    # component...
    if 'sequence' in lib1:
        lib1 = lib1.replace('sequence', '')
        lib2 = lib2.replace('sequence', '')
        lib1 = lib1.strip('_.-')
        lib2 = lib2.strip('_.-')
    
    # heuristic: hope for common sense in file naming...
    lib1 = lib1.rstrip('1')
    lib2 = lib2.rstrip('2')
    lib1 = lib1.strip('_.-')
    lib2 = lib2.strip('_.-')

    if not lib1 == lib2:
        raise ValueError(f'Library IDs of mate pairs do not match: {mate1} / {mate2} -> {lib1} / {lib2}')
    
    lib_id = lib1

    lib1 = lib_id + '_1.fastq.gz'
    lib2 = lib_id + '_2.fastq.gz'

    return lib_id, lib1, lib2


def collect_strandseq_libraries(root_path):

    import os
    import sys
    import collections as col
    import functools as fnt
    import re

    if not os.path.isdir(root_path):
        raise ValueError(f'Specified Strand-seq data path does not exist: {root_path}')

    allowed_chars = re.compile('^[0-9A-Za-z_\-\.]+$')

    LIBS_PER_SAMPLE = int(config.get('libs_per_sample', 96))

    debug = bool(config.get('debug', False))

    sample_ids = []
    library_ids = dict()
    library_link_sources = col.defaultdict(list)

    fastq_extensions = config.get('fastq_extensions', ('.fastq.gz', 'fq.gz', '.txt.gz'))

    keep_file = fnt.partial(
        _keep_file,
        fastq_extensions
    )

    for root, dirs, files in os.walk(os.path.abspath(root_path), followlinks=False):
        sample_libraries = sorted([f for f in files if keep_file(f)])
        if not sample_libraries:
            if debug:
                sys.stderr.write(f'Skipping folder (no Strand-seq libraries): {root}\n')
            continue
        num_libs = len(sample_libraries)
        if num_libs % 2 != 0:
            raise ValueError(f'Number of collected Strand-seq libraries is not even (no mate-pairing): {root} / {num_libs}')
        sample_id = os.path.split(root)[1]
        if allowed_chars.match(sample_id) is None:
            raise ValueError(f'Sample ID contains invalid characters (allowed: [0-9A-Za-z_\-\.]): {sample_id}')
        if not all(sample_id in lib for lib in sample_libraries):
            raise ValueError(f'Not all library file names contain sample ID: {root} / {sample_id}')

        sample_library_ids = []

        for i in range(0, num_libs, 2):
            mate1 = sample_libraries[i]
            if allowed_chars.match(mate1) is None:
                raise ValueError(f'Mate 1 ID contains invalid characters (allowed: [0-9A-Za-z_\-\.]): {mate1}')
            mate2 = sample_libraries[i+1]
            if allowed_chars.match(mate2) is None:
                raise ValueError(f'Mate 2 ID contains invalid characters (allowed: [0-9A-Za-z_\-\.]): {mate2}')
            lib_id, mate1_norm, mate2_norm = annotate_libraries(
                mate1,
                mate2,
                fastq_extensions
            )
            mate1_source = os.path.join(root, mate1)
            mate1_target = os.path.join('input', 'fastq', sample_id, mate1_norm)
            library_link_sources[sample_id].append((mate1_source, mate1_target))

            mate2_source = os.path.join(root, mate2)
            mate2_target = os.path.join('input', 'fastq', sample_id, mate2_norm)
            library_link_sources[sample_id].append((mate2_source, mate2_target))

            sample_library_ids.append(lib_id)
        multi_lib_ids = [n for n, c in col.Counter(sample_library_ids).most_common() if c > 1]
        if multi_lib_ids:
            raise ValueError(f'Library IDs for sample {sample_id} encountered at least twice: {multi_lib_ids}')

        sample_ids.append(sample_id)
        if len(sample_library_ids) % LIBS_PER_SAMPLE != 0:
            raise ValueError(f'Number of libraries for sample {sample_id} is not a multiple of {LIBS_PER_SAMPLE}: {len(sample_library_ids)}')
        library_ids[sample_id] = sorted(sample_library_ids)

    multi_sample_ids = [n for n, c in col.Counter(sample_ids).most_common() if c > 1]
    if multi_sample_ids:
        raise ValueError(f'Sample IDs encountered at least twice: {multi_sample_ids}')

    return sorted(sample_ids), library_ids, library_link_sources


SSEQ_DATA_ROOT = config.get('strandseq_data_root', None)
if SSEQ_DATA_ROOT is None:
    raise RuntimeError('Configuration does not contain path to Strand-seq read data files "strandseq_data_root"')

SSEQ_SAMPLE_IDS, SSEQ_SMP_LIB_MAP, SSEQ_LINK_PAIRS = collect_strandseq_libraries(SSEQ_DATA_ROOT)

CONSTRAINT_SSEQ_SAMPLE_IDS = '(' + '|'.join(SSEQ_SAMPLE_IDS) + ')'


rule link_strandseq_libraries:
    """
    This rule creates a symbolic link for each Strand-seq
    library file underneath the root data path that has
    to be specified in the configuration file.

    The following assumptions are checked:
    - one sample per folder under the root data path
    - each library ID contains the sample ID
    - library IDs of pairs/mates are at consecutive positions after lex-sorting
    """
    output:
        directory('input/fastq/{sample}')
    log:
        'log/input/fastq/{sample}.links.txt'
    wildcard_constraints:
        sample = CONSTRAINT_SSEQ_SAMPLE_IDS
    run:
        import os
        os.makedirs(output[0], exist_ok=True)

        LIBS_PER_SAMPLE = int(config.get('libs_per_sample', 96))

        sample = wildcards.sample
        if not sample in SSEQ_SAMPLE_IDS:
            raise ValueError(f'Sample value {sample} is invalid')

        library_link_pairs = SSEQ_LINK_PAIRS[sample]
        if len(library_link_pairs) != LIBS_PER_SAMPLE * 2:
            raise ValueError(f'Only {len(library_link_pairs)} libraries for symlinking for sample {sample}')

        with open(log[0], 'w') as logfile:
            _ = logfile.write(f'#Linking {len(library_link_pairs)} data source files for sample: {sample}\n')

            for src, trg in library_link_pairs:
                os.symlink(src, trg)
                _ = logfile.write(f'SRC\t{src}\n')
                _ = logfile.write(f'TRG\t{trg}\n')
    ### END OF RUN BLOCK


REFERENCE_GENOME_PATH = config['reference_file']
assert os.path.isfile(REFERENCE_GENOME_PATH)
REFERENCE_GENOME_NAME = config['reference_name']
assert REFERENCE_GENOME_NAME in REFERENCE_GENOME_PATH


rule generate_bwa_index:
    input:
        reference = ancient(REFERENCE_GENOME_PATH)
    output:
        multiext(
            f'input/bwa_index/{REFERENCE_GENOME_NAME}',
            '.amb',
            '.ann',
            '.bwt',
            '.pac',
            '.sa'
        )
    log:
        f'log/input/bwa_index/{REFERENCE_GENOME_NAME}.bwa.log'
    benchmark:
        f'rsrc/input/bwa_index/{REFERENCE_GENOME_NAME}.bwa.rsrc'
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    conda:
         '../../environment/conda_biotools.yml'
    params:
        prefix = os.path.join('input', 'bwa_index', REFERENCE_GENOME_NAME)
    shell:
        'bwa index -p {params.prefix} {input.reference} &> {log}'


rule bwa_strandseq_to_reference_alignment:
    input:
        mate1 = 'input/fastq/{sample}/{library_id}_1.fastq.gz',
        mate2 = 'input/fastq/{sample}/{library_id}_2.fastq.gz',
        ref_index = 'input/bwa_index/{reference}.bwt'
    wildcard_constraints:
        sample = CONSTRAINT_SSEQ_SAMPLE_IDS
    output:
        bam = temp(os.path.join(
            'output',
            'alignments',
            '{reference}',
            'tmp',
            '{sample}/{library_id}.psort.bam'
        ))
    log:
        bwa = os.path.join(
            'log',
            'output',
            'alignments',
            '{reference}',
            '{sample}/{library_id}.bwa.log'
        )
    benchmark:
        os.path.join(
            'rsrc',
            'output',
            'alignments',
            '{reference}',
            '{sample}/{library_id}.bwa.rsrc'
        )
    conda:
        '../../environment/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 12288 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        idx_prefix = lambda wildcards, input: input.ref_index.rsplit('.', 1)[0],
        discard_flag = config['discard_sseq_aln']
    shell:
        'bwa mem -t {threads}'
            ' -R "@RG\\tID:{wildcards.library_id}\\tPL:Illumina\\tSM:{wildcards.sample}"'
            ' -v 2 {params.idx_prefix} {input.mate1} {input.mate2} 2> {log.bwa} | '
            ' samtools view -u -F {params.discard_flag} /dev/stdin | '
            ' samtools sort -@ 2 -m 512M -l 6 -O BAM > {output.bam}'


rule mark_duplicate_reads_strandseq:
    input:
        bam = os.path.join(
            'output',
            'alignments',
            '{reference}',
            'tmp',
            '{sample}/{library_id}.psort.bam'
        )
    output:
        bam = os.path.join(
            'output',
            'alignments',
            '{reference}',
            '{sample}/{library_id}.psort.mdup.bam'
        ),
        bai = os.path.join(
            'output',
            'alignments',
            '{reference}',
            '{sample}/{library_id}.psort.mdup.bam.bai'
        )
    log:
        os.path.join(
            'log',
            'output',
            'alignments',
            '{reference}',
            '{sample}/{library_id}.psort.mdup.log'
        )
    benchmark:
        os.path.join(
            'rsrc',
            'output',
            'alignments',
            '{reference}',
            '{sample}/{library_id}.psort.mdup.rsrc'
        )
    conda:
        '../../environment/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: config['num_cpu_low'] * 512 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'sambamba markdup -t {threads} --overflow-list-size 600000 {input.bam} {output.bam} &> {log}'
            ' && '
        'samtools index {output.bam} &>> {log}'


rule compute_file_checksum:
    input:
        'input/fastq/{sample}/{library_id}_1.fastq.gz',
        'input/fastq/{sample}/{library_id}_2.fastq.gz',
    output:
        'output/checksums/{sample}/{library_id}_1.file-chk.md5',
        'output/checksums/{sample}/{library_id}_2.file-chk.md5'
    wildcard_constraints:
        sample = CONSTRAINT_SSEQ_SAMPLE_IDS
    resources:
        mem_mb = lambda wildcards, attempt: 92 * attempt,
        runtime_hrs = 0,
        runtime_min = 10
    shell:
        'md5sum {input[0]} | cut -d " " -f 1 > {output[0]}'
            ' && '
        'md5sum {input[1]} | cut -d " " -f 1 > {output[1]}'


rule compute_sequence_checksum:
    """
    Compute checksum over concatenated read sequences after sorting
    """
    input:
        'input/fastq/{sample}/{library_id}_{mate}.fastq.gz'
    output:
        'output/checksums/{sample}/{library_id}_{mate}.seq-chk.md5'
    wildcard_constraints:
        sample = CONSTRAINT_SSEQ_SAMPLE_IDS,
        mate = '(1|2)'
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt
    run:
        import gzip
        import hashlib

        reads = []
        with gzip.open(input[0], 'rb') as fastq:
            for ln, line in enumerate(fastq, start=1):
                if ln % 4 == 0:
                    continue
                if ln % 2 == 0:
                    rd = line.strip()
                    assert rd, f'FASTQ line iter failed in line {ln} for file {input[0]}'
                    reads.append(rd)
        
        reads = b''.join(sorted(reads))
        chk = hashlib.md5(reads).hexdigest()
        with open(output[0], 'w') as dump:
            _ = dump.write(chk + '\n')
    # END OF RUN BLOCK


rule collect_sample_library_checksums:
    """
    Convenience rule for housekeeping, i.e.
    collect all library checksums per sample
    """
    input:
        file_chk = lambda wildcards: expand(
            'output/checksums/{{sample}}/{library_id}_{mate}.file-chk.md5',
            library_id=SSEQ_SMP_LIB_MAP[wildcards.sample],
            mate=[1,2]
        ),
        seq_chk = lambda wildcards: expand(
            'output/checksums/{{sample}}/{library_id}_{mate}.seq-chk.md5',
            library_id=SSEQ_SMP_LIB_MAP[wildcards.sample],
            mate=[1,2]
        )
    output:
        'output/checksums/{sample}.libchk.tsv'
    run:
        import os
        import io

        out_buffer = io.StringIO()
        header = '\t'.join([
            'sample',
            'library',
            'file_md5_1',
            'seq_md5_1',
            'file_md5_2',
            'seq_md5_2'
        ])
        _ = out_buffer.write(header + '\n')
        sample = wildcards.sample
        if len(input.file_chk) != len(input.seq_chk):
            raise ValueError(f'Missing checksum files: {len(input.file_chk)} file MD5 vs {len(input.seq_chk)} sequence MD5')

        for fchk, schk in zip(sorted(input.file_chk), sorted(input.seq_chk)):
            fchk_fname = os.path.basename(fchk)
            lib_id = fchk_fname.rsplit('_', 1)[0]
            mate = 1 if '_1.' in fchk_fname else 2
            with open(fchk, 'r') as dump:
                file_md5 = dump.readline().strip()
            with open(schk, 'r') as dump:
                seq_md5 = dump.readline().strip()
            if mate == 1:
                _ = out_buffer.write(f'{sample}\t{lib_id}\t{file_md5}\t{seq_md5}\t')
            else:
                _ = out_buffer.write(f'{file_md5}\t{seq_md5}\n')
        
        with open(output[0], 'w') as tsv:
            _ = tsv.write(out_buffer.getvalue())
    # END OF RUN BLOCK
            
  
rule collect_sample_alignments:
    """
    Convenience rule to trigger alignments
    for just one sample (testing etc.)
    """
    input:
        bam = lambda wildcards: expand(
            'output/alignments/{{reference}}/{{sample}}/{library_id}.psort.mdup.bam',
            library_id=SSEQ_SMP_LIB_MAP[wildcards.sample]
        )
    output:
        'output/alignments/{reference}_{sample}.libraries.txt'
    run:
        if len(input.bam) == 0:
            raise ValueError(f'No alignments collected for {wildcards.reference} and sample {wildcards.sample}')
        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(sorted(input.bam)) + '\n')
