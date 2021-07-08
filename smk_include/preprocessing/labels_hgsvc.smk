import pathlib
import re

localrules: hgsvc_annotation_format_conversion


def find_annotation_path():
    init_dir = workflow.basedir
    current_dir = pathlib.Path(init_dir)
    while 1:
        p = pathlib.Path(current_dir, 'annotation')
        if p.is_dir():
            return p
        else:
            current_dir = current_dir.parent
            if current_dir.is_mount():
                break
    raise RuntimeError(f'Could not determine annotation folder starting at {init_dir}')


HGSVC_IDENTIFIER_REGEXP = re.compile('^[A-Z0-9]+$')


def normalize_hgsvc_sample_name(sample_name):
    init_name = sample_name
    if HGSVC_IDENTIFIER_REGEXP.match(sample_name) is not None:
        return sample_name
    sample_name = sample_name.split('_')[1]
    if HGSVC_IDENTIFIER_REGEXP.match(sample_name) is not None:
       return sample_name
    sample_name = sample_name.split('x')[0]
    if HGSVC_IDENTIFIER_REGEXP.match(sample_name) is not None:
       return sample_name
    raise ValueError(f'Cannot normalize HGSVC sample name: {init_name}')


def normalize_hgsvc_cell_name(cell_name):
    init_name = cell_name
    if HGSVC_IDENTIFIER_REGEXP.match(cell_name) is not None:
        return cell_name
    cell_name = cell_name.split('x')[1]
    if HGSVC_IDENTIFIER_REGEXP.match(cell_name) is not None:
       return cell_name
    raise ValueError(f'Cannot normalize HGSVC cell name: {init_name}')


ANNOTATION_FOLDER = find_annotation_path()


rule hgsvc_annotation_format_conversion:
    input:
        txt_tables = [
            pathlib.Path(ANNOTATION_FOLDER, 'hgsvc', 'LibSelection_Cell_scoring_Information_ReplacementSamples_July2020.txt'),
            pathlib.Path(ANNOTATION_FOLDER, 'hgsvc', '20200528_ASanders_QCselect_HGSVClibs.tsv'),
        ]
    output:
        table = 'metadata/annotation/hgsvc_sample_library_annotation.tsv',
        comments = 'metadata/supplementary/hgsvc_annotation_comments.counts.tsv',
        samples = 'metadata/supplementary/hgsvc_annotation_samples.counts.tsv',
    run:
        import pandas as pd

        concat = []

        for txt_table in input.txt_tables:
            df = pd.read_csv(
                txt_table,
                sep='\t',
                header=0,
                comment='#',
                index_col=None
            )
            # the BAM column is unnecessary
            try:
                df.drop('bam', axis=1, inplace=True)
            except KeyError:
                print(df.columns)
                raise
            df['sample'] = df['sample'].apply(normalize_hgsvc_sample_name)
            df['cell'] = df['cell'].apply(normalize_hgsvc_cell_name)
            concat.append(df)

        df = pd.concat(concat, axis=0, ignore_index=False)
        df.to_csv(
            output.table,
            sep='\t',
            header=True,
            index=False
        )
        sample_counts = df['sample'].value_counts()
        with open(output.samples, 'w') as dump:
            _ = dump.write('\n'.join([f'{label}\t{count}' for label, count in sample_counts.items()]) + '\n')

        comment_counts = df['comments_unified'].value_counts()
        with open(output.comments, 'w') as dump:
            _ = dump.write('\n'.join([f'{comment}\t{count}' for comment, count in comment_counts.items()]) + '\n')


