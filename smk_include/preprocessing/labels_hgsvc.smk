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
HGSVC_SAMPLE_ALIAS_REGEXP = re.compile('^[A-Z]+[0-9]+')

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


def set_sample_alias(sample_name):
    init_sample = sample_name
    sample_name = re.search('^[A-Z]+[0-9]+', sample_name).group(0)
    if 'GM' in sample_name:
        sample_name = sample_name.replace('GM', 'NA')
    assert HGSVC_IDENTIFIER_REGEXP.match(sample_name) is not None, \
        f'Cannot set sample alias for sample {init_sample}'
    return sample_name


ANNOTATION_FOLDER = find_annotation_path()


def add_rank_information(df):
    import pandas as pd
    import numpy as np
    bins = np.arange(0, 1.01, 0.1, dtype=np.float16).round(1)
    rank_sets = []
    for (sample, sample_alias), lib_infos in df.groupby(['sample', 'sample_alias']):
        mapped_rnk = lib_infos['mapped'].rank(ascending=True, pct=True).round(3)
        mapped_rnk_bin = pd.cut(mapped_rnk, bins=bins, right=True, retbins=False, labels=False)
        good_rnk = lib_infos['good'].rank(ascending=True, pct=True).round(3)
        good_rnk_bin = pd.cut(good_rnk, bins=bins, right=True, retbins=False, labels=False)
        subset = pd.DataFrame(
            {
                'mapped_rank_pct': mapped_rnk,
                'mapped_rank_bin': mapped_rnk_bin,
                'good_rank_pct': good_rnk,
                'good_rank_bin': good_rnk_bin
            },
            index=lib_infos.index
        )
        rank_sets.append(subset)
    rank_sets = pd.concat(rank_sets, axis=0, ignore_index=False)
    df = df.join(rank_sets, how='inner')
    return df


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
            df.drop('bam', axis=1, inplace=True)
            df['sample'] = df['sample'].apply(normalize_hgsvc_sample_name)
            df['sample_alias'] = df['sample'].apply(set_sample_alias)
            df['cell'] = df['cell'].apply(normalize_hgsvc_cell_name)
            df.columns = [c.lower() for c in df.columns]
            concat.append(df)

        df = pd.concat(concat, axis=0, ignore_index=True)
        df = add_rank_information(df)
        df.sort_values(['sample', 'cell'], ascending=True, inplace=True)

        if pd.isna(df).any(axis=0).any():
            raise ValueError(f'N/A introduced during format conversion')

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
    # END OF RUN BLOCK

