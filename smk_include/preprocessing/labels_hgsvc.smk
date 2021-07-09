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

        # for consistency, map comments to a
        # "label" similar to the binary quality labels
        # (see rule: hgsvc_assign_binary_quality_labels)
        comment_to_label = {
            'OK': 'label_ok',
            'good': 'label_good',
            'high_Q': 'label_good',
            'over_BrdU': 'label_overBrdU',
            'under_BrdU': 'label_underBrdU',
            'low_Q': 'label_lowQ',
            'low_cov': 'label_lowCov',
            'wgs': 'label_wgs'
        }

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
        
        # label libraries on the basis of unified comments
        for comment in df['comments_unified'].unique():
            label = comment_to_label[comment]
            if label not in df:
                df[label] = 0
            df.loc[df['comments_unified'] == comment, label] = 1

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


def hgsvc_set_column_order(columns):

    sort_order = [
        'sample',
        'cell',
        'sts_score',
        'final',
        'comments_unified',
        'comments_ads',
        'label_highQ',
        'label_medQ',
        'label_clustQ',
        'label_badQ',
        'label_good',
        'label_ok',
        'label_overBrdU',
        'label_underBrdU',
        'label_lowQ',
        'label_lowCov',
        'label_wgs'
    ]
    lex_sort = sorted(c for c in columns if c not in sort_order and not c.startswith('nb_'))
    sort_order.extend(lex_sort)
    sort_order.extend(['nb_p', 'nb_r', 'nb_a'])
    assert all(c in columns for c in sort_order), 'Typo in sort_order entry'
    assert all(c in sort_order for c in columns), 'Missing original column from sort_order'
    return sort_order


rule hgsvc_assign_binary_quality_labels:
    input:
        norm_annotation = 'metadata/annotation/hgsvc_sample_library_annotation.tsv'
    output:
        labeled_annotation = 'metadata/annotation/hgsvc_annotation_labeled.tsv'
    run:
        import pandas as pd

        # lower threshold for number of "good" reads
        # for low-quality libraries to be considered
        # for clustering applications (e.g. PGAS)
        rc_threshold = 5e4
        # for unified categories "low_Q" and "low_cov",
        # filter column "comments_ads" to remove libraries
        # even unwanted for clustering
        bad_words = ['high_bg', 'chr', 'loss', 'gain', 'patchy', 'gap', 'ugly']
        remove_library = lambda x: any(word in x for word in bad_words)

        md = pd.read_csv(input.norm_annotation, sep='\t', header=0)
        md['label_highQ'] = 0
        md['label_medQ'] = 0
        md['label_clustQ'] = 0
        md['label_badQ'] = 0

        md.loc[md['sts_score'] == 1, 'label_highQ'] = 1
        md.loc[(md['sts_score'] == 0.5) & ((md['final'] == 1) | (md['final'] == 0)), 'label_medQ'] = 1
        md.loc[(md['sts_score'] == 0) & (md['good'] < rc_threshold), 'label_badQ'] = 1

        # 0 - label all control probes (everything with wgs characteristics)
        selector = md['comments_unified'] == 'wgs'
        md.loc[selector, 'label_badQ'] = 1

        # what is left now: lower quality libraries with acceptable number of good reads
        # 1 - under_BrdU
        selector = (md['sts_score'] == 0) & (md['good'] >= rc_threshold) & (md['comments_unified'] == 'under_BrdU')
        md.loc[selector, 'label_badQ'] = 1

        # 2 - over_BrdU
        selector = (md['sts_score'] == 0) & (md['good'] >= rc_threshold) & (md['comments_unified'] == 'over_BrdU')
        md.loc[selector, 'label_clustQ'] = 1

        # 3 - filter categories low_cov and low_Q
        selector = (md['sts_score'] == 0) & (md['good'] >= rc_threshold) & (md['comments_unified'].isin(['low_cov', 'low_Q']))
        # for simplicity, scan all comments to have compatible shape
        bad_comments = md['comments_ads'].apply(remove_library)

        select_bad = selector & bad_comments
        # sanity check - don't process values twice
        assert (md.loc[select_bad, 'label_badQ'] == 0).all(), 'Select bad: value is already set'
        md.loc[select_bad, 'label_badQ'] = 1

        select_good = selector & ~bad_comments
        assert (md.loc[select_good, 'label_clustQ'] == 0).all(), 'Select good/clustQ: value is already set'
        md.loc[select_good, 'label_clustQ'] = 1

        # check that all libraries have a label assigned
        total_assigned_labels = md[['label_highQ', 'label_medQ', 'label_clustQ', 'label_badQ']].sum(axis=1).sum()
        if total_assigned_labels != md.shape[0]:
            unassigned = md[['label_highQ', 'label_medQ', 'label_clustQ', 'label_badQ']].sum(axis=1) == 0
            print(md.loc[unassigned, ])
            raise ValueError(f'Assigned only {total_assigned_labels} for a total of {md.shape[0]} libraries.')
        
        # now complete labeling by logic
        select_highq = md['label_highQ'] == 1
        md.loc[select_highq, ['label_medQ', 'label_clustQ']] = 1, 1
        select_medq = md['label_medQ'] == 1
        md.loc[select_medq, 'label_clustQ'] = 1

        reordered = hgsvc_set_column_order(md.columns)
        md = md[reordered]

        md.to_csv(
            output.labeled_annotation,
            sep='\t',
            header=True,
            index=False
        )
    # END OF RUN BLOCK
