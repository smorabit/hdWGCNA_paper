import pandas as pd
import pyranges as pr
import cerberus
import swan_vis as swan

def make_reports(gname):
    gb='annotation'
    sg.gen_report(gname,
                  prefix='figures/{}'.format(gname),
                  layer='pi',
                  cmap='magma',
                  groupby=gb,
                  transcript_col='tid',
                  metadata_cols=[gb],
                  display_numbers=True,
                  datasets={'annotation': ['RGL', 'NIPCs', 'ASC', 'Ependymal']},
                  browser=True)

    sg.gen_report(gname,
              prefix='figures/{}'.format(),
              layer='tpm',
              cmap='viridis',
              groupby=gb,
              transcript_col='tid',
              metadata_cols=[gb],
              display_numbers=True,
              indicate_novel=True,
              datasets={'annotation': ['RGL', 'NIPCs', 'ASC', 'Ependymal']})

def make_sg(ref, annot, ab, meta, o):
    o = o.split('.')[0]

    # make swangraph
    sg = swan.SwanGraph(sc=True)
    sg.add_annotation(ref)
    sg.add_transcriptome(annot)
    sg.add_abundance(ab)
    sg.add_metadata(meta)

    # add metadata colors
    c_dict, order = get_celltype_colors(sg)
    sg.set_metadata_colors('annotation', c_dict)

    # save
    sg.save_graph(o)

def get_celltype_colors(sg):
    c_dict = {'ASC': '#04a3bd',
              'RGL': '#f0be3d',
              'Ependymal': '#931e18',
              'NIPCs': '#da7901'}
    order = ['RGL', 'NIPCs', 'ASC', 'Ependymal']

    celltypes = sg.adata.obs.annotation.unique().tolist()
    for c in celltypes:
        if c not in c_dict.keys():
            c_dict[c] = 'gray'
            order.append(c)

    return c_dict, order


def get_biotype_map():
    """
    Get a dictionary mapping each gene type to a more general biotype
    """
    map = {'protein_coding': ['protein_coding'],
           'lncRNA': ['lincRNA',
                      'lncRNA',
                      'processed_transcript',
                      'sense_intronic',
                      '3prime_overlapping_ncRNA',
                      'bidirectional_promoter_lncRNA',
                      'sense_overlapping',
                      'non_coding',
                      'macro_lncRNA',
                      'antisense'],
           'pseudogene': ['unprocessed_pseudogene',
                          'translated_unprocessed_pseudogene',
                          'transcribed_unprocessed_pseudogene',
                          'processed_pseudogene',
                          'transcribed_processed_pseudogene',
                          'transcribed_unitary_pseudogene',
                          'unitary_pseudogene',
                          'polymorphic_pseudogene',
                          'pseudogene',
                          'translated_processed_pseudogene'],
           'miRNA': ['miRNA'],
           'other': ['snRNA', 'vault_RNA',
                     'misc_RNA', 'TEC',
                     'snoRNA', 'scaRNA',
                     'rRNA_pseudogene', 'rRNA',
                     'IG_V_pseudogene',
                     'scRNA', 'IG_V_gene',
                     'IG_C_gene', 'IG_J_gene',
                     'sRNA', 'ribozyme',
                     'vaultRNA', 'TR_C_gene',
                     'TR_J_gene', 'TR_V_gene',
                     'TR_V_pseudogene', 'TR_D_gene',
                     'IG_C_pseudogene', 'IG_D_gene',
                     'IG_pseudogene', 'Mt_tRNA',
                     'Mt_rRNA', 'TR_J_pseudogene',
                     'IG_J_pseudogene']}
    return map

def get_gene_info(gtf, o):
    df = pr.read_gtf(gtf, as_df=True)

    # remove sirvs and erccs
    print(len(df.index))
    df = df.loc[(~df.Chromosome.str.contains('SIRV'))&~(df.Chromosome.str.contains('ERCC'))]
    print(len(df.index))

    # only gene
    df = df.loc[df.Feature == 'gene'].copy(deep=True)

    # rename some columns
    m = {'gene_id': 'gid',
         'gene_name': 'gname',
         'transcript_id': 'tid',
         'gene_type': 'biotype'}
    df.rename(m, axis=1, inplace=True)

    map = get_biotype_map()

    beeps = []
    for key, item in map.items():
        beeps += item

    set(df.biotype.unique().tolist())-set(beeps)

    # pivot map
    biotype_map = {}
    for key, biotypes in map.items():
        for biotype in biotypes:
            biotype_map[biotype] = key

    # then add map to df
    df['biotype_category'] = df.biotype.map(biotype_map)

    # gene length
    df['length'] = (df.Start-df.End).abs()

    # and save
    df = df[['gid', 'gname', 'length', 'biotype', 'biotype_category']]
    df.to_csv(o, sep='\t', index=False)


def make_hier_entry(df, how='t'):
    """
    kind {'g','t'}
    """
    agg_dict = {'min_coord': 'min', 'max_coord': 'max'}
    t_df = df.copy(deep=True)
    t_df['min_coord'] = t_df[['Start', 'End']].min(axis=1)
    t_df['max_coord'] = t_df[['Start', 'End']].max(axis=1)
    if how == 't':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id', 'transcript_id', 'transcript_name',
                   'tss_id', 'tes_id',
                   'new_transcript_id', 'original_transcript_id',
                   'original_transcript_name', 'ag1', 'ag2']
        gb_cols = list(set(gb_cols)&(set(t_df.columns)))
    elif how == 'g':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id']

    cols = gb_cols + ['min_coord', 'max_coord']
    t_df = t_df[cols]
    t_df = t_df.groupby(gb_cols).agg(agg_dict).reset_index()
    t_df.rename({'min_coord': 'Start', 'max_coord': 'End'}, axis=1, inplace=True)
    if how == 't':
        t_df['Feature'] = 'transcript'
    elif how == 'g':
        t_df['Feature'] = 'gene'

    return t_df

def make_gtf_from_exons(df):

    # make transcript entries
    t_df = make_hier_entry(df, how='t')

    # make gene entries
    g_df = make_hier_entry(df, how='g')

    # concat everything and sort by gene id, transcript id, feature rank (gene =0, t =1, exon=2), then start coords
    df = pd.concat([df, t_df, g_df])
    df = cerberus.sort_gtf(df)
    return df

def isoquant_struct_to_gtf(iq_struct,
                           gene_meta,
                           ofile):
    df = pd.read_csv(iq_struct, sep='\t', header=None,
                names=['transcript_id', 'original_transcript_id',
                       'tss', 'ic', 'tes'])

    # pull chromosome and strand out from tss
    df[['Chromosome', 'Strand']] = df.tss.str.split('_', expand=True)[[0,3]]

    # get intron chain coordinates
    df['coords'] = df.ic.str.split(';%;')
    df.drop('ic', axis=1, inplace=True)
    df = df.explode('coords', ignore_index=True)

    # remove weird blank entries
    df = df.loc[df.coords != '']

    # get start and stop of each intron
    df[['Start', 'End']] = df.coords.str.split('_', expand=True)[[1,2]]
    df[['Start', 'End']] = df[['Start', 'End']].astype(int)
    df.drop('coords', axis=1, inplace=True)

    # increment / decrement coords to turn to exons
    df.Start = df.Start-1
    df.End = df.End+1

    # melt to get coords individually
    df = df.melt(id_vars=['transcript_id', 'original_transcript_id',
                     'tss', 'tes', 'Chromosome', 'Strand'],
                     value_vars=['Start', 'End'],
                     value_name='Coord')
    df.drop('variable', axis=1, inplace=True)
    df.Coord = df.Coord.astype(int)

    # sort based on strand
    fwd, rev = cerberus.get_stranded_gtf_dfs(df)
    fwd.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, True], inplace=True)
    rev.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, False], inplace=True)

    # get start / end coords of tss / tes regions
    df[['tss_Start', 'tss_End']] = df.tss.str.split('_', expand=True)[[1,2]].astype(int)
    df[['tes_Start', 'tes_End']] = df.tes.str.split('_', expand=True)[[1,2]].astype(int)

    fwd, rev = cerberus.get_stranded_gtf_dfs(df)

    # loop through strands and ends
    ends = pd.DataFrame()
    for mode in ['tss', 'tes']:
        for strand, temp in zip(['+', '-'], [fwd, rev]):
            end_cols = ['{}_Start'.format(mode),
                        '{}_End'.format(mode)]
            keep_cols = ['transcript_id', 'original_transcript_id',
                         'Chromosome', 'Strand']
            keep_cols += end_cols
            temp = temp[keep_cols]
            old_end, new_end, func = cerberus.get_update_ends_settings(strand, mode)
            temp['Coord'] = temp[end_cols].apply(func, axis=1)
            temp.drop(end_cols, axis=1, inplace=True)
            temp.drop_duplicates(inplace=True)
            ends = pd.concat([ends, temp])


    # drop unnecessary columns
    df.drop(['tss', 'tes', 'tss_Start',
             'tes_Start', 'tss_End', 'tes_End'],
            axis=1, inplace=True)

    # concatenate and then sort based on strand
    df = pd.concat([df, ends])
    fwd, rev = cerberus.get_stranded_gtf_dfs(df)
    fwd.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, True], inplace=True)
    rev.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, False], inplace=True)
    df = pd.concat([fwd, rev])

    # add gene id and gene name
    df['gene_id'] = df.transcript_id.str.split('_', expand=True)[0]
    gene_info = pd.read_csv(gene_meta, sep='\t')
    gene_info = gene_info[['gid', 'gname']]
    gene_info.rename({'gid': 'gene_id', 'gname': 'gene_name'}, axis=1, inplace=True)
    df = df.merge(gene_info, how='left', on='gene_id')

    # pivot from individual coordinates to pairs of coordinates for each exon
    fwd, rev = cerberus.get_stranded_gtf_dfs(df)
    df = pd.DataFrame()
    for strand, temp in zip(['+', '-'], [fwd, rev]):
        if strand == '+':
            start_ind = 0
            end_ind = 1
        elif strand == '-':
            start_ind = 1
            end_ind = 0

        exons = pd.DataFrame({'Start':temp['Coord'].iloc[start_ind::2].values,
                              'Start_tid':temp['transcript_id'].iloc[start_ind::2].values,
                              'End':temp['Coord'].iloc[end_ind::2].values,
                              'End_tid':temp['transcript_id'].iloc[end_ind::2].values})
        temp = temp.iloc[::2]
        temp.drop('Coord', axis=1, inplace=True)
        exons.drop(['Start_tid','End_tid'], axis=1, inplace=True)

        temp['Start'] = exons['Start'].tolist()
        temp['End'] = exons['End'].tolist()

        df = pd.concat([df, temp])


    # add exon feature tag
    df['Feature'] = 'exon'

    # add gene and transcript entries
    df = make_gtf_from_exons(df)

    # decrement the start coordinate to convert bed style coordinates
    # to gtf style coordinates
    df['Start'] = df['Start']-1

    # write
    df = pr.PyRanges(df)
    df.to_gtf(ofile)
