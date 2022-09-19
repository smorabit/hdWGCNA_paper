import swan_vis as swan
import pandas as pd


def make_reports(gname):
    gb='annotation'
    sg.gen_report(gname,
                  prefix='figures/{}'.format(gname.lower()),
                  layer='pi',
                  cmap='magma',
                  groupby=gb,
                  transcript_col='tid',
                  metadata_cols=[gb],
                  display_numbers=True,
                  datasets={'annotation': ['RGL', 'NIPCs', 'ASC', 'Ependymal']},
                  browser=True)

    sg.gen_report(gname,
              prefix='figures/{}'.format(gname.lower()),
              layer='tpm',
              cmap='viridis',
              groupby=gb,
              transcript_col='tid',
              metadata_cols=[gb],
              display_numbers=True,
              indicate_novel=True,
              datasets={'annotation': ['RGL', 'NIPCs', 'ASC', 'Ependymal']})


sg = swan.read('swan.p')

# make reports for all of the genes that sam gave me
df = pd.read_csv('multimod_genes.csv')
gnames = df.gene_name.unique().tolist()
for g in gnames:
  if g in sg.t_df.gname.tolist():
      make_reports(g)
