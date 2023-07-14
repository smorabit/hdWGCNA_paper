import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
import anndata
import os
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy import io

os.chdir("~/swaruplab/smorabit/analysis/scWGCNA/t1d")

data_dir = 'data/'
fig_dir = 'figures/'

sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor = 'white', figsize=(6,6), format='pdf')


################################################################################
# load the Parse Bio 1M cell dataset:
################################################################################

adata = sc.read_mtx(data_dir + 'DGE_1M_PBMC.mtx')


# reading in gene and cell data
gene_data = pd.read_csv(data_dir + 'all_genes_1M_PBMC.csv')
cell_meta = pd.read_csv(data_dir + 'cell_metadata_1M_PBMC.csv')

# find genes with nan values and filter
gene_data = gene_data[gene_data.gene_name.notnull()]
gene_data.index = gene_data['gene_name']
notNa = gene_data.index
notNa = notNa.to_list()

# remove genes with nan values and assign gene names
#adata = adata[:,notNa]
adata.var = gene_data
adata.var.set_index('gene_name', inplace=True)
adata.var.index.name = None
adata.var_names_make_unique()

# add cell meta data to anndata object
adata.obs = cell_meta
adata.obs.set_index('bc_wells', inplace=True)
adata.obs.index.name = None
adata.obs_names_make_unique()

sc.pp.filter_cells(adata, min_counts=100)
sc.pp.filter_genes(adata, min_cells=5)
adata.shape

# add the disease status as a column:
adata.obs['Condition'] = adata.obs['sample'].str.split('_', n=1, expand=True).iloc[:,0]

# save the unprocessed dataset:
adata.write_h5ad("{}t1d_pb_1M_unprocessed.h5ad".format(data_dir))


################################################################################
# QC
################################################################################

adata = sc.read_h5ad("{}t1d_pb_1M_unprocessed.h5ad".format(data_dir))

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter out
adata[(adata.obs.pct_counts_mt < 25) &
(adata.obs.n_genes_by_counts < 5000) &
(adata.obs.total_counts < 25000),:].shape

# Do the filtering
adata = adata[(adata.obs.pct_counts_mt < 25) &
(adata.obs.n_genes_by_counts < 5000) &
(adata.obs.total_counts < 25000),:]

################################################################################
# Processing
################################################################################

# preserve the counts as its own layer
adata.layers['counts'] = adata.X

# normalize data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# call HVGs
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25)

# subset by HVGs
adata.raw = adata
adata = adata[:,adata.var.highly_variable].copy()

sc.pp.scale(adata, max_value=10)


# run PCA and Harmony:

sc.tl.pca(adata, svd_solver='arpack')
sc.external.pp.harmony_integrate(adata, key = 'sample')

# neighbors
sc.pp.neighbors(adata, use_rep = 'X_pca_harmony', n_neighbors=20, n_pcs=30, metric='cosine')
sc.tl.umap(adata, min_dist=0.35, method='umap')

sc.tl.leiden(adata, resolution=1)

#compute dendrogram for the dotplot:
sc.tl.dendrogram(adata, 'leiden')

# save the processed data
adata.write_h5ad("{}t1d_pb_1M_processed2.h5ad".format(data_dir))



# plot the umap
sc.pl.umap(
    adata,
    color=['leiden'],
    frameon=False,
    legend_loc='on data',
    legend_fontoutline=1,
    legend_fontsize=9,
    title='',
    save='_leiden2.pdf'
)

# how many cells of each type?
adata.obs.leiden.value_counts()



marker_dict = {
'ASDC_mDC': ['AXL', 'LILRA4', 'SCN9A', 'CLEC4C'],
'ASDC_pDC':['SCT', 'PROC', 'LTK', 'SCN9A'],
'cDC': ['CD14', 'FCER1A', 'CLEC10A', 'ENHO', 'CD1C'],
'Eryth': ['HBM', 'HBD', 'SNCA', 'ALAS2'],
'Mature B': ['MS4A1', 'IGKC', 'IGHM', 'CD24', 'CD22'],
'Memory B': ['BANK1', 'CD79A', 'SSPN'],
'Naive B': ['TCL1A', 'IGHD', 'IL4R', 'CD79B'],
'CD14 Mono': ['CD14', 'LYZ', 'CTSD', 'VCAN', 'CTSS'],
'CD16 Mono': ['LST1', 'YBX1', 'AIF1', 'MS4A7'],
'CD4 CTL': ['GZMH', 'CD4', 'GNLY', 'IL7R', 'CCL5'],
'CD4+ Naive T': ['TCF7', 'LEF1', 'NUCB2', 'LDHB'],
'Proliferating T': ['PCLAF', 'CD8B', 'CD3D', 'TRAC', 'CLSPN'],
'CD8 TCM': ['CD8A', 'SELL', 'CCL5', 'ITGB1'],
'gDT': ['TRDC', 'TRGV9', 'TRDV2', 'KLRD1'],
'MAIT': ['KLRB1', 'NKG7', 'GZMK', 'SLC4A10', 'NCR3', 'CTSW', 'IL7R', 'KLRG1', 'CEBPD', 'DUSP2'],
'NK': ['STMN1', 'KLRC2', 'GNLY', 'S100A4', 'CD3E', 'CST7', 'NKG7'],
'Plasmablast': ['TYMS', 'TK1', 'ASPM', 'SHCBP1', 'TPX2'],
'Platelet': ['GNG11', 'PPBP', 'NRGN', 'PF4', 'TUBB1', 'CLU'],
'Treg': ['B2M', 'FOXP3', 'RTKN2', 'TIGIT', 'CTLA4', 'TRAC'],
'Basophil': ['KIT', 'CPA3', 'TPSAB1', 'CD44', 'GATA2', 'GRAP2']
}



sc.pl.dotplot(
    adata, marker_dict, 'leiden',
    dendrogram=True,
    standard_scale='var', swap_axes=True,
    save = '_leiden_markers_dict2.pdf'
)


sc.pl.violin(adata,[ 'n_counts'], inner='box', size=0,  groupby='leiden', multi_panel=False, rotation=30,  save='_umis.pdf')
sc.pl.violin(adata,[ 'n_genes_by_counts'], inner='box', size=0,  groupby='leiden', multi_panel=False, rotation=30, save='_genes.pdf')
sc.pl.violin(adata,[ 'pct_counts_mt'], inner='box', size=0,  groupby='leiden', multi_panel=False, rotation=30, save='_mt.pdf')





# load cluster anno df
anno_df = pd.read_csv('data/t1d_leiden.csv')

anno_df.leiden = anno_df.leiden.astype(str)
temp = adata.obs.merge(
    anno_df, how='left',
    on = 'leiden')

adata.obs['annotation'] = temp['annotation'].astype(str).to_list()


clusters_to_remove = ['Unknown', 'Doublet']
adata = adata[~adata.obs.annotation.isin(clusters_to_remove)].copy()

# run umap again without the doublet

sc.tl.umap(adata, min_dist=0.35, method='umap')


# plot the umap
sc.pl.umap(
    adata,
    color=['annotation'],
    frameon=False,
    legend_loc='on data',
    legend_fontoutline=1,
    legend_fontsize=9,
    title='',
    save='_annotation.pdf'
)


# write anndata:
adata.write_h5ad('{}/t1d_processed.h5ad'.format(data_dir))

# write obs table
adata.obs['UMAP_1'] = adata.obsm['X_umap'][:,0]
adata.obs['UMAP_2'] = adata.obsm['X_umap'][:,1]
adata.obs.to_csv('{}/t1d_processed_obs.csv'.format(data_dir))

# write var table:
adata.raw.var.to_csv('{}/t1d_processed_var.csv'.format(data_dir))

# save the sparse matrix for Seurat:
X = adata.raw.X
X = scipy.sparse.csr_matrix.transpose(X)
io.mmwrite('{}/t1d_counts.mtx'.format(data_dir), X)

# save the Harmony Matrix for Seurat:
pd.DataFrame(adata.obsm['X_pca_harmony']).to_csv('{}t1d_processed_harmony.csv'.format(data_dir))
