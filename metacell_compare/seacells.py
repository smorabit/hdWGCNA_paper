
# conda activate seacells
# have to type "python3.8" for it to actually get the right python version for this
# following this tutorial: https://github.com/dpeerlab/SEACells/blob/main/notebooks/SEACell_computation.ipynb

import numpy as np
import pandas as pd
import scanpy as sc
import SEACells

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import scipy
from scipy import io

adata = sc.read_h5ad('data/Swarup_2021_processed_counts.h5ad')

# add harmony to adata
h = pd.read_csv('data/Swarup_2021_processed_harmony.csv', index_col=0)
h.index = adata.obs.index
adata.obsm['X_harmony'] = h.to_numpy()


# adata = sc.read_h5ad('data/cd34_multiome_rna.h5ad')


raw_ad = sc.AnnData(adata.X)
raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
adata.raw = raw_ad

# standard pre-processing
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1500)

#sc.tl.pca(adata, n_comps=50, use_highly_variable=True)


##################################################################################
# Running SEACells
##################################################################################

# they recommend one metacell for every 75 real cells
n_SEACells = int(np.floor(adata.obs.shape[0] / 75))

#build_kernel_on = 'X_pca' # key in ad.obsm to use for computing metacells
build_kernel_on = 'X_harmony' # key in ad.obsm to use for computing metacells
                          # This would be replaced by 'X_svd' for ATAC data

## Additional parameters
n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells
waypoint_proportion = 0.9 # Proportion of metacells to initialize using waypoint analysis,
                        # the remainder of cells are selected by greedy selection


model = SEACells.core.SEACells(adata,
                  build_kernel_on=build_kernel_on,
                  n_SEACells=n_SEACells,
                  n_waypoint_eigs=n_waypoint_eigs,
                  waypt_proportion=waypoint_proportion,
                  convergence_epsilon = 1e-5)

# Initialize archetypes
model.initialize_archetypes()

model.fit(n_iter=20)


adata.obs[['SEACell']].head()

SEACell_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')
SEACell_ad

# save the results:
adata.obs.to_csv('data/Swarup_2021_singlecell.csv')
SEACell_ad.write_h5ad('data/Swarup_2021_SEACells.h5ad')


adata.write_h5ad('data/tutorial_singlecell.h5ad')
SEACell_ad.write_h5ad('data/tutorial_seacells.h5ad')


# write the components of the seacells object:

# re-load
adata = sc.read_h5ad('data/tutorial_singlecell.h5ad')
SEACell_ad = sc.read_h5ad('data/tutorial_seacells.h5ad')

################################################################################
# Save components of single-cell dataset
################################################################################

# directory to save the output files
data_dir = './'

# write obs table
adata.obs['UMAP_1'] = adata.obsm['X_umap'][:,0]
adata.obs['UMAP_2'] = adata.obsm['X_umap'][:,1]
adata.obs.to_csv('{}/tutorial_singlecell_obs.csv'.format(data_dir))

# write var table:
adata.var.to_csv('{}/tutorial_singlecell_var.csv'.format(data_dir))

# save the sparse matrix for Seurat:
X = adata.raw.X
X = scipy.sparse.csr_matrix.transpose(X)
io.mmwrite('{}/tutorial_singlecell.mtx'.format(data_dir), X)

pd.DataFrame(adata.obsm['X_pca']).to_csv('{}/tutorial_singlecell_pca.csv'.format(data_dir))


################################################################################
# Save components of metacell dataset
################################################################################

# write obs table
SEACell_ad.obs.to_csv('{}/tutorial_seacells_obs.csv'.format(data_dir))


# save the sparse matrix for Seurat:
X = SEACell_ad.X
X = scipy.sparse.csr_matrix.transpose(X)
io.mmwrite('{}/tutorial_seacells.mtx'.format(data_dir), X)
