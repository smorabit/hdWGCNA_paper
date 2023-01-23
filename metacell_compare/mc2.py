# following this tutorial: https://metacells.readthedocs.io/en/latest/Metacells_Vignette.html


import numpy as np
import pandas as pd
import scanpy as sc
import SEACells

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import scipy
from scipy import io
import os

import anndata as ad
import metacells as mc
import scipy.sparse as sp
import seaborn as sb
from math import hypot



adata = sc.read_h5ad('data/cd34_multiome_rna.h5ad')

# need to run these utilities functions to fix the counts matrix
X = adata.X
mc.utilities.typing.sum_duplicates(X)
mc.utilities.typing.sort_indices(X)
adata.X = X

# set the raw counts matrix
raw_ad = sc.AnnData(adata.X)
raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
adata.raw = raw_ad

# name of the dataset
mc.ut.set_name(adata, 'cd34')


################################################################################
# Cleaning genes
################################################################################

excluded_gene_names = ['IGHMBP2', 'IGLL1', 'IGLL5', 'IGLON5', 'NEAT1', 'TMSB10', 'TMSB4X']
excluded_gene_patterns = ['MT-.*']

mc.pl.analyze_clean_genes(adata,
                          excluded_gene_names=excluded_gene_names,
                          excluded_gene_patterns=excluded_gene_patterns,
                          random_seed=123456)

# combine into a clean gene mask
mc.pl.pick_clean_genes(adata)

################################################################################
# Clean cells
################################################################################

full = adata


properly_sampled_min_cell_total = 800
properly_sampled_max_cell_total = 15000

total_umis_of_cells = mc.ut.get_o_numpy(full, name='__x__', sum=True)

too_small_cells_count = sum(total_umis_of_cells < properly_sampled_min_cell_total)
too_large_cells_count = sum(total_umis_of_cells > properly_sampled_max_cell_total)

too_small_cells_percent = 100.0 * too_small_cells_count / len(total_umis_of_cells)
too_large_cells_percent = 100.0 * too_large_cells_count / len(total_umis_of_cells)

print(f"Will exclude %s (%.2f%%) cells with less than %s UMIs"
      % (too_small_cells_count,
         too_small_cells_percent,
         properly_sampled_min_cell_total))
print(f"Will exclude %s (%.2f%%) cells with more than %s UMIs"
      % (too_large_cells_count,
         too_large_cells_percent,
         properly_sampled_max_cell_total))


properly_sampled_max_excluded_genes_fraction = 0.1

excluded_genes_data = mc.tl.filter_data(full, var_masks=['~clean_gene'])[0]
excluded_umis_of_cells = mc.ut.get_o_numpy(excluded_genes_data, name='__x__', sum=True)
excluded_fraction_of_umis_of_cells = excluded_umis_of_cells / total_umis_of_cells

too_excluded_cells_count = sum(excluded_fraction_of_umis_of_cells > properly_sampled_max_excluded_genes_fraction)

too_excluded_cells_percent = 100.0 * too_excluded_cells_count / len(total_umis_of_cells)

print(f"Will exclude %s (%.2f%%) cells with more than %.2f%% excluded gene UMIs"
      % (too_excluded_cells_count,
         too_excluded_cells_percent,
         100.0 * properly_sampled_max_excluded_genes_fraction))


mc.pl.analyze_clean_cells(
    full,
    properly_sampled_min_cell_total=properly_sampled_min_cell_total,
    properly_sampled_max_cell_total=properly_sampled_max_cell_total,
    properly_sampled_max_excluded_genes_fraction=properly_sampled_max_excluded_genes_fraction)


mc.pl.pick_clean_cells(full)

clean = mc.pl.extract_clean_data(full)


################################################################################
# Forbidden genes
################################################################################

suspect_gene_names = ['PCNA', 'MKI67', 'TOP2A', 'HIST1H1D',
                      'FOS', 'JUN', 'HSP90AB1', 'HSPA1A',
                      'ISG15', 'WARS' ]
suspect_gene_patterns = [ 'MCM[0-9]', 'SMC[0-9]', 'IFI.*' ]
suspect_genes_mask = mc.tl.find_named_genes(clean, names=suspect_gene_names,
                                            patterns=suspect_gene_patterns)
suspect_gene_names = sorted(clean.var_names[suspect_genes_mask])


mc.pl.relate_genes(clean, random_seed=123456)

# which groups of genes contain sus genes
module_of_genes = clean.var['related_genes_module']
suspect_gene_modules = np.unique(module_of_genes[suspect_genes_mask])
suspect_gene_modules = suspect_gene_modules[suspect_gene_modules >= 0]
print(suspect_gene_modules)


similarity_of_genes = mc.ut.get_vv_frame(clean, 'related_genes_similarity')
for gene_module in suspect_gene_modules:
    module_genes_mask = module_of_genes == gene_module
    similarity_of_module = similarity_of_genes.loc[module_genes_mask, module_genes_mask]
    similarity_of_module.index = \
    similarity_of_module.columns = [
        '(*) ' + name if name in suspect_gene_names else name
        for name in similarity_of_module.index
    ]
    ax = plt.axes()
    sb.heatmap(similarity_of_module, vmin=0, vmax=1, xticklabels=True, yticklabels=True, ax=ax, cmap="YlGnBu")
    ax.set_title(f'Gene Module {gene_module}')
    plt.savefig('figures/module_heatmap_{}.pdf'.format(gene_module))
    plt.clf()



# genes that are correlated with the known forbidden genes
forbidden_genes_mask = suspect_genes_mask
for gene_module in [17, 19, 24, 78, 113, 144, 149]:
    module_genes_mask = module_of_genes == gene_module
    forbidden_genes_mask |= module_genes_mask

forbidden_gene_names = sorted(clean.var_names[forbidden_genes_mask])
print(len(forbidden_gene_names))
print(' '.join(forbidden_gene_names))




################################################################################
# computing metacells
################################################################################

max_parallel_piles = mc.pl.guess_max_parallel_piles(clean)
print(max_parallel_piles)
mc.pl.set_max_parallel_piles(max_parallel_piles)



with mc.ut.progress_bar():
    mc.pl.divide_and_conquer_pipeline(clean,
                                      forbidden_gene_names=forbidden_gene_names,
                                      #target_metacell_size=...,
                                      random_seed=123456)


metacells = mc.pl.collect_metacells(clean, name='cd34.metacells')

# output dir
data_dir = './'

# write the h5ad file
metacells.write_h5ad('{}/tutorial_MC2.h5ad'.format(data_dir))


# write obs/var tables
metacells.obs.to_csv('{}/tutorial_MC2_obs.csv'.format(data_dir))
metacells.var.to_csv('{}/tutorial_MC2_var.csv'.format(data_dir))

# save the sparse matrix for Seurat:
X = metacells.X
X = scipy.sparse.csr_matrix(np.transpose(X).astype(int))
io.mmwrite('{}/tutorial_MC2.mtx'.format(data_dir), X)
