

#' BinPseudotime
#'
#' Makes evenly-sized bins of cells along a pseudotime trajectory.
#'
#' @param seurat_obj A Seurat object
#' @param n_bins number of hub genes to use in the UMAP computation
#' @param pseudotime_col
#' @keywords scRNA-seq
#' @export
#' @examples
BinPseudotime <- function(
  seurat_obj,
  n_bins = 50,
  pseudotime_col = 'pseudotime'
){
  # cut pseudotime into bins
  for(cur_bins in n_bins){
    bin_name <- paste0(pseudotime_col, '_bins_', cur_bins)
    seurat_obj@meta.data[[bin_name]] <- as.numeric(seurat_obj@meta.data[[pseudotime_col]]) %>% Hmisc::cut2(g=cur_bins)

    # make slot in the seurat misc related to these bins:
    seurat_obj@misc[[bin_name]] <- list()

  }

  # return seurat obj
  seurat_obj
}
