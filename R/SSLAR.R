#' SSLAR (Spatial Single-cell Labeling and Annotation using Reference) Function
#'
#' This function performs spatial single-cell labeling and annotation using a reference dataset.
#' It integrates single-cell RNA sequencing (scRNA-seq) data with spatial transcriptomics (ST) data
#' to infer cell type proportions in each spatial location (spot).
#'
#' @param sc_count A matrix or data frame containing single-cell RNA sequencing counts.
#'                 Rows represent cells, and columns represent genes.
#' @param sc_label A data frame containing cell labels. It should have one column named "subclass".
#' @param st_count A matrix or data frame containing spatial transcriptomics counts.
#'                 Rows represent spots, and columns represent genes.
#' @param st_position Optional. A matrix or data frame containing spatial coordinates of the spots.
#'                    If provided, it should have two columns: "x" and "y".
#' @param maker_path Optional. Path to a precomputed marker gene list. If not provided, marker genes will be identified from the scRNA-seq data.
#' @param min_cont Minimum proportion threshold for cell types. Cell types with proportions below this threshold will be set to zero.
#' @param numiters Number of iterations for the NMF algorithm.
#' @param nmf Method for non-negative matrix factorization (NMF). Default is "ssNMF".
#' @param regr Regression method for deconvolution. Default is "lar" (Least Angle Regression).
#' @param regr_type Type of regression for deconvolution. Default is "nnls" (Non-negative Least Squares).
#'
#' @return A matrix or data frame containing the inferred cell type proportions for each spot.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' # sc_count, sc_label, st_count, and st_position are example datasets
#' result <- run.SSLAR(sc_count = sc_count,
#'                     sc_label = sc_label,
#'                     st_count = st_count,
#'                     st_position = st_position,
#'                     min_cont = 0.01)
#'
run.SSLAR <- function(sc_count,
                      sc_label,
                      st_count,
                      st_position = NULL,
                      maker_path = NULL,
                      min_cont = 0.01,
                      numiters = 100,
                      nmf = "ssNMF",
                      regr = "lar",
                      regr_type = "nnls") {
  # usethis::use_package("Seurat", min_version = "5.1.0")
  suppressMessages(require(Seurat))
  # Ensure the column name of sc_label is "subclass"
  # sc_label <- data.frame(subclass = sc_label)

  # Extract gene names from sc_count and st_count
  sc_gene = colnames(sc_count)
  st_gene = colnames(st_count)

  # Find common genes between sc_count and st_count
  intersect_gene = intersect(sc_gene, st_gene)

  # Subset the count matrices to include only the common genes
  sc_count = sc_count[, intersect_gene]
  st_count = st_count[, intersect_gene]

  # Create a Seurat object for the single-cell data
  sc_data = Seurat::CreateSeuratObject(counts = as(t(sc_count), "dgCMatrix"),
                                       meta.data = sc_label)

  # Perform SCTransform normalization on the single-cell data
  sc_data <- Seurat::SCTransform(sc_data, verbose = FALSE)
  # browser()
  # Set the active identifications to the subclass labels
  Seurat::Idents(object = sc_data) <- sc_data@meta.data$subclass

  # Load or identify marker genes
  if (!is.null(maker_path) && file.exists(maker_path)) {
    marker_genes = readRDS(maker_path)
  } else {
    marker_genes = Seurat::FindAllMarkers(object = sc_data,
                                          assay = "SCT",
                                          only.pos = TRUE)
    if (!is.null(maker_path)) {
      saveRDS(marker_genes, maker_path)
    }
  }

  # Create a Seurat object for the spatial transcriptomics data
  st_data = Seurat::CreateSeuratObject(counts = as(t(st_count), "dgCMatrix"))

  # Perform SCTransform normalization on the spatial transcriptomics data
  st_data = Seurat::SCTransform(st_data)

  # Extract the normalized count matrix from the spatial transcriptomics data
  st_count = as.matrix(st_data@assays[["SCT"]]@counts)

  # Perform topic identification using NMF
  nmf_mod_ls <- topic_identification(cluster_markers = marker_genes,
                                     se_sc = sc_data,
                                     mtrx_spatial = st_count,
                                     ntop = NULL,
                                     transf = "uv",
                                     clust_vr = 'subclass',
                                     nmf = nmf,
                                     numiters = numiters,
                                     assay = "SCT",
                                     slot = "counts")

  # Calculate the topic profiles for each cell type
  ct_topic_profiles <- topic_profile_per_cluster_nmf(h = coef(nmf_mod_ls[[1]]),
                                                     train_cell_clust = nmf_mod_ls[[2]])

  # Deconvolve the spots to infer cell type proportions
  pred_comp <- cell_type_annotation(nmf_mod = nmf_mod_ls[[1]],
                                    mixture_transcriptome = st_count,
                                    transf = "uv",
                                    reference_profiles = ct_topic_profiles,
                                    st_position = st_position,
                                    min_cont = min_cont,
                                    regr = regr,
                                    regr_type = regr_type)

  # Return the inferred cell type proportions for each spot
  return(pred_comp)
}
