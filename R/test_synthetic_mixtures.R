#' Test Synthetic Mixtures Performance
#'
#' This function takes in a Seurat object with several tuning parameters and assesses its performance on synthetic test spots.
#'
#' @param se_sc Object of class Seurat with the scRNAseq data.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @param cluster_markers Object of class dataframe obtained from the function Seurat::FindAllMarkers().
#' @param n_syn_mixt Object of class integer specifying how many synthetic mixtures to generate.
#' @param ntop Object of class "numeric" or NULL; number of unique markers per cluster used to seed the model, by default NULL. If NULL it uses all of them.
#' @param transf Transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), sct (Seurat::SCTransform), raw (no transformation applied). By default CPM.
#' @param method Object of class character; Type of method to use to find W and H. Look at NMF package for the options and specifications, by default nsNMF.
#' @param assay Character string specifying the assay to use.
#' @param slot Character string specifying the slot to use.
#' @param numiters Integer specifying the number of iterations for NMF.
#' @param lam Numeric specifying the regularization parameter for NMF.
#' @param min_cont Object of class numeric; Indicates the minimum contribution we expect from a cell in that spot. Since we're working with proportions by setting 0.01, by default, means that we will accept those cell types whose weight coefficient is at least 1\% of the total.
#'
#' @return This function returns a list where the first element is a list with the NMF model trained and the cell labels, the second is a list with the raw_statistics.
#' @export
#' @examples
#' # Load necessary libraries
#' library(Seurat)
#' library(NMF)
#'
#'
#' # Find all markers
#' cluster_markers <- FindAllMarkers(se_sc)
#'
#' # Run the function
#' results <- test_synthetic_mixtures(
#'   se_sc = se_sc,
#'   clust_vr = "seurat_clusters",
#'   cluster_markers = cluster_markers,
#'   n_syn_mixt = 1000,
#'   ntop = 10,
#'   transf = "uv",
#'   method = "nsNMF",
#'   min_cont = 0.07,
#'   assay = "RNA",
#'   slot = "counts",
#'   numiters = 100,
#'   lam = 100
#' )
#'
#' # View the results
#' print(results$model)
#' print(results$performance)
#'

test_synthetic_mixtures <- function(se_sc,
                                    clust_vr,
                                    cluster_markers,
                                    n_syn_mixt = 1000,
                                    ntop = NULL,
                                    transf = "uv",
                                    method = "ssNMF",
                                    min_cont = 0.07,
                                    assay = "RNA",
                                    slot = "counts",
                                    numiters = 100,
                                    lam = 100) {

  # Generate test data
  print(sprintf("Generating %s synthetic test mixtures", n_syn_mixt))
  test_spot_ls <- test_spot_fun(se_obj = se_sc,
                                clust_vr = clust_vr,
                                n = n_syn_mixt)

  test_spot_counts <- as.matrix(test_spot_ls$topic_profiles)
  colnames(test_spot_counts) <- paste("mixt", 1:ncol(test_spot_counts), sep = "_")
  test_spot_metadata <- test_spot_ls$cell_composition

  ###################
  #### Train NMF ####
  ###################
  print("Train NMF")
  # Train the NMF model
  nmf_mod_ls <- topic_identification(cluster_markers = cluster_markers,
                                     se_sc = se_sc,
                                     mtrx_spatial = test_spot_counts,
                                     ntop = ntop,
                                     transf = transf,
                                     clust_vr = clust_vr,
                                     nmf = method,
                                     numiters = numiters,
                                     lam = lam,
                                     assay = assay,
                                     slot = slot)

  #################################
  #### Get mixture composition ####
  #################################
  # Run test spots through the basis to get the pertinent coefficients. To do this for every spot we are going to set up a system of linear equations where we need to find the coefficient, we will use non-negative least squares to determine the best coefficient fit.
  # Get cell type specific topic profiles
  ct_topic_profiles <- topic_profile_per_cluster_nmf(h = coef(nmf_mod_ls[[1]]),
                                                     train_cell_clust = nmf_mod_ls[[2]])

  print("Deconvolute synthetic spots")
  # Perform deconvolution of the capture location mixtures
  pred_comp <- cell_type_annotation(nmf_mod = nmf_mod_ls[[1]],
                                    mixture_transcriptome = test_spot_counts,
                                    transf = transf,
                                    reference_profiles = ct_topic_profiles,
                                    min_cont = min_cont,
                                    regr = 'lar',
                                    regr_type = 'nnls')
  pred_comp = pred_comp$pred_comp
  ################################
  #### Performance statistics ####
  ################################
  # Performance
  ct_cols <- colnames(pred_comp)[which(colnames(pred_comp) != "res_ss")]
  performance <- test_performances(test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
                                   spot_composition_mtrx = pred_comp[, ct_cols])

  return(list(model = nmf_mod_ls, performance = performance))
}
