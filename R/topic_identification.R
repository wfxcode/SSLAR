#' This function carries out the ssNMF and returns an NMF object.
#'
#' @param se_sc Object of class Seurat with the scRNAseq data.
#' @param mtrx_spatial Object of class matrix of shape GENESxSPOT.
#' @param cluster_markers Object of class dataframe obtained from the function Seurat::FindAllMarkers().
#' @param clust_vr Object of class character; Name of the variable containing the cell clustering.
#' @param ntop Object of class "numeric" or NULL; number of unique markers per cluster used to seed the model, by default NULL. If NULL it uses all of them.
#' @param transf Transformation to normalize the count matrix: "uv" (unit variance), "raw" (no transformation applied). By default "uv".
#' @param nmf Object of class character; Type of method to use to find W and H. Look at NMF package for the options and specifications, by default "ssNMF".
#' @param numiters Number of iterations for the NMF algorithm, by default 100.
#' @param lam Regularization parameter for the NMF algorithm, by default 1.
#' @param assay Object of class character; From which assay to grab the expression data to train the model, by default "RNA".
#' @param slot Object of class character; From which slot to grab the expression data to train the model, by default "counts".
#'
#' @return This function returns a list with the initialized matrices H and W.
#' @export
#' @examples
#' # Example usage
#' # Assuming se_sc, mtrx_spatial, and cluster_markers are already defined
#' nmf_result <- topic_identification(se_sc, mtrx_spatial, cluster_markers, clust_vr = "seurat_clusters")
#' print(nmf_result)

topic_identification <- function(se_sc,
                                 mtrx_spatial,
                                 cluster_markers,
                                 clust_vr,
                                 ntop = NULL,
                                 transf = "uv",
                                 nmf = "ssNMF",
                                 numiters = 100,
                                 lam = 1,
                                 assay = "RNA",
                                 slot = "counts") {

  # Check variables
  if (!is.data.frame(cluster_markers)) stop("ERROR: cluster_markers must be a data frame object returned from Seurat::FindAllMarkers()!")
  if (!inherits(se_sc, "Seurat")) stop("ERROR: se_sc must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!(is.numeric(ntop) | is.null(ntop))) stop("ERROR: ntop must be numeric or NULL!")
  if (!is.character(transf)) stop("ERROR: transf must be a character string!")
  if (!is.character(nmf)) stop("ERROR: nmf must be a character string!")
  if (!is.numeric(numiters)) stop("ERROR: numiters must be a numeric value!")
  if (!is.numeric(lam)) stop("ERROR: lam must be a numeric value!")
  if (!is.character(assay)) stop("ERROR: assay must be a character string!")
  if (!is.character(slot)) stop("ERROR: slot must be a character string!")

  # Loading libraries
  suppressMessages(require(NMF))
  suppressMessages(require(Seurat))
  suppressMessages(require(Matrix))

  print("Preparing Gene set")

  # Only train the model with genes shared between the scRNAseq and spatial data
  mtrx_sc <- as.matrix(Seurat::GetAssayData(se_sc, assay = assay, slot = slot))
  genes_0_sc <- which(! rowSums(mtrx_sc == 0) == ncol(mtrx_sc))
  se_sc <- se_sc[genes_0_sc, ]
  # browser()
  genes_0_sp <- which(!rowSums(as.matrix(mtrx_spatial) == 0) == ncol(mtrx_spatial))
  mtrx_spatial <- mtrx_spatial[genes_0_sp, ]

  ## Remove non intersecting genes from the scRNAseq data
  genes_spatial <- rownames(mtrx_spatial)
  genes_sc <- rownames(Seurat::GetAssayData(se_sc, assay = assay, slot = slot))

  if (length(intersect(genes_sc, genes_spatial)) < 10) stop("Not enough genes in common between the single-cell and mixture dataset.")
  se_sc <- se_sc[intersect(genes_sc, genes_spatial), ]

  # Update mtrx_sc with the intersecting genes only
  mtrx_sc <- as.matrix(Seurat::GetAssayData(se_sc, assay = assay, slot = slot))

  ## Remove non intersecting genes from the marker list
  cluster_markers <- cluster_markers[cluster_markers$gene %in% rownames(se_sc), ]

  # Normalize count matrix
  print("Normalizing count matrix")
  count_mtrx <- mtrx_sc
  if (transf == "uv") {
    count_mtrx_t <- scale(t(mtrx_sc), center = FALSE, scale = apply(mtrx_sc, 1, sd, na.rm = TRUE))
    count_mtrx <- t(count_mtrx_t)
  } else if (transf == "raw") {
    count_mtrx <- mtrx_sc
  } else {
    stop("Error: Invalid parameter passed for transf!")
  }

  # Rank of the model equals the number of cell types
  k <- length(unique(se_sc@meta.data[, clust_vr]))

  #################################
  ###### Initialize the model #####
  #################################
  # Train NMF model
  start_t <- Sys.time()

  # Define initial seeding model and set the right type
  if (nmf == "nsNMF") {
    mod <- "NMFns"
  } else if (nmf == "ssNMF") {
    mod <- "NMFstd"
  } else {
    stop("Error: Invalid parameter passed for nmf!")
  }

  print("Seeding initial matrices")
  # Get init seeding matrices
  init_mtrx <- seed_init_mtrx_nmf(cluster_markers = cluster_markers,
                                  se_obj = se_sc,
                                  ntop = ntop,
                                  clust_vr = clust_vr)

  W <- init_mtrx[["W"]]
  H <- init_mtrx[["H"]]

  print("Training...")
  # Use seeds to initialize the NMF model
  if (nmf == "nsNMF") {
    nmf_init <- NMF::nmfModel(W = W, H = H, model = mod)
    nmf_mod <- NMF::nmf(x = count_mtrx, rank = k, seed = nmf_init, method = nmf)
  } else if (nmf == "ssNMF") {
    nmf_mod <- torch_ssNMF(X = count_mtrx, k = k, A = W, S = H, Y = H, N = as.integer(numiters), lam = as.numeric(lam), mult = TRUE)
    s_W <- nmf_mod$A
    rownames(s_W) <- rownames(W)
    # write.table(s_W, file = "SSLAR_left.csv", sep = ",")
    s_H <- nmf_mod$S
    colnames(s_H) <- colnames(H)
    # write.table(s_H, file = "SSLAR_right.csv", sep = ",")
    nmf_mod <- NMF::nmfModel(W = s_W, H = s_H)
  }

  total_t <- round(difftime(Sys.time(), start_t, units = "mins"), 2)
  print(sprintf("Time to train NMF model was %smins", total_t))

  return(list(nmf_mod = nmf_mod, clusters = as.vector(se_sc@meta.data[, clust_vr])))
}

#' This functions seeds initialization matrices H and W to perform NMF
#'
#' @param cluster_markers Object of class dataframe obtained from the function Seurat::FindAllMarkers().
#' @param se_obj Object of class Seurat with the data of interest.
#' @param ntop Object of class "numeric" or NULL; number of unique markers per cluster used to seed the model, by default NULL. If NULL it uses all of them.
#' @return This function returns a list with the initialized matrices H and W.
seed_init_mtrx_nmf <- function(cluster_markers,
                               se_obj,
                               clust_vr,
                               ntop = NULL) {

  #### Get dataset ready ####
  se_nmf_ready <- prep_seobj_topic_fun(se_obj = se_obj)

  # Rank of the model equals the number of cell types
  k <- length(unique(se_obj@meta.data[, clust_vr]))


  # Select all marker genes for each cluster AND compute their Z score
  if (is.null(ntop)) ntop <- max(table(cluster_markers$cluster))
  cluster_markers_cut <- suppressMessages(cut_markers(markers = cluster_markers,
                                                       ntop = ntop))

  # Select unique markers from each cluster, if there are common markers between clusters the model gets confused and could classify different clusters as belonging to the same topic just because the seeding induced it!
  cluster_markers_uniq <- lapply(unique(cluster_markers_cut$cluster), function(clust) {
    ls1 <- cluster_markers_cut[cluster_markers_cut$cluster == clust, "gene"]
    ls2 <- cluster_markers_cut[cluster_markers_cut$cluster != clust, "gene"]
    ls1_unique <- ls1[! ls1 %in% ls2]

    return(cluster_markers_cut[cluster_markers_cut$cluster == clust & cluster_markers_cut$gene %in% ls1_unique, ])
  }) %>%
    bind_rows()

  # Set seedwords from top markers. Here we are setting the weights for each topic, the words that are weighted positively are those belonging to the list of top markers for a cluster.
  # In the seedgenes matrix each row represents a topic and each column represents a gene.

  # To the NMF model we need to pass a matrix with k rows and ngene columns, where each cell has the weight of that gene for that topic. The weight we're assigning is the logFC

  # initialize matrix
  seedgenes <- matrix(nrow = k, ncol = ncol(se_nmf_ready), data = 1e-10)
  colnames(seedgenes) <- colnames(se_nmf_ready)

  # Add seeds to model, if a cluster-topic has 0 unique markers its row will be set to all 0
  for (i in seq_len(k)) {
    # print(i)
    clust_row <- cluster_markers_uniq$cluster == as.character(unique(se_obj@meta.data[, clust_vr])[[i]])
    seedgenes[i, as.character(cluster_markers_uniq[clust_row, "gene"])] = cluster_markers_uniq[clust_row, "weight"]
  }
  W <- t(seedgenes)
  ###################
  #### Seeding H ####
  ###################
  H <- matrix(data = 1e-10,
              nrow = k,
              ncol = nrow(se_nmf_ready))

  for (i in seq_len(nrow(se_nmf_ready))) {
    h_row <- which(unique(se_obj@meta.data[, clust_vr]) == se_obj@meta.data[i, clust_vr])
    H[h_row, i] <- 1
  }

  rownames(W) <- rownames(se_obj@assays$SCT@counts)

  # get matrix H
  colnames(H) <- colnames(se_obj@assays$SCT@counts)

  return(list("W" = W, "H" = H))
}

#' This functions takes in a seurat object and returns the transposed count matrix with the highly variable genes selected
#'
#' @param se_obj Object of class Seurat with the data of interest
#' @return This function returns a sparse matrix object.
prep_seobj_topic_fun <- function(se_obj) {

  # Check variables
  if (is(se_obj) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")

  #load required packages
  suppressMessages(require(Seurat))
  suppressMessages(require(Matrix))

  # 1st get from the counts matrix from the SCT pocket, raw counts
  count_mtrx <- t(as.matrix(se_obj@assays$SCT@counts))

  # 2nd reconvert the matrix to sparse format again
  count_mtrx <- Matrix::Matrix(count_mtrx, sparse = T)

  return(count_mtrx)
}
