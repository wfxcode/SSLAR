#' This functions carries out the ssNMF and returns and NMF object.
#'
#' @param se_sc Object of class Seurat with the scRNAseq data.
#' @param mtrx_spatial Object of class matrix of shape GENESxSPOT.
#' @param cluster_markers Object of class dataframe obtained from the function Seurat::FindAllMarkers()
#' @param clust_vr Object of class character; Name of the variable containing the cell clustering.
#' @param ntop Object of class "numeric" or NULL; number of unique markers per cluster used to seed the model, by default NULL If NULL it uses all of them.
#' @param transf Transformation to normalize the count matrix: uv (unit variance), raw (no transformation applied). By default UV.
#' @param nmf Object of class character; Type of method to us to find W and H. Look at NMF package for the options and specifications, by default nsNMF.
#' @param hvg Object of class numeric or NULL; Number of highly variable genes to use on top of the marker genes, if NULL then it is completely unsupervised and use top 3000 HVG.
#' @param assay Object of class character; From which assay to grab the expression data to train the model, by default "RNA".
#' @param slot Object of class character; From which slot to grab the expression data to train the model, by default "counts".
#' @return This function returns a list with the initialized matrices H and W.
#' @export
#' @examples
#'

topic_identification <- function(se_sc,
                      mtrx_spatial,
                      cluster_markers,
                      clust_vr,
                      ntop = NULL,
                      transf = "uv",
                      nmf = "ssNMF",
                      numiters = 100,
                      lam = 100,
                      assay = "RNA",
                      slot = "counts") {
  # print(11)
  # Check variables
  if (!is.data.frame(cluster_markers)) stop("ERROR: cluster_markers_all must be a data frame object returned from Seurat::FindAllMarkers()!")
  if (is(se_sc) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  # if (!is(mtrx_spatial, "matrix")) stop("ERROR: mtrx_spatial must be a matrix object!")
  if (!(is.numeric(ntop) | is.null(ntop))) stop("ERROR: ntop must be numeric or NULL!")
  if (!is.character(transf)) stop("ERROR: transf must be a character string!")
  if (!(nmf == "nsNMF" | nmf == "ssNMF")) stop("ERROR: method must be 'nsNMF' or 'ssNMF'!")
  # if (!(is.numeric(hvg))) stop("ERROR: hvg must be a numeric ")
  if (!(is.numeric(numiters))) stop("ERROR: numiters must be a numeric ")

  # Loading libraries
  suppressMessages(require(NMF))
  suppressMessages(require(Seurat))
  suppressMessages(require(Matrix))
  suppressMessages(require(dplyr))
  suppressMessages(require(edgeR))

  print("Preparing Gene set")

  # Only train the model with genes shared between the scRNAseq and spatial data
  ## Remove rows with all gene counts 0
  mtrx_sc <- as.matrix(Seurat::GetAssayData(se_sc,
                                            assay = assay,
                                            slot = slot))
  # 19819个基因，1404个单细胞
  genes_0_sc <- which(! rowSums(mtrx_sc == 0) == ncol(mtrx_sc))
  se_sc <- se_sc[genes_0_sc, ]

  genes_0_sp <- which(! rowSums(as.matrix(mtrx_spatial) == 0) == ncol(mtrx_spatial))
  mtrx_spatial <- mtrx_spatial[genes_0_sp, ]
  # browser()
  ## Remove non intersecting genes from the scRNAseq data
  genes_spatial <- rownames(mtrx_spatial)
  genes_sc <- rownames(Seurat::GetAssayData(se_sc,
                                            assay = assay,
                                            slot = slot))

  if (length(intersect(genes_sc, genes_spatial)) < 10) stop("Not enough genes in common between the single-cell and mixture dataset.")
  se_sc <- se_sc[intersect(genes_sc, genes_spatial), ]

  # Update mtrx_sc with the intersecting genes only
  mtrx_sc <- as.matrix(Seurat::GetAssayData(se_sc,
                                            assay = assay,
                                            slot = slot))

  ## Remove non intersecting genes from the marker list
  cluster_markers <- cluster_markers[cluster_markers$gene %in% rownames(se_sc), ]

  # Normalize count matrix
  print("Normalizing count matrix")
  if (transf == "uv") {
    count_mtrx_t <- scale(t(mtrx_sc),
                          center = FALSE,
                          scale = apply(mtrx_sc, 1, sd, na.rm = TRUE))
    count_mtrx <- t(count_mtrx_t)

  } else if (transf == "raw") {
    count_mtrx <- mtrx_sc
  }

  # Rank of the model equals the number of cell types
  k <- length(unique(se_sc@meta.data[, clust_vr]))

  #################################
  ###### Initialize the model #####
  #################################
  # Train NMF model
  start_t <- Sys.time()

  # Define initial seeding model and set the right type
  # if (method == "nsNMF") mod <- "NMFns" else mod <- "NMFstd"
  if (nmf == "nsNMF") mod <- "NMFns"

  print("Seeding initial matrices")
  # Get init seeding matrices
  # 获得种子矩阵

  init_mtrx <- seed_init_mtrx_nmf(cluster_markers = cluster_markers,
                                  se_obj = se_sc,
                                  ntop = ntop,
                                  clust_vr = clust_vr)
  # 通过maker基因计算得到W
  W = init_mtrx[["W"]]
  # 通过cell type得到H
  H = init_mtrx[["H"]]
  # browser()
  print("Training...")
  # 用种子矩阵初始化矩阵
  if(nmf == "nsNMF"){
    # Initialize the matrix with the seeded matrices
    nmf_init <- NMF::nmfModel(W = init_mtrx[["W"]],
                              H = init_mtrx[["H"]],
                              model = mod)
    nmf_mod <- NMF::nmf(x = count_mtrx,
                        rank = k,
                        seed = nmf_init,
                        method = nmf)
    # W0 = getW(count_mtrx,k)
  }else if(nmf == "ssNMF"){
    # browser()
    # switch(modelNum,print(1),print(2),print(3),print(4),print(5),print(6),print(7))
    # 1、3、5matrix，2、4、6torch
    nmf_mod = torch_ssNMF(X = count_mtrx,
                    k = k,
                    A = W,
                    S = H,
                    Y = H,
                    N = as.integer(numiters),
                    lam = as.numeric(lam),
                    mult = TRUE)

    # 行名为基因名
    s_W = nmf_mod$A
    rownames(s_W) <- rownames(W)
    # 列名为细胞名
    s_H = nmf_mod$S
    colnames(s_H) <- colnames(H)
    nmf_mod = nmfModel(W=s_W,H=s_H)
  }
  total_t <- round(difftime(Sys.time(), start_t, units = "mins"), 2)
  print(sprintf("Time to train NMF model was %smins", total_t))

  return(list(nmf_mod, as.vector(se_sc@meta.data[, clust_vr])))
}

