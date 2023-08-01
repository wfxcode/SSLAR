#' Run test spots through the basis to get the pertinent coefficients. To do this for every spot we are going to set up a system of linear equations where we need to find the coefficient, we will use non-negative least squares to determine the best coefficient fit. Returns a matrix of KxnSPOTS dimensions.
#' 已知空间转录组数据V’与每个topic的基因分布W，通过NNLS计算出每个spot的topic H’
#' @param nmf_mod Object of class dataframe obtained from the function Seurat::FindAllMarkers().
#' @param mixture_transcriptome Object of class matric of dimensions GENESxSPOTS
#' @param transf Transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), raw (no transformation applied).
#' @return This function returns a matrix with the coefficients of the spatial mixtures.
#' @export
#' @examples
#'

topic_distribution_analysis <- function(nmf_mod,
                                         mixture_transcriptome,
                                         transf,
                                         regr) {

  # Check variables
  # if (!is(h, "matrix")) stop("ERROR: h must be a matric object!")
  # if (! is(train_cell_clust, "vector")) stop("ERROR: train_cell_clust must be a vector/list object!")
  # if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")

  # Loading libraries
  suppressMessages(require(nnls))
  suppressMessages(require(edgeR))

  ## Extract genes used in w, if there are genes not present add them with all 0
  keep_genes <- rownames(basis(nmf_mod))[rownames(basis(nmf_mod)) %in% rownames(mixture_transcriptome)]
  # fill_genes <- rownames(basis(nmf_mod))[! rownames(basis(nmf_mod)) %in% rownames(mixture_transcriptome)]

  mixture_transcriptome_subs <- as.matrix(mixture_transcriptome[keep_genes, ])

  # # Add 0s to those genes not detected
  # mtrx_0 <- matrix(data = 0,
  #                  nrow = length(fill_genes),
  #                  ncol = ncol(mixture_transcriptome_subs))
  # rownames(mtrx_0) <- fill_genes
  #
  # mixture_transcriptome_subs <- rbind(mixture_transcriptome_subs, mtrx_0)

  if (transf == "cpm") {
    count_mtrx <- edgeR::cpm(mixture_transcriptome_subs,
                             normalized.lib.sizes = FALSE)

  } else if (transf == "uv") {
    count_mtrx <- scale(t(mixture_transcriptome_subs),
                        center = FALSE,
                        scale = apply(mixture_transcriptome_subs, 1, sd, na.rm = TRUE))
    count_mtrx <- t(count_mtrx)
    # Set all 0s if the row is NA
    pos_0 <- which(rowSums(is.na(count_mtrx)) == ncol(count_mtrx))
    count_mtrx[pos_0, ] <- 0

  } else if (transf == "raw") {
    count_mtrx <- mixture_transcriptome_subs

  } else stop("Error non specified parameter passed for transf!")

  ##### Extract Basis matrix W #####
  W <- basis(nmf_mod)

  coef_pred <- matrix(data = NA,
                      nrow = ncol(W),
                      ncol = ncol(count_mtrx))
  colnames(coef_pred) <- colnames(count_mtrx)
  ##### Perform NNLS to get coefficients #####
  print("get coefficients!")
  total <- ncol(count_mtrx)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  # pb <- txtProgressBar(min = 0, max = total, style = 3)
  if(regr == "tmp"){
    # library(tlars)
    for (i in seq_len(ncol(count_mtrx))) {


      # X <- W
      # y <- count_mtrx[, i]
      # p <- ncol(X)
      # n <- nrow(X)
      # num_dummies <- p
      # dummies <- matrix(stats::rnorm(n * num_dummies), nrow = n, ncol = num_dummies)
      # XD <- cbind(X, dummies)
      # mod_tlars <- tlars_model(X = XD, y = y, num_dummies = num_dummies)
      # tlars(model = mod_tlars, early_stop = FALSE,info = FALSE)
      # tmp <- mod_tlars$get_beta()[1:p]



      coef_pred[, i] = pmax(tmp, 0)
      # coef_pred[, i] = arborescent(M = count_mtrx[, i],realW = W)
      # print(dim(coef_pred)[1])
      # print(dim(coef_pred)[2])
      # tmp = sbgd(y=count_mtrx[, i],X=W)
      # coef_pred[, i] = pmax(tmp, 0)
      # print(dim(tmp))
      # print(dim(W))
      # nnls_pred <- nnls::nnls(A = W, b = count_mtrx[, i])
      # coef_pred[, i] <- nnls_pred$x
      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }else if(regr == "nnls"){
    print("nnls")
    for (i in seq_len(ncol(count_mtrx))) {
      nnls_pred <- nnls::nnls(A = W, b = count_mtrx[, i])
      # browser()
      coef_pred[, i] <- nnls_pred$x

      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }else if(regr == "lar"){

    for (i in seq_len(ncol(count_mtrx))) {
      # browser()
      object <- lars(W,count_mtrx[, i],type=regr)
      # browser()
      tmp = object$beta[which.min(object$Cp),]
      coef_pred[, i] = pmax(tmp, 0)


      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }

  close(pb)
  return(coef_pred)
}

