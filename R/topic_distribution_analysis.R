#' Run test spots through the basis to get the pertinent coefficients. To do this for every spot we are going to set up a system of linear equations where we need to find the coefficient, we will use non-negative least squares to determine the best coefficient fit. Returns a matrix of KxnSPOTS dimensions.
#'
#' @param nmf_mod Object of class NMFModel obtained from the NMF package.
#' @param mixture_transcriptome Matrix of dimensions GENESxSPOTS.
#' @param regr Regression method to use: "nnls" (non-negative least squares), "lar" (least angle regression).
#' @param transf Transformation to normalize the count matrix: "cpm" (Counts per million), "uv" (unit variance), "raw" (no transformation applied).
#'
#' @return This function returns a matrix with the coefficients of the spatial mixtures.
#' @export
#' @examples
#' # Example usage
#' # Assuming nmf_mod and mixture_transcriptome are already defined
#' coef_matrix <- topic_distribution_analysis(nmf_mod, mixture_transcriptome, transf = "uv", regr = "nnls")
#' print(coef_matrix)

topic_distribution_analysis <- function(nmf_mod,
                                        mixture_transcriptome,
                                        transf = "uv",
                                        regr = "lar") {

  # Check variables
  if (!inherits(nmf_mod, "NMF")) stop("ERROR: nmf_mod must be an NMFModel object!")
  if (!is.matrix(mixture_transcriptome)) stop("ERROR: mixture_transcriptome must be a matrix object!")

  # Loading libraries
  suppressMessages(require(nnls))
  suppressMessages(require(lars))

  ## Extract genes used in W, if there are genes not present add them with all 0
  keep_genes <- rownames(basis(nmf_mod))[rownames(basis(nmf_mod)) %in% rownames(mixture_transcriptome)]
  mixture_transcriptome_subs <- as.matrix(mixture_transcriptome[keep_genes, ])

  # Transform the count matrix
  if (transf == "uv") {
    count_mtrx <- scale(t(mixture_transcriptome_subs), center = FALSE, scale = apply(mixture_transcriptome_subs, 1, sd, na.rm = TRUE))
    count_mtrx <- t(count_mtrx)
    # Set all 0s if the row is NA
    pos_0 <- which(rowSums(is.na(count_mtrx)) == ncol(count_mtrx))
    count_mtrx[pos_0, ] <- 0
  } else if (transf == "raw") {
    count_mtrx <- mixture_transcriptome_subs
  } else {
    stop("Error: Invalid parameter passed for transf!")
  }

  ##### Extract Basis matrix W #####
  W <- basis(nmf_mod)

  # Initialize coefficient matrix

  coef_pred <- matrix(data = NA, nrow = ncol(W), ncol = ncol(count_mtrx))
  colnames(coef_pred) <- colnames(count_mtrx)

  ##### Perform NNLS to get coefficients #####
  print("Getting coefficients!")
  total <- ncol(count_mtrx)
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  if (regr == "nnls") {
    for (i in seq_len(ncol(count_mtrx))) {
      nnls_pred <- nnls::nnls(A = W, b = count_mtrx[, i])
      coef_pred[, i] <- nnls_pred$x
      setTxtProgressBar(pb, i)
    }
  } else if (regr == "lar") {
    for (i in seq_len(ncol(count_mtrx))) {
      object <- lars(W, count_mtrx[, i], type = "lar")
      tmp <- object$beta[which.min(object$Cp), ]
      coef_pred[, i] <- pmax(tmp, 0)
      if (sum(coef_pred[, i]) == 0) {
        nnls_pred <- nnls::nnls(A = W, b = count_mtrx[, i])
        coef_pred[, i] <- nnls_pred$x
      }
      setTxtProgressBar(pb, i)
    }
  } else {
    stop("Error: Invalid parameter passed for regr!")
  }

  close(pb)
  return(coef_pred)
}
