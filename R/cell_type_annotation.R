#' @title
#' Run mixtures through the NMF model to get the cell type composition.
#'
#' @param nmf_mod Object of class NMF containing the trained nNMF model.
#' @param mixture_transcriptome Matrix of dimensions GENES x SPOTS containing the spatial transcriptomics data.
#' @param transf Character string specifying the transformation to normalize the count matrix: "uv" (unit variance) or "raw" (no transformation).
#' @param reference_profiles Matrix of dimensions TOPICS x CELLS containing the cell type profiles.
#' @param st_position Optional. Matrix containing the spatial coordinates of the spots. If provided, it should have two columns: "x" and "y".
#' @param min_cont Numeric value indicating the minimum contribution expected from a cell in that spot. Default is 0.01, meaning that cell types with a weight coefficient of at least 1% of the total will be accepted.
#' @param regr Character string specifying the regression method for deconvolution: "lar" (Least Angle Regression) or "nnls" (Non-negative Least Squares). Default is "lar".
#' @param regr_type Character string specifying the type of regression: "nnls", "lar", or "none". Default is "nnls".
#'
#' @return A matrix with the coefficients of the spatial mixtures, representing the inferred cell type proportions for each spot.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' # nmf_mod, mixture_transcriptome, reference_profiles, and st_position are example datasets
#' result <- cell_type_annotation(nmf_mod = nmf_mod,
#'                               mixture_transcriptome = mixture_transcriptome,
#'                               transf = "uv",
#'                               reference_profiles = reference_profiles,
#'                               st_position = st_position,
#'                               min_cont = 0.01,
#'                               regr = "lar",
#'                               regr_type = "nnls")
#'
cell_type_annotation <- function(nmf_mod,
                                 mixture_transcriptome,
                                 transf,
                                 reference_profiles,
                                 st_position = NULL,
                                 min_cont = 0.01,
                                 regr = "lar",
                                 regr_type = "nnls") {

  # Check input variables
  if (!is(nmf_mod, "NMF")) stop("ERROR: nmf_mod must be an NMF object!")
  if (!is.character(transf)) stop("ERROR: transf must be a character string!")
  if (!is.matrix(reference_profiles)) stop("ERROR: reference_profiles must be a matrix!")
  if (!is.numeric(min_cont)) stop("ERROR: min_cont must be numeric!")

  # Load required libraries
  suppressMessages(require(nnls))
  suppressMessages(require(lars))

  # Compute the topic distribution for the spatial transcriptomics data
  profile_mtrx <- topic_distribution_analysis(nmf_mod = nmf_mod,
                                              mixture_transcriptome = mixture_transcriptome,
                                              transf = transf,
                                              regr = regr)

  # Initialize the deconvolution matrix with NA values
  decon_mtrx <- matrix(data = NA,
                       nrow = ncol(profile_mtrx),
                       ncol = ncol(reference_profiles))
  colnames(decon_mtrx) <- colnames(reference_profiles)

  # Create a progress bar
  print("Deconvoluting spots")
  total <- ncol(profile_mtrx)
  pb <- txtProgressBar(min = 0, max = total, style = 3)

  # Perform deconvolution based on the specified regression type
  if (regr_type == "nnls") {
    for (i in seq_len(ncol(profile_mtrx))) {
      # Perform non-negative least squares regression
      nnls_pred <- nnls::nnls(A = reference_profiles, b = profile_mtrx[, i])
      weights <- nnls_pred$x
      weights[is.na(weights)] <- 0
      decon_mtrx[i, ] <- weights

      # Update the progress bar
      setTxtProgressBar(pb, i)
    }
  } else if (regr_type == "lar") {
    # Scale the reference profiles
    reference_profiles_t <- scale(t(reference_profiles),
                                  center = FALSE,
                                  scale = apply(reference_profiles, 1, sd, na.rm = TRUE))
    reference_profiles <- t(reference_profiles_t)

    for (i in seq_len(ncol(profile_mtrx))) {
      # Normalize the profile vector
      tmp <- profile_mtrx[, i]
      profile_mtrx[, i] <- tmp / sum(tmp)

      # Perform least angle regression
      object <- lars::lars(reference_profiles, profile_mtrx[, i], type = "lar", eps = 0)

      if (length(which.min(object$cp)) == 0) {
        tmp <- profile_mtrx[, i]
      } else {
        tmp <- object$beta[which.min(object$cp), ]
      }

      # Apply non-negative constraint
      weights <- pmax(tmp, 0)
      weights[is.na(weights)] <- 0
      decon_mtrx[i, ] <- weights

      # Update the progress bar
      setTxtProgressBar(pb, i)
    }
  } else if (regr_type == "none") {
    for (i in seq_len(ncol(profile_mtrx))) {
      # Apply non-negative constraint directly
      tmp <- profile_mtrx[, i]
      weights <- pmax(tmp, 0)
      weights[is.na(weights)] <- 0
      decon_mtrx[i, ] <- weights

      # Update the progress bar
      setTxtProgressBar(pb, i)
    }
  }
  # Close the progress bar
  close(pb)
  pred_comp <- list()
  # Refine the labels based on spatial position and minimum contribution
  if(!is.null(st_position)){
    pred_comp <- refine_label(st_position = st_position, pred_comp = decon_mtrx, k = 20, min_cont = min_cont)
  }else{
    pred_comp$pred_comp <- decon_mtrx
    for (spot in seq(nrow(pred_comp))) {
      weight <- pred_comp$pred_comp[spot, ]
      weight[weight < min_cont] <- 0
      pred_comp$pred_comp[spot,] <- weight / sum(weight)
    }
  }


  # Return the refined cell type proportions
  return(pred_comp)
}
