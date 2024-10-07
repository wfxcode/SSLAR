#' Refine Cell Proportions Using k-Nearest Neighbors
#'
#' This function refines cell proportions by considering the k-nearest neighbors in both spatial and compositional spaces.
#'
#' @param st_position A matrix of spatial positions.
#' @param pred_comp A matrix of predicted cell proportions.
#' @param k The number of nearest neighbors to consider.
#' @param metric The distance metric to use (default is 'manhattan').
#' @param min_cont The minimum proportion threshold (default is 0.01).
#' @param verbose Whether to print process information (default is TRUE).
#'
#' @return A list containing:
#' - `refine_comp`: The refined cell proportions.
#' - `pred_comp`: The original predicted cell proportions.
#'
#' @export
#'
#' @examples
#' # Example usage
#' st_position <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2, byrow = TRUE)
#' pred_comp <- matrix(runif(12), ncol = 3)
#' result <- refine_label(st_position, pred_comp, k = 3, verbose = TRUE)
#' print(result$refine_comp)
refine_label <- function(
    st_position, # Spatial positions
    pred_comp,   # Predicted cell proportions
    k = 6,       # k-nearest neighbors parameter
    metric = 'manhattan', # Distance metric
    min_cont = 0.01,
    verbose = TRUE  # Whether to print process information
) {

  # Load required libraries
  suppressMessages(require(FNN))
  suppressMessages(require(stats))
  suppressMessages(require(purrr))


  ## Calculate kNN spots per spot
  distances <- as.matrix(dist(st_position))
  differences <- as.matrix(dist(pred_comp, metric))

  # Get k-nearest neighbors indices using FNN::get.knn
  knn.indexes <- FNN::get.knn(distances, k = nrow(distances) - 1)[[1]]
  knn.diff <- FNN::get.knn(differences, k = nrow(st_position) - 1)[[1]]
  intersections <- get_top_k_intersection(knn.indexes, knn.diff, k)

  # Set row names of knn.indexes to match st_position
  rownames(knn.indexes) <- rownames(st_position)

  ## Sum cell labels for each cell and its k-nearest neighbors
  sum_comp <- matrix(0, nrow = nrow(pred_comp), ncol = ncol(pred_comp))
  for (i in seq_len(nrow(pred_comp))) {
    idx <- intersections[[i]]
    if (purrr::is_empty(idx)) {
      sum_comp[i,] <- pred_comp[i,]
    } else if (length(idx) == 1) {
      sum_comp[i,] <- pred_comp[idx,]
    } else {
      sum_comp[i,] <- colSums(pred_comp[idx,])
    }
  }

  colnames(sum_comp) <- colnames(pred_comp)

  # Normalize the summed and original label matrices
  for (i in seq_len(nrow(sum_comp))) {
    weight <- sum_comp[i,]
    sum_comp[i,] <- weight / sum(weight)

    weight <- pred_comp[i,]
    pred_comp[i,] <- weight / sum(weight)
  }

  ## Calculate distances between actual and contextual matrices
  if (verbose) message("\n=== Calculating distances in transcriptome space\n")
  dist_mm <- dist_calc(x = sum_comp, y = pred_comp)

  # Calculate alpha regularization factors based on distances
  if (verbose) message("=== Calculating alpha factors based on distances\n")
  dist.mod <- as.vector(dist_mm) <= mean(as.vector(dist_mm))
  dist.rescaled <- rescale_function(
    as.vector(dist_mm)[dist.mod],
    to = c(0, 0.5)
  )

  # Calculate intrinsic and extrinsic alpha factors
  alpha.extrin <- ifelse(test = dist.mod, yes = dist.rescaled, no = 0)
  alpha.intric <- 1 - alpha.extrin

  # Initialize the refined cell proportion matrix
  refine_comp <- matrix(0, nrow = nrow(sum_comp), ncol = ncol(sum_comp))
  for (spot in seq(nrow(pred_comp))) {
    # Weighted average of cell proportions
    refine_comp[spot, ] <- sapply(
      X = colnames(pred_comp),
      FUN = \(x) {
        stats::weighted.mean(
          c(pred_comp[spot, x], sum_comp[spot, x]),
          w = c(alpha.intric[spot], alpha.extrin[spot])
        )
      }
    )
    weight <- pred_comp[spot, ]
    weight[weight < min_cont] <- 0
    pred_comp[spot,] <- weight / sum(weight)

    weight <- refine_comp[spot, ]
    weight[weight < min_cont] <- 0
    refine_comp[spot,] <- weight / sum(weight)
  }

  # Set column and row names for the refined cell proportion matrix
  colnames(refine_comp) <- colnames(pred_comp)
  rownames(refine_comp) <- rownames(pred_comp)

  # Return the result list
  return(
    list(
      refine_comp, pred_comp
    ) %>% setNames(
      c("refine_comp", "pred_comp")
    )
  )

}

#' Get Top k Intersections
#'
#' This function gets the intersection of the top k indices from two sets of k-nearest neighbors.
#'
#' @param knn.indexes Matrix of k-nearest neighbor indices.
#' @param knn.diff Matrix of k-nearest neighbor indices based on differences.
#' @param k Number of nearest neighbors to consider.
#'
#' @return A list of intersections for each row.
#'
#' @export
#'
#' @examples
#' # Example usage
#' knn.indexes <- matrix(1:12, ncol = 4)
#' knn.diff <- matrix(13:24, ncol = 4)
#' intersections <- get_top_k_intersection(knn.indexes, knn.diff, k = 3)
#' print(intersections)
get_top_k_intersection <- function(knn.indexes, knn.diff, k) {
  # Get the top k indices for each row
  top_k_indexes <- t(apply(knn.indexes, 1, function(row) head(row, k)))
  top_k_diff <- t(apply(knn.diff, 1, function(row) head(row, k * 2)))

  # Get the intersection of top k indices
  intersections <- lapply(1:nrow(top_k_indexes), function(i) {
    intersect(top_k_indexes[i, ], top_k_diff[i, ])
  })

  return(intersections)
}

#' Calculate Distances Between Matrices
#'
#' This function calculates the distances between corresponding rows of two matrices.
#'
#' @param x First matrix.
#' @param y Second matrix.
#'
#' @return A vector of distances.
#'
#' @export
#'
#' @examples
#' # Example usage
#' x <- matrix(1:6, ncol = 2)
#' y <- matrix(7:12, ncol = 2)
#' dists <- dist_calc(x, y)
#' print(dists)
dist_calc <- function(x, y) {
  dist_mm <- matrix(0, ncol = 1, nrow = nrow(x))
  for (i in seq(nrow(x))) {
    dist_mm[i, 1] <- stats::dist(
      t(cbind(x[i, ], y[i, ])), method = "manhattan"
    )
  }
  return(dist_mm)
}

#' Rescale Values
#'
#' This function rescales values to a specified range.
#'
#' @param x Vector of values to rescale.
#' @param to Range to rescale to (default is c(0, 1)).
#' @param from Range to rescale from (default is the range of x).
#'
#' @return Rescaled values.
#'
#' @export
#'
#' @examples
#' # Example usage
#' x <- c(1, 2, 3, 4, 5)
#' rescaled_x <- rescale_function(x, to = c(0, 10))
#' print(rescaled_x)
rescale_function <- function(
    x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)
) {
  (x - from[1]) / diff(from) * diff(to) + to[1]
}
