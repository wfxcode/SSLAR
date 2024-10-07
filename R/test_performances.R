#' Test the performance of the model on synthetic generated test spots.
#'
#' @param test_spots_metadata_mtrx Object of class matrix containing the ground truth composition of each spot as obtained from the function syn_spot_comb_topic_fun.R, 2nd element.
#' @param spot_composition_mtrx Object of class matrix with the predicted topic probability distributions for each spot.
#' @return This function returns a list with TP, TN, FP, FN, and the Jensen-Shannon Divergence index.
#' @export
#' @examples
#' # Example usage
#' # Assuming test_spots_metadata_mtrx and spot_composition_mtrx are already defined
#' performance_metrics <- test_performances(test_spots_metadata_mtrx, spot_composition_mtrx)
#' print(performance_metrics)

test_performances <- function(test_spots_metadata_mtrx,
                              spot_composition_mtrx) {
  # Check variables
  if (!is.matrix(test_spots_metadata_mtrx)) stop("ERROR: test_spots_metadata_mtrx must be a matrix object!")
  if (!is.matrix(spot_composition_mtrx)) stop("ERROR: spot_composition_mtrx must be a matrix object!")

  colnames(spot_composition_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                          x = colnames(spot_composition_mtrx),
                                          perl = TRUE)
  colnames(test_spots_metadata_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                             x = colnames(test_spots_metadata_mtrx),
                                             perl = TRUE)
  # Load required packages
  suppressMessages(require(philentropy))

  ##### Get TRUE JSD between real-predicted proportions #####
  true_rmse_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), ncol = 1)
  true_jsd_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), ncol = 1)
  true_ps_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), ncol = 1)
  tp <- 0; tn <- 0; fp <- 0; fn <- 0
  true_tmp = sweep(test_spots_metadata_mtrx, 1, rowSums(test_spots_metadata_mtrx), '/')

  for (i in seq_len(nrow(test_spots_metadata_mtrx))) {
    # Create matrix to feed to JSD
    x <- rbind(test_spots_metadata_mtrx[i, ],
               spot_composition_mtrx[i, ])

    true_rmse_mtrx[i, 1] <- RMSE(true_tmp[i, ], spot_composition_mtrx[i, ])

    # Calculate JSD and save it in true_JSD_mtrx
    if(sum(spot_composition_mtrx[i, ]) > 0) {
      true_jsd_mtrx[i, 1] <- suppressMessages(philentropy::JSD(x = x, unit = "log2",
                                                               est.prob = "empirical"))
    } else {
      true_jsd_mtrx[i, 1] <- 1
    }
    true_ps_mtrx[i, 1] <- philentropy::lin.cor(spot_composition_mtrx[i, ], test_spots_metadata_mtrx[i, ])

    #### Calculate TP-TN-FP-FN ####
    for (index in colnames(test_spots_metadata_mtrx)) {
      if (x[1, index] > 0 & x[2, index] > 0) {
        tp <- tp + 1
      } else if (x[1, index] == 0 & x[2, index] == 0) {
        tn <- tn + 1
      } else if (x[1, index] > 0 & x[2, index] == 0) {
        fn <- fn + 1
      } else if (x[1, index] == 0 & x[2, index] > 0) {
        fp <- fp + 1
      }
    }; rm(index)
  }; rm(i)

  #### Performance metrics ####
  accuracy <- round((tp + tn) / (tp + tn + fp + fn), 4)
  sensitivity <- round(tp / (tp + fn), 4)
  specificity <- round(tn / (tn + fp), 4)
  precision <- round(tp / (tp + fp), 4)
  recall <- round(tp / (tp + fn), 4)
  F1 <- round(2 * ((precision * recall) / (precision + recall)), 4)

  quants_rmse <- round(quantile(matrixStats::rowMins(true_rmse_mtrx,
                                                     na.rm = TRUE),
                                c(0.25, 0.5, 0.75)), 4)

  quants_jsd <- round(quantile(matrixStats::rowMins(true_jsd_mtrx,
                                                    na.rm = TRUE),
                               c(0.25, 0.5, 0.75)), 4)
  quants_ps <- round(quantile(matrixStats::rowMins(true_ps_mtrx,
                                                   na.rm = TRUE),
                              c(0.25, 0.5, 0.75)), 4)

  cat(sprintf("The following summary statistics are obtained:
              Accuracy: %s,
              Sensitivity: %s,
              Specificity: %s,
              Precision: %s,
              F1 score: %s,
              RMSE: %s[%s-%s],
              JSD quantiles: %s[%s-%s],
              Pearson: %s[%s-%s]",
              accuracy, sensitivity, specificity, precision, F1,
              quants_rmse[[2]], quants_rmse[[1]], quants_rmse[[3]],
              quants_jsd[[2]], quants_jsd[[1]], quants_jsd[[3]],
              quants_ps[[2]], quants_ps[[1]], quants_ps[[3]]), sep = "\n")

  rmse_range = sprintf("[%s-%s]", quants_rmse[[1]], quants_rmse[[3]])
  jsd_range = sprintf("[%s-%s]", quants_jsd[[1]], quants_jsd[[3]])
  pearson_range = sprintf("[%s-%s]", quants_ps[[1]], quants_ps[[3]])
  return(list(Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity,
              Precision = precision, F1_score = F1, RMSE = quants_rmse,
              pearson = quants_ps, JSD = quants_jsd,
              RMSE_list = true_rmse_mtrx, JSD_list = true_jsd_mtrx, pearson_list = true_ps_mtrx
  ))
}
#' Test the performance of the model for each class on synthetic generated test spots.
#'
#' @param test_spots_metadata_mtrx Object of class matrix containing the ground truth composition of each spot as obtained from the function syn_spot_comb_topic_fun.R, 2nd element.
#' @param spot_composition_mtrx Object of class matrix with the predicted topic probability distributions for each spot.
#' @return This function returns a list with performance metrics for each class.
#' @export
#' @examples
#' class_performance_metrics <- each_class_performances(test_spots_metadata_mtrx, spot_composition_mtrx)
#' print(class_performance_metrics)

each_class_performances <- function(test_spots_metadata_mtrx,
                                    spot_composition_mtrx) {
  # Check variables
  if (!is.matrix(test_spots_metadata_mtrx)) stop("ERROR: test_spots_metadata_mtrx must be a matrix object!")
  if (!is.matrix(spot_composition_mtrx)) stop("ERROR: spot_composition_mtrx must be a matrix object!")

  colnames(spot_composition_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                          x = colnames(spot_composition_mtrx),
                                          perl = TRUE)
  colnames(test_spots_metadata_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                             x = colnames(test_spots_metadata_mtrx),
                                             perl = TRUE)
  # Load required packages
  suppressMessages(require(philentropy))

  ##### Get TRUE JSD between real-predicted proportions #####
  true_accuracy_mtrx <- matrix(nrow = ncol(test_spots_metadata_mtrx), ncol = 1)
  true_sensitivity_mtrx <- matrix(nrow = ncol(test_spots_metadata_mtrx), ncol = 1)
  true_specificity_mtrx <- matrix(nrow = ncol(test_spots_metadata_mtrx), ncol = 1)
  true_precision_mtrx <- matrix(nrow = ncol(test_spots_metadata_mtrx), ncol = 1)
  true_F1_mtrx <- matrix(nrow = ncol(test_spots_metadata_mtrx), ncol = 1)
  rownames(true_accuracy_mtrx) <- rownames(true_sensitivity_mtrx) <- rownames(true_specificity_mtrx) <-
    rownames(true_precision_mtrx) <- rownames(true_F1_mtrx) <- colnames(test_spots_metadata_mtrx)

  true_rmse_mtrx <- matrix(nrow = ncol(test_spots_metadata_mtrx), ncol = 1)
  true_jsd_mtrx <- matrix(nrow = ncol(test_spots_metadata_mtrx), ncol = 1)
  true_ps_mtrx <- matrix(nrow = ncol(test_spots_metadata_mtrx), ncol = 1)
  rownames(true_rmse_mtrx) <- rownames(true_jsd_mtrx) <-
    rownames(true_ps_mtrx) <- colnames(test_spots_metadata_mtrx)

  true_tmp <- sweep(test_spots_metadata_mtrx, 1, rowSums(test_spots_metadata_mtrx), '/')

  for (i in seq_len(ncol(test_spots_metadata_mtrx))) {
    # Create matrix to feed to JSD
    x <- rbind(test_spots_metadata_mtrx[, i],
               spot_composition_mtrx[, i])

    # Calculate JSD and save it in true_JSD_mtrx
    if(sum(spot_composition_mtrx[, i]) > 0) {
      true_jsd_mtrx[i, 1] <- round(suppressMessages(philentropy::JSD(x = x, unit = "log2",
                                                                     est.prob = "empirical")), 4)
    } else {
      true_jsd_mtrx[i, 1] <- 1
    }
    true_rmse_mtrx[i, 1] <- round(RMSE(true_tmp[, i], spot_composition_mtrx[, i]), 4)
    true_ps_mtrx[i, 1] <- round(philentropy::lin.cor(spot_composition_mtrx[, i], test_spots_metadata_mtrx[, i]), 4)

    #### Calculate TP-TN-FP-FN ####
    tp <- 0; tn <- 0; fp <- 0; fn <- 0
    for (index in rownames(test_spots_metadata_mtrx)) {
      if (x[1, index] > 0 & x[2, index] > 0) {
        tp <- tp + 1
      } else if (x[1, index] == 0 & x[2, index] == 0) {
        tn <- tn + 1
      } else if (x[1, index] > 0 & x[2, index] == 0) {
        fn <- fn + 1
      } else if (x[1, index] == 0 & x[2, index] > 0) {
        fp <- fp + 1
      }
    }; rm(index)

    accuracy <- round((tp + tn) / (tp + tn + fp + fn), 4)
    sensitivity <- round(tp / (tp + fn), 4)
    specificity <- round(tn / (tn + fp), 4)
    precision <- round(tp / (tp + fp), 4)
    F1 <- round(2 * ((precision * sensitivity) / (precision + sensitivity)), 4)
    true_accuracy_mtrx[i, 1] <- accuracy
    true_sensitivity_mtrx[i, 1] <- sensitivity
    true_specificity_mtrx[i, 1] <- specificity
    true_precision_mtrx[i, 1] <- precision
    true_F1_mtrx[i, 1] <- F1

    cell_type <- colnames(test_spots_metadata_mtrx)[i]
    cat(sprintf("The %s cell type performance are obtained:
              F1 score: %s,
              RMSE: %s,
              JSD quantiles: %s,
              Pearson: %s",
                cell_type, true_F1_mtrx[i, 1], true_rmse_mtrx[i, 1],
                true_jsd_mtrx[i, 1], true_ps_mtrx[i, 1]), sep = "\n")
  }; rm(i)

  #### Performance metrics ####
  return(list(Accuracy = true_accuracy_mtrx, Sensitivity = true_sensitivity_mtrx,
              Specificity = true_specificity_mtrx, Precision = true_precision_mtrx,
              F1 = true_F1_mtrx, PCC = true_ps_mtrx,
              RMSE = true_rmse_mtrx, JSD = true_jsd_mtrx
  ))
}

RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

