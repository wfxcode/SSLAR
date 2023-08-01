#' If you wanna test the performance of your model on synthetic generated test spots you can use this function to benchmark and get a sense of the model's performance.
#  如果您想在合成生成的测试点上测试模型的性能，您可以使用此函数进行基准测试并了解模型的性能。
#' @param test_spots_metadata Object of class matrix containing the ground truth composition of each spot as obtained from the function syn_spot_comb_topic_fun.R, 2nd element.
#  类矩阵的对象，从函数 syn_spot_comb_topic_fun.R的第二个元素中获得的每个spot的真实组成。
#' @param spot_composition_mtrx Object of class matrix with the predicted topic probability distributions for each spot.
#  具有每个spot的预测主题概率分布的类矩阵对象。
#' @return This function returns a list with TP, TN, FP, FN and the Jensen-Shannon Divergence index.
#  此函数返回一个包含 TP、TN、FP、FN 和 Jensen-Shannon 散度指数的列表。
#' @export
#' @examples
#'

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
  #load required packages
  suppressMessages(require(philentropy))

  ##### Get TRUE JSD between real-predicted proportions #####
  true_rmse_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), ncol = 1)
  true_jsd_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), ncol = 1)
  true_ps_mtrx <- matrix(nrow = nrow(test_spots_metadata_mtrx), ncol = 1)
  tp <- 0; tn <- 0; fp <- 0; fn <- 0
  true_tmp = sweep(test_spots_metadata_mtrx, 1, rowSums(test_spots_metadata_mtrx), '/')
  for (i in seq_len(nrow(test_spots_metadata_mtrx))) {
    # print(i)
    # Create matrix to feed to JSD
    x <- rbind(test_spots_metadata_mtrx[i, ],
               spot_composition_mtrx[i, ])
    true_rmse_mtrx[i, 1] <- RMSE(true_tmp[i, ], spot_composition_mtrx[i, ])

    # Calculate JSD and save it in true_JSD_mtrx
    if(sum(spot_composition_mtrx[i, ]) > 0) {
      true_jsd_mtrx[i, 1] <- suppressMessages(JSD(x = x, unit = "log2",
                                                  est.prob = "empirical"))
    } else {
      true_jsd_mtrx[i, 1] <- 1
    }
    true_ps_mtrx[i, 1] <- lin.cor(spot_composition_mtrx[i, ], test_spots_metadata_mtrx[i, ])
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
  # 保存JSD

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

  # cat("raw statistics are returned in the list - TP, TN, FP, FN, JSD quantiles",
  #     sep = "\n")
  rmse_range = sprintf("[%s-%s]", quants_rmse[[1]], quants_rmse[[3]])
  jsd_range = sprintf("[%s-%s]", quants_jsd[[1]], quants_jsd[[3]])
  pearson_range = sprintf("[%s-%s]", quants_ps[[1]], quants_ps[[3]])
  return(list(Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity,
              Precision = precision, F1_score = F1,RMSE = quants_rmse[[2]],
              pearson = quants_ps[[2]], JSD = quants_jsd[[2]],
              # JSD_range = jsd_range, pearson_range = pearson_range,
              RMSE_list = true_rmse_mtrx,JSD_list = true_jsd_mtrx, pearson_list = true_ps_mtrx
  ))
}
RMSE = function(m, o){
  sqrt(mean((m - o)^2))
}

