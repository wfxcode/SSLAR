#' Run mixtures through the NMF model to get the cell type composition.
#' 通过 NMF 模型得到spot的cell type比例
#' @param nmf_mod Object of class NMF containing the trained nNMF model.
#' @param mixture_transcriptome Object of class matric of dimensions GENESxSPOTS
#' 空间转录组数据
#' @param transf Transformation to normalize the count matrix: uv (unit variance), raw (no transformation applied). By default UV.
#' 对计数矩阵进行归一化的转换：uv（单位方差），raw（未应用转换）。默认情况下为 UV。
#' @param reference_profiles Object of class matrix containing the TOPICSxCELLS Coefficient matrix from where want to get the weights. It can be cell type profiles or cell specific profiles.
#' 矩阵Q，包含每个cell type的topic
#' @param min_cont Object of class numeric; Indicates the minimum contribution we expect from a cell in that spot. Since we're working with proportions by setting 0.01, by default, means that we will accept those cell types whose weight coefficient is at least 1\% of the total.
#' @return This function returns a matrix with the coefficients of the spatial mixtures.
#' @export
#' @examples
#'

cell_type_annotation <- function(nmf_mod,
                                      mixture_transcriptome,
                                      transf,
                                      reference_profiles,
                                      min_cont = 0.01,
                                      regr = "lar",
                                      regr_type = "nnls") {

  # Check variables
  if (!is(nmf_mod, "NMF")) stop("ERROR: nmf_mod must be an NMF object!")
  if (!is.character(transf)) stop("ERROR: transf must be a character string!")
  if (!is.matrix(reference_profiles)) stop("ERROR: reference_profiles must be a matrix!")
  if (!is.numeric(min_cont)) stop("ERROR: min_cont must be numeric!")

  # Loading libraries
  suppressMessages(require(nnls))

  profile_mtrx <- topic_distribution_analysis(nmf_mod = nmf_mod,
                               mixture_transcriptome = mixture_transcriptome,
                               transf = transf,
                               regr = regr)

  # We add 1 extra column to add the residual error
  # 添加一个额外的列添加残差
  decon_mtrx <- matrix(data = NA,
                       nrow = ncol(profile_mtrx),
                       ncol = ncol(reference_profiles) + 1)
  colnames(decon_mtrx) <- c(colnames(reference_profiles), "res_ss")

  # create progress bar
  print("Deconvoluting spots")
  total <- ncol(profile_mtrx)
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  # reference_profiles = t(reference_profiles)
  if(regr_type == "nnls"){
    for (i in seq_len(ncol(profile_mtrx))) {
      # weights = profile_mtrx[, i]
      # print(i)
      # tmp = R_stepwise(y = profile_mtrx[, i],X = reference_profiles,sle=0.001, sls=0.002)
      # weights = pmax(tmp, 0)
      ## NNLS to get cell type composition
      ## NNLS获得cell type比例
      # tmp = bgd(y = profile_mtrx[, i], X = reference_profiles)
      # weights = pmax(tmp, 0)
      # print(Matrix::norm(profile_mtrx[, i]-weights))
      # weights <- pmax(weights, 0)
      # weights <- solve(reference_profiles,profile_mtrx[, i])
      # weights <- pls(reference_profiles,profile_mtrx[, i])
      nnls_pred <- nnls::nnls(A = reference_profiles, b = profile_mtrx[, i])

      weights <- nnls_pred$x
      # weights = profile_mtrx[, i]

      ## get proportions of each cell type
      comp <- weights / sum(weights)

      ## Remove cell types not contributing the minimum
      comp[comp < min_cont] <- 0

      ### Updated proportions after filtering out minimum contributions
      comp_prop <- comp / sum(comp)
      comp_prop[is.na(comp_prop)] <- 0

      ## Get Total sum of squares
      fit_null <- 0
      tot_ss <- sum((profile_mtrx[, i] - fit_null) ^ 2)

      ## Get % of unexplained residuals
      unexpl_ss <- nnls_pred$deviance / tot_ss
      # unexpl_ss <- arbor[[2]] / tot_ss
      decon_mtrx[i, 1:(ncol(decon_mtrx) - 1)] <- comp_prop
      decon_mtrx[i, ncol(decon_mtrx)] <- unexpl_ss

      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }else if(regr_type == "lar"){
    count = 0
    for (i in seq_len(ncol(profile_mtrx))) {
      object <- lars(reference_profiles,profile_mtrx[, i],type="lar", eps = 0)

      if(length(which.min(object$cp)) == 0){
        tmp = profile_mtrx[, i]
        count = count+1
      }else{
        tmp = object$beta[which.min(object$cp),]
      }

      # tmp = profile_mtrx[, i]
      weights = pmax(tmp, 0)

      ## get proportions of each cell type
      comp <- weights / sum(weights)

      ## Remove cell types not contributing the minimum
      comp[comp < min_cont] <- 0
      weights[comp < min_cont] <- 0

      ### Updated proportions after filtering out minimum contributions
      comp_prop <- comp / sum(comp)
      comp_prop[is.na(comp_prop)] <- 0


      unexpl_ss = 0
      # unexpl_ss <- arbor[[2]] / tot_ss
      decon_mtrx[i, 1:(ncol(decon_mtrx) - 1)] <- comp_prop
      decon_mtrx[i, ncol(decon_mtrx)] <- unexpl_ss
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    print(count)
  }
  # Close progress bar
  close(pb)

  # print('replace1')
  return(decon_mtrx = decon_mtrx)
}
