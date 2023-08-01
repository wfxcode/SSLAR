# 加载python
# Sys.setenv(RETICULATE_PYTHON = "/data/anaconda3/envs/wfx38/bin/python")

{
  library(reticulate)
  # py_module_available("tables")
  # Choose your python environment
  use_condaenv(condaenv = 'tensor36', required = TRUE)
  # use_python("D:/path/Anaconda3/envs/nmf38/python.exe", required = TRUE)
  py_config() # 安装的python版本环境查看，显示anaconda和numpy的详细信息
  source_python('R/method/nmf_method.py')
  # source('R/R_STRIDE.R')
}
# 加载包
# library(Matrix)
# library(data.table)
library(Seurat)
library(SeuratData)
# library(dplyr)
# library(gt)
library(SPOTlight)

# library(igraph)
# library(RColorBrewer)
# library(lmtest)
# library(modEvA)
# library(magrittr)
# library(parallel)
library(lars)
# library(RCTD)
# library(stringr)
# # install_julia()
# 加载julia
# library(JuliaCall)
# julia <- julia_setup('/home/wangfx/.julia/pythoncall/julia-1.6.4/bin')
# 更新函数
source("R/method/topic_identification.R")
source("R/method/cell_type_annotation.R")
source("R/method/topic_distribution_analysis.R")
source("R/method/test_synthetic_mixtures.R")
source("R/util/test_performances.R")
# source("R/util/utils.R")
# source("R/util/contrast_fun.R")
# source("R/My.stepwise.R")
# source("R/downsample_se_obj_fun.R")
# 加载函数

# options(java.parameters = "-Xmx15000m")
# library(xlsx)

# use_fake_coords = function (counts)
# {
#   coords <- data.frame(as.matrix(Matrix(0, nrow = dim(counts)[2],
#                                         ncol = 2)))
#   colnames(coords) <- c("x", "y")
#   rownames(coords) <- colnames(counts)
#   return(coords)
# }
# source_python('R/gradient_descent.py')

# init = function(){
#   # 加载python
#   Sys.setenv(RETICULATE_PYTHON = "/opt/anaconda3/envs/wfx38/bin/python")
#   library(reticulate)
#   use_condaenv(condaenv = 'wfx38', required = TRUE)
#   py_config() #安装的python版本环境查看，显示anaconda和numpy的详细信息
#
#   # # install_julia()
#   # 加载julia
#   # library(JuliaCall)
#   # julia <- julia_setup('/home/wangfx/.julia/pythoncall/julia-1.6.4/bin')
# }
# renew = function(){
#   # 更新函数
#   source("R/train_nmf1.R")
#   source("R/mixture_deconvolution_nmf_fun1.R")
#   source("R/spatial_decon_syn_assessment_nmf_fun1.R")
#   source("R/predict_spatial_mixtures_nmf_fun1.R")
#
#   # 加载函数
#   source_python('R/nmf_method.py')
#   source_python('R/r_mbpls.py')
# }
