source("R/import.R")
# source("R/utils.R")
source("R/parameter_benchmarking/params.R")
{
  cl_n=20; hvg="all"; min_cont=0.07; numiters=100; lam=100; nmf = "ssNMF"; regr = "lar"; regr_type = "nnls"
}
data_name_list = c("PBMC_Smart-Seq2","PBMC_C1HT-medium", "PBMC_C1HT-small", "PBMC_CEL-Seq2",
                   "PBMC_Chromium(sn)", "PBMC_Chromium", "PBMC_ddSEQ", "PBMC_Drop-Seq", "PBMC_ICELL8",
                   "PBMC_inDrop", "PBMC_MARS-Seq", "PBMC_mcSCRB-Seq", "PBMC_Quartz-Seq2")


params = list()
for(i in seq(13)){
  data_name = data_name_list[[i]]
  ajust_path = "/data/wangfx/test_result/SPOTlight/adjust_params/"

  # ajust_path = "/data/wangfx/test_result/SPOTlight/adjust_params/"

  if (!dir.exists(ajust_path))
    dir.create(ajust_path)

  # change 1 and 1.1
  format = "protocol=%s/"
  folder = sprintf(format, substr(paste(data_name_list, collapse = "_"), 1 , 20))

  save_path = paste(ajust_path, folder, sep = "")

  # change 2
  save_files = get_save_files(path = save_path, params_array = data_name_list)
  data_files = get_data_files(path = data_path, data_name = data_name)

  {
    save_file = save_files[[i]]
    params[[i]] = list(data_name=data_name, hvg=hvg, cl_n=cl_n, min_cont=min_cont,
                       numiters=numiters, lam = lam, sle=sle, sls=sls, nmf=nmf, regr=regr, regr_type=regr_type,
                       data_files = data_files, save_files = save_file)
  }
}
{

  # memory.limit() / memory.size()
  cl <- makeCluster(13)

  parLapply(cl, params, contrast)

  stopCluster(cl)
  # for (sl in sls) {
  #   params$sl = sl
  #   # print(paste(params,sep = " "))
  #   print(sprintf(format, params, params, params, params, params))
  # }
}
