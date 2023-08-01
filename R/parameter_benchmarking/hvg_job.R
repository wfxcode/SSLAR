source("R/import.R")
# source("R/utils.R")

source("R/parameter_benchmarking/params.R")
{
  cl_n=9999999;  min_cont=0.07; numiters=100; lam=100; nmf = "ssNMF"; regr = "lar"; regr_type = "nnls"
}
for(data_name in data_name_list){

  ajust_path = "/data/wangfx/test_result/SPOTlight/adjust_params/"

  # ajust_path = "/data/wangfx/test_result/SPOTlight/adjust_params/"


  ajust_path = paste(ajust_path, data_name, "/", sep = "")

  if (!dir.exists(ajust_path))
    dir.create(ajust_path)

  # change 1 and 1.1
  format = "hvg=%s/"
  folder = sprintf(format, substr(paste(hvg_list, collapse = "_"), 1 , 20))

  save_path = paste(ajust_path, folder, sep = "")

  # change 2
  save_files = get_save_files(path = save_path, params_array = hvg_list)
  data_files = get_data_files(path = data_path, data_name = data_name)

  {
    if(exists("params")) temp = params
    else rm(temp)

    params = list()
    # change 3
    for(i in seq_len(length(hvg_list))){
      # change 4 and 4.1
      hvg = hvg_list[i]
      save_file = save_files[[i]]
      params[[i]] = list(data_name=data_name, hvg=hvg, cl_n=cl_n, min_cont=min_cont,
                         numiters=numiters, lam = lam, sle=sle, sls=sls, nmf=nmf, regr=regr, regr_type=regr_type,
                         data_files = data_files, save_files = save_file)
      # contrast(params[[i]])
      # break
    }
    if(exists("temp")){
      params = c(temp,params)
    }
  }
}
{

  # memory.limit() / memory.size()
  cl <- makeCluster(12)

  parLapply(cl, params, contrast)

  stopCluster(cl)
  # for (sl in sls) {
  #   params$sl = sl
  #   # print(paste(params,sep = " "))
  #   print(sprintf(format, params, params, params, params, params))
  # }
}
