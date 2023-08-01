source("R/import.R")
# source("R/utils.R")
source("R/parameter_benchmarking/min_cont/cont_params.R")
# No 2 time
{
  hvg = 3000; cl_n=100;
}
for(data_name in data_name_list){

  ajust_path = "/data/wangfx/test_result/SPOTlight/adjust_params/"
  nmf = "nsNMF"; regr = "lar"; regr_type = "nnls"
  ajust_path = paste(ajust_path, "/min_cont/", data_name,"/LARS", sep = "")

  # change 1 and 1.1
  format = "min_cont=%s/"
  folder = sprintf(format, substr(paste(min_cont_list, collapse = "_"), 1 , 20))

  save_path = paste(ajust_path, folder, sep = "")

  # change 2
  save_files = get_save_files(path = save_path, params_array = min_cont_list)
  data_files = get_data_files(path = data_path, data_name = data_name)

  {
    if(exists("params")) temp = params
    else rm(temp)

    params = list()
    # change 3
    for(i in seq_len(length(min_cont_list))){
      # change 4 and 5
      min_cont = min_cont_list[i]
      save_file = save_files[[i]]
      params[[i]] = list(data_name=data_name, hvg=hvg, cl_n=cl_n, min_cont=min_cont,
                         numiters=numiters, lam = lam, sle=sle, sls=sls, nmf=nmf, regr=regr, regr_type=regr_type,
                         data_files = data_files, save_files = save_file)
      # contrast(params[[i]])
      # break
    }
    #   params = list(data_name=data_name,  hvg=hvg, cl_n=cl_n, min_cont=min_conts[1],
    #                 numiters=numiters, sl=sl, nmf=nmf, regr=regr, regr_type=regr_type)
    #   save_files = save_files[[1]]
    #   # contrast(params = params, data_path = data_path, save_files = save_files[[1]])
    if(exists("temp")){
      params = c(temp,params)
    }

  }
}
{
  # memory.limit() / memory.size()
  cl <- makeCluster(15)
  # lams
  parLapply(cl, params, contrast)

  stopCluster(cl)
  # for (sl in sls) {
  #   params$sl = sl
  #   # print(paste(params,sep = " "))
  #   print(sprintf(format, params, params, params, params, params))
  # }
}
# params[[1]] = list(data_name=data_name, hvg=hvg, cl_n=cl_n, min_cont=min_cont,
#                    numiters=numiters, lam = lam, sle="0.00001", sls="0.00002", nmf=nmf, regr=regr, regr_type=regr_type,
#                    data_files = data_files, save_files = save_file)
# contrast(params[[1]])
