{
  {
  source("R/import.R")
  source("R/contrast_depth/depth_params.R")

  cl_n=100; hvg=3000; min_cont=0.09; numiters=100; lam=100;
  nmf = "nsNMF"; regr = "nnls"; regr_type = "nnls"; model="spotlight"
  data_name = "PBMC_Smart-Seq2"

  if(exists("params"))
    temp = params
  else rm(temp)
    params = list()
}
  {
  # change 1
  for(i in seq_len(length(depth_list))){

    save_path = "/data/wangfx/test_result/SPOTlight/depth/"



    save_path = paste(save_path, data_name, "/",model,"/", sep = "")

    if (!dir.exists(ajust_path))
      dir.create(ajust_path)


    # change 2
    save_files = get_save_files(path = save_path, params_array = depth_list)
    data_files = get_depth_data_files(path = data_path, data_name = data_name,depth = depth_list[[i]])

    {
      # change 3 and 4
      save_file = save_files[[i]]
      params[[i]] = list(data_name=data_name, hvg=hvg, cl_n=cl_n, min_cont=min_cont,
                         numiters=numiters, lam = lam, sle=sle, sls=sls, nmf=nmf, regr=regr, regr_type=regr_type,
                         data_files = data_files, save_files = save_file)
    }
  }

  if(exists("temp"))
    params = c(temp,params)
}
}
{
  # memory.limit() / memory.size()
  cl <- makeCluster(15)

  parLapply(cl, params, contrast)

  stopCluster(cl)
  # for (sl in sls) {
  #   params$sl = sl
  #   # print(paste(params,sep = " "))
  #   print(sprintf(format, params, params, params, params, params))
  # }
}
