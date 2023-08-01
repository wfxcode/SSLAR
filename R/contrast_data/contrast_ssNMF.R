{
  source("R/import.R")
  # source("R/utils.R")
  # change 0

  source("R/contrast_data/contrast_params.R")

  # data_path = "/data/wangfx/data/my_deconvolution_spot/"
  # data_name="PBMC_Quartz-Seq2"; hvg=3000; cl_n=100; min_cont=0.05; numiters=120; lam=200
  # sle = 0.0001; sls=0.0002; nmf = "ssNMF"; regr = "stepwise"; regr_type = "nnls"
  # ajust_path = "/data/wangfx/test_result/SPOTlight/adjust_params/"
  {
  if(exists("params")) temp = params
  else rm(temp)
  hvg="all"; cl_n="all"; min_cont=0.07
  nmf = "ssNMF"; regr = "nnls"; regr_type = "nnls"; model="ssNMF"

  params = list()

  for(i in seq_len(length(data_name_list))){


    save_path = "/data/wangfx/test_result/SPOTlight/contrast/"


    data_name = data_name_list[[i]]
    save_path = paste(save_path, data_name, "/",model,"/", sep = "")




    save_file = get_save_files(path = save_path)


    data_files = get_data_files(path = data_path, data_name = data_name)


    params[[i]] = list(data_name=data_name, hvg=hvg, cl_n=cl_n, min_cont=min_cont,
                       numiters=numiters, lam = lam, sle=sle, sls=sls, nmf=nmf, regr=regr, regr_type=regr_type,
                       data_files = data_files, save_files = save_file)
    # contrast(params[[i]])
    # break

    #   params = list(data_name=data_name,  hvg=hvg, cl_n=cl_n, min_cont=min_conts[1],
    #                 numiters=numiters, sl=sl, nmf=nmf, regr=regr, regr_type=regr_type)
    #   save_files = save_files[[1]]
    #   # contrast(params = params, data_path = data_path, save_files = save_files[[1]])

  }


  if(exists("temp")){
    params = c(temp,params)
  }
  }
}

contrast(params[[1]])
{
  # memory.limit() / memory.size()
  cl <- makeCluster(5)

  parLapply(cl, params, contrast)

  stopCluster(cl)
  # for (sl in sls) {
  #   params$sl = sl
  #   # print(paste(params,sep = " "))
  #   print(sprintf(format, params, params, params, params, params))
  # }
}
