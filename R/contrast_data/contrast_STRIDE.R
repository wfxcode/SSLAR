{
  source("R/import.R")
  # source("R/utils.R")
  # change 0

  source("R/contrast_data/contrast_params.R")
  data_name_list = c("PBMC_CEL-Seq2","PBMC_Chromium(sn)",
        "PBMC_ddSEQ","PBMC_Drop-Seq","PBMC_inDrop",
        "PBMC_MARS-Seq","PBMC_Quartz-Seq2","PBMC_Smart-Seq2","Young")
  {
  if(exists("params")) temp = params
  else rm(temp)

  model="STRIDE"
  params = list()

  for(i in seq_len(length(data_name_list))){


    save_path = "/data/wangfx/test_result/SPOTlight/contrast/"

    model="STRIDE2"
    data_name = data_name_list[[i]]
    save_path = paste(save_path, data_name, "/",model,"/", sep = "")




    save_file = get_save_files(path = save_path)


    data_files = get_data_files(path = data_path, data_name = data_name)

    out_dir = paste("Result/",data_name,sep = "")
    model="STRIDE"
    params[[i]] = list(model = model,data_name=data_name,out_dir=out_dir,
                       data_files = data_files, save_files = save_file)

  }


  if(exists("temp")){
    params = c(temp,params)
  }
  }
}
# contrast(params[[1]])
{
  # memory.limit() / memory.size()
  cl <- makeCluster(10)

  parLapply(cl, params, contrast)

  stopCluster(cl)
  # for (sl in sls) {
  #   params$sl = sl
  #   # print(paste(params,sep = " "))
  #   print(sprintf(format, params, params, params, params, params))
  # }
}
contrast(params = params[[1]])
for(i in seq_len(length(params))){
  # if(i == 1)
  #   next
  contrast(params = params[[i]])
}
