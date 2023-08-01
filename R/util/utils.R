

get_dstg_data = function(data_name){
  # 空间转录组数据
  spots_data = paste("E:/XX/language/RLanguage/ST/SPOTlight-master/data/",data_name,"/DSTG_",data_name,"_",spot_number,"_data.rds",sep = "")
  # 单细胞数据
  sc_data_name = paste("E:/XX/language/RLanguage/ST/SPOTlight-master/data/",data_name,"/DSTG_",data_name,"_scRNAseq_data.rds",sep = "")

  sc.count <- as.matrix(readRDS(sc_data_name))
  st.count <- as.matrix(readRDS(spots_data))

  # 找到共有基因
  intersect.genes <- intersect(rownames(sc.count),rownames(st.count))
  sc.count <- sc.count[intersect.genes,]
  st.count <- st.count[intersect.genes,]
  count.list <- list(sc.count,st.count)

  # 单细胞数据标签
  sc_label_name = paste("E:/XX/language/RLanguage/ST/SPOTlight-master/data/",data_name,"/DSTG_",data_name,"_scRNAseq_label.rds",sep = "")
  label.list <- list(data.frame(readRDS(sc_label_name),stringsAsFactors=F))
  # 空间转录组数据标签
  spots_label = paste("E:/XX/language/RLanguage/ST/SPOTlight-master/data/",data_name,"/DSTG_",data_name,"_",spot_number,"_label.rds",sep = "")
  true <- readRDS(spots_label)

  return(list(count = count.list, sc_label = label.list, true = true))
}

evaluation = function(sheet, true , time, path, p){
  # 计算结果
  predict <- as.matrix(read.csv('/predict_output.csv',header=F))
  # 筛选spot中占比低于p的细胞类型
  predict[predict < p] = 0

  # 计算结果
  result = test_synthetic_performance1(
    test_spots_metadata_mtrx = as.matrix(true),
    spot_composition_mtrx = as.matrix(predict)
  )
  JSD_list = result[12]
  result[12] = NULL

  # 计算AUC与AUPR
  true_temp = as.numeric(true)
  predict = as.numeric(predict)
  for (i in seq_len(length(true_temp))){
    if(true_temp[i] != 0 && true_temp[i] != 1){
      predict[i] = predict[i] / true_temp[i]
      true_temp[i] = 1
    }
  }

  res_list = auc_aupr(predict, true_temp)
  auc_data = cbind(res_list[[1]],res_list[[2]])
  aupr_data = cbind(res_list[[3]],res_list[[4]])

  result = append(result,values = time,after = 0)
  # 将结果指标添加至excel

  JSD_file = paste(path, "JSD.xlsx",sep = "")
  auc_file = paste(path, "DSTG_auc.xlsx",sep = "")
  aupr_file = paste(path, "DSTG_aupr.xlsx",sep = "")

  if (sheet$getLastRowNum() != 0) {
    addDataFrame(result, sheet, row.names = F, col.names = F, startRow = sheet$getLastRowNum() + 2)
    write.xlsx(auc_data, file = auc_file, append = T, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",time,sep = ""))
    write.xlsx(aupr_data, file = aupr_file, append = T, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",time,sep = ""))
  } else {
    # 第一行
    addDataFrame(result,sheet,row.names = F)
    write.xlsx(auc_data, file = auc_file, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",time,sep = ""))
    write.xlsx(aupr_data, file = aupr_file, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",time,sep = ""))
  }



  return(list(sheet = sheet, auc_data = auc_data, aupr_data = aupr_data))
}

get_data_files = function(path, data_name){

  # 原始单细胞数据地址
  orgin_path = paste(path, "orgin/", data_name, ".rds", sep = "")

  path = paste(path, "contrast/", data_name, sep = "")
  # 在创建数据时，不存在就创建文件夹
  if (!dir.exists(path))
    dir.create(path, recursive=TRUE)

  # dstg single cell data adress
  dstg_sc_data_name = paste(path, "/DSTG_",data_name,"_scRNAseq_data.rds",sep = "")
  dstg_sc_label_name = paste(path, "/DSTG_",data_name,"_scRNAseq_label.rds",sep = "")

  dstg_sc_data = list(data = dstg_sc_data_name, label = dstg_sc_label_name)

  # spotlight single cell data adress
  spotlight_sc_data_name = paste(path, "/SPOTlight_", data_name, "_scRNAseq.rds", sep = "")
  spotlight_sc_marker_name = paste(path, "/SPOTlight_", data_name, "_scRNAseq_maker.rds", sep = "")
  spotlight_sc_data = list(data = spotlight_sc_data_name, maker = spotlight_sc_marker_name)

  # dstg and spotlight spacial data address

  dstg_spacial_data_names = list()
  spotlight_spacial_data_names = list()
  for(time in seq_len(10)){

    # dstg spacial data address
    dstg_spots_data = paste(path, "/DSTG_", data_name, "_",time,"_data.rds",sep = "")
    dstg_spots_label = paste(path, "/DSTG_", data_name, "_",time,"_label.rds",sep = "")

    dstg_spacial = list(data = dstg_spots_data, label = dstg_spots_label)

    dstg_spacial_data_names[[time]] = dstg_spacial

    # spotlight spacial data address
    spotlight_spacial = paste(path, "/SPOTlight_",data_name,"_",time,".rds",sep = "")

    spotlight_spacial_data_names[[time]] = spotlight_spacial
  }

  # for(time in seq_len(20)){
  #   print(dstg_spacial_data_names[[time]])
  #   print(spotlight_spacial_data_names[[time]])
  # }
  return(list(dstg_sc = dstg_sc_data, dstg_spacial = dstg_spacial_data_names
              , spotlight_sc = spotlight_sc_data, spotlight_spacial = spotlight_spacial_data_names
              ,orgin_sc = orgin_path))
}

# 不同测序深度数据
get_depth_data_files = function(path, data_name,depth){

  # 原始单细胞数据地址
  orgin_path = paste(path, "orgin/", data_name, ".rds", sep = "")

  spacial_path = paste(path, "contrast/",depth ," " ,data_name, sep = "")

  path = paste(path, "contrast/", data_name, sep = "")


  # 在创建数据时，不存在就创建文件夹
  if (!dir.exists(path))
    dir.create(path)
  if (!dir.exists(spacial_path))
    dir.create(spacial_path)
  # dstg single cell data adress
  dstg_sc_data_name = paste(path, "/DSTG_",data_name,"_scRNAseq_data.rds",sep = "")
  dstg_sc_label_name = paste(path, "/DSTG_",data_name,"_scRNAseq_label.rds",sep = "")

  dstg_sc_data = list(data = dstg_sc_data_name, label = dstg_sc_label_name)

  # spotlight single cell data adress
  spotlight_sc_data_name = paste(path, "/SPOTlight_", data_name, "_scRNAseq.rds", sep = "")
  spotlight_sc_marker_name = paste(path, "/SPOTlight_", data_name, "_scRNAseq_maker.rds", sep = "")
  spotlight_sc_data = list(data = spotlight_sc_data_name, maker = spotlight_sc_marker_name)

  # dstg and spotlight spacial data address

  dstg_spacial_data_names = list()
  spotlight_spacial_data_names = list()
  for(time in seq_len(10)){

    # dstg spacial data address
    dstg_spots_data = paste(spacial_path, "/DSTG_",depth ," " , data_name, "_",time,"_data.rds",sep = "")
    dstg_spots_label = paste(spacial_path, "/DSTG_",depth ," " , data_name, "_",time,"_label.rds",sep = "")

    dstg_spacial = list(data = dstg_spots_data, label = dstg_spots_label)

    dstg_spacial_data_names[[time]] = dstg_spacial

    # spotlight spacial data address
    spotlight_spacial = paste(spacial_path, "/SPOTlight_",depth ," " ,data_name,"_",time,".rds",sep = "")

    spotlight_spacial_data_names[[time]] = spotlight_spacial
  }

  # for(time in seq_len(20)){
  #   print(dstg_spacial_data_names[[time]])
  #   print(spotlight_spacial_data_names[[time]])
  # }
  return(list(dstg_sc = dstg_sc_data, dstg_spacial = dstg_spacial_data_names
              , spotlight_sc = spotlight_sc_data, spotlight_spacial = spotlight_spacial_data_names
              ,orgin_sc = orgin_path))
}

get_save_files = function(path, params_array = NULL){

  if (!dir.exists(path))
    dir.create(path, recursive=TRUE)
  # Record AUC AUPR and JSD


  address = list()
  # 没有参数列表时
  if(is.null(params_array)){

    address$indicator = paste(path, "indicator.xlsx", sep = "")
    address$rmse = paste(path, "rmse.xlsx", sep = "")
    address$jsd = paste(path, "jsd.xlsx", sep = "")
    address$pearson = paste(path, "pearson.xlsx", sep = "")

  # 参数列表
  }else{
    analysis = paste(path, "analysis/", sep = "")
    if (!dir.exists(analysis))
      dir.create(analysis, recursive=TRUE)
    for(i in seq_len(length(params_array))){
      temp = list()
      temp$indicator = paste(path, params_array[i],"-indicator.xlsx", sep = "")
      # temp$auc = paste(analysis, params_array[i], "-auc.xlsx", sep = "")
      # temp$aupr = paste(analysis, params_array[i], "-aupr.xlsx", sep = "")
      temp$rmse = paste(analysis, params_array[i], "-rmse.xlsx", sep = "")

      temp$jsd = paste(analysis, params_array[i], "-jsd.xlsx", sep = "")
      temp$pearson = paste(analysis, params_array[i], "-pearson.xlsx", sep = "")
      address[[i]] = temp
    }
  }

  return(address)
}

aver_stan = function(count){
  aver = list()
  stan = list()
  for(i in seq_len(length(count))){
    aver[[i]] = round(mean(count[[i]]),4)
    sd = round(sd(count[[i]]),4)
    stan[[i]] = paste(aver[[i]],"±",sd,sep="")
  }
  return(list(aver = aver, stan = stan))
}
