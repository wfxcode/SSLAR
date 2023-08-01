# 由参数跑20次实验
# params 参数列表
# data_path 数据地址
# save_files 保存数据地址
contrast = function(params){
  # library(Seurat)
  source("R/import.R")

  data_files = params$data_files
  save_files = params$save_files
  # data_files = get_data_files(path = data_path, data_name = params$data_name)


  # 保存指标
  wb <- xlsx::createWorkbook()
  sheet <- createSheet(wb,sheetName = paste("sheet1", sep = ""))
  # browser()
  for(i in seq_len(length(data_files$spotlight_spacial))){

    # 空间数据
    st_data = readRDS(data_files$spotlight_spacial[[i]])
    # 空间数据标签
    true <- st_data[[2]]
    # saveRDS(true,"R/debug/true.RDS")
    st_counts <- as.matrix(st_data[[1]])
    # break

    # 计时开始
    start_t <- Sys.time()
    if("model" %in% names(params)){
      if(params$model == "RCTD"){
        print("执行RCTD")
        pred_comp = RCTD_decon(params,data_files$dstg_spacial[[i]])
      }else if(params$model == "STRIDE"){
        print("执行STRIDE")
        pred_comp = STRIDE_decon(params,data_files$dstg_spacial[[i]],i)
        # browser()
      }else if(params$model == "SCDC"){
        print("执行SCDC")
        pred_comp = SCDC_decon(params,st_counts)
      }else if(params$model == "MuSiC"){
        print("执行MuSiC")
        pred_comp = MuSiC_decon(params,st_counts)
      }
    }else{
      # SSLAR类似方法
      pred_comp = my_decon(params = params,
                           st_counts = st_counts)
      # 筛选
      ct_cols <- colnames(pred_comp)[which(colnames(pred_comp) != "res_ss")]
      pred_comp = pred_comp[, ct_cols]
      pred_comp = pred_comp[,sort(colnames(pred_comp))]
    }


    # 时间
    total_t <- round(difftime(Sys.time(), start_t, units = "mins"), 2)
    # total_t <- round(difftime(Sys.time(), start_t, units = "secs"), 2)

    true = true[,sort(colnames(true))]
    # colnames(pred_comp) = colnames(true)
    # 得到指标
    result <- test_synthetic_performance1(test_spots_metadata_mtrx = as.matrix(true),
                                          spot_composition_mtrx = as.matrix(pred_comp))
    # 拆分指标并保存
    RMSE_list = result$RMSE_list
    JSD_list = result$JSD_list
    pearson_list = result$pearson_list
    result = result[1:8]

    result$time_min = as.numeric(total_t)


    if (sheet$getLastRowNum() != 0) {
      count = mapply(c, count, result, SIMPLIFY=FALSE)
      addDataFrame(result, sheet, row.names = F, col.names = F, startRow = sheet$getLastRowNum() + 2)
      write.xlsx(RMSE_list, file = save_files$rmse, append = T, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
      write.xlsx(JSD_list, file = save_files$jsd, append = T, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
      write.xlsx(pearson_list, file = save_files$pearson, append = T, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
      # write.xlsx(auc_data, file = save_files$auc, append = T, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
      # write.xlsx(aupr_data, file = save_files$aupr, append = T, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
    } else {
      count = result
      # 第一行
      addDataFrame(result,sheet,row.names = F)
      write.xlsx(RMSE_list, file = save_files$rmse, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
      write.xlsx(JSD_list, file = save_files$jsd, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
      write.xlsx(pearson_list, file = save_files$pearson, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
      # write.xlsx(auc_data, file = save_files$auc, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
      # write.xlsx(aupr_data, file = save_files$aupr, row.names = FALSE, col.names = FALSE, sheetName = paste("Sheet",i,sep = ""))
    }
    # if(i == 1)
    #   break
  }
  # 计算均值与标准差
  analysis = aver_stan(count = count)
  addDataFrame(analysis$aver, sheet, row.names = F, col.names = F, startRow = sheet$getLastRowNum() + 3)
  addDataFrame(analysis$stan, sheet, row.names = F, col.names = F, startRow = sheet$getLastRowNum() + 2)

  saveWorkbook(wb, file = save_files$indicator)
  print(save_files$indicator)
}

my_decon = function(params, st_counts){
  source("R/downsample_se_obj1.R")
  data_files = params$data_files
  # load sc data
  sc_data = readRDS(data_files$spotlight_sc$data)
  sc_maker = readRDS(data_files$spotlight_sc$maker)

  colnames(st_counts) <- paste("mixt", 1:ncol(st_counts), sep = "_")
  if(params$cl_n == 9999999 && params$hvg == 9999999){
    sc_down = sc_data
  }else{
    # Downsample number of genes and number of samples
    sc_down <- downsample_se_obj1(se_obj = sc_data,
                                  clust_vr = 'subclass',
                                  cluster_markers = sc_maker,
                                  cl_n = params$cl_n,
                                  hvg = params$hvg)
  }
  # Train the NMF model
  nmf_mod_ls <- topic_identification(cluster_markers = sc_maker,
                                     se_sc = sc_down,
                                     mtrx_spatial = st_counts,
                                     ntop = NULL,
                                     transf = "uv",
                                     clust_vr = 'subclass',
                                     nmf = params$nmf,
                                     numiters = params$numiters,
                                     lam = params$lam,
                                     assay = "RNA",
                                     slot = "counts")

  # 计算出矩阵Q，包含每个cell type的topic
  ct_topic_profiles <- topic_profile_per_cluster_nmf(h = coef(nmf_mod_ls[[1]]),
                                                     train_cell_clust = nmf_mod_ls[[2]])
  # 解卷积spot，通过NMF模型得到spot的cell type比例
  # source('R/import.R')

  pred_comp <- cell_type_annotation(nmf_mod = nmf_mod_ls[[1]],
                                    mixture_transcriptome = st_counts,
                                    transf = "uv",
                                    reference_profiles = ct_topic_profiles,
                                    min_cont = params$min_cont,
                                    regr = params$regr,
                                    regr_type = params$regr_type)
  return(pred_comp)
  # return(list(nmf_mod_ls, pred_comp))
}

RCTD_decon = function(params,spacial_path){
  data_files = params$data_files
  counts = readRDS(data_files$dstg_sc[["data"]])

  cell_info = readRDS(data_files$dstg_sc[["label"]])
  cell_types = cell_info[,"subclass"]
  names(cell_types) = row.names(cell_info)
  cell_types <- as.factor(cell_types) # convert to factor data type

  ### Create the Reference object
  ## 单细胞数据，nUMI可以不写，会自动计算
  # as(counts, "dgCMatrix")
  reference <- Reference(as(counts, "dgCMatrix"), cell_types)
  counts <- readRDS(spacial_path[["data"]])

  ### Create SpatialRNA object
  puck <- RCTD::SpatialRNA(NULL, as(counts, "dgCMatrix"), use_fake_coords = TRUE)

  myRCTD <- RCTD::create.RCTD(puck, reference, max_cores = 4, CELL_MIN_INSTANCE = 1)
  myRCTD <- RCTD::run.RCTD(myRCTD, doublet_mode = 'full')

  results <- myRCTD@results

  pred_comp <- as.matrix(sweep(results$weights, 1, rowSums(results$weights), '/'))
  return(pred_comp)
}

STRIDE_decon = function(params,spacial_path,i){
  data_files = params$data_files
  sc_count_file = data_files$dstg_sc[["data"]]
  sc_anno_file = data_files$dstg_sc[["label"]]
  st_count_file = spacial_path[["data"]]

  # out_dir = "Result/Muraro"
  # out_prefix = "Muraro"
  out_dir = params$out_dir
  out_prefix = params$data_name
  # browser()
  pred_comp = R_STRIDE(sc_count_file,sc_anno_file,st_count_file,out_dir,out_prefix)
  return(pred_comp)
}

SCDC_decon = function(params,st_counts){
  library(Biobase)
  data_files = params$data_files

  sc_data = readRDS(data_files$spotlight_sc$data)
  sc_data@meta.data$batch = "pseudo"

  sc_data[["barcode"]] <- colnames(sc_data)
  expr_sc <- Biobase::ExpressionSet(assayData = as.matrix(sc_data@assays$RNA@counts),
                                    phenoData = AnnotatedDataFrame(data = data.frame(sc_data@meta.data)))

  expr_st <- Biobase::ExpressionSet(assayData = as.matrix(st_counts))

  scdc_deconv <- SCDC::SCDC_prop_ONE(bulk.eset = expr_st,
                                     sc.eset = expr_sc,
                                     ct.varname = "subclass",
                                     sample = "batch",
                                     ct.sub = unique(expr_sc$subclass))

  pred_comp = scdc_deconv$prop.est.mvw

  pred_comp = pred_comp[,sort(colnames(pred_comp))]
  return(pred_comp)
}

MuSiC_decon = function(params,st_counts){
  library(Biobase)
  library(MuSiC)
  start_t <- Sys.time()
  data_files = params$data_files

  sc_data = readRDS(data_files$spotlight_sc$data)

  sc_data[["barcode"]] <- colnames(sc_data)
  expr_sc <- Biobase::ExpressionSet(assayData = as.matrix(sc_data@assays$RNA@counts),
                                    phenoData = Biobase::AnnotatedDataFrame(data = sc_data@meta.data))

  expr_st <- Biobase::ExpressionSet(assayData = as.matrix(st_counts))
  print("计算music")
  music_deconv <- MuSiC::music_prop(bulk.eset = expr_st,
                                    sc.eset = expr_sc,
                                    markers = NULL,
                                    clusters = "subclass",
                                    samples = "barcode",
                                    select.ct = NULL,
                                    verbose = FALSE)

  pred_comp <- music_deconv$Est.prop.weighted

  pred_comp = pred_comp[,sort(colnames(pred_comp))]
  total_t <- round(difftime(Sys.time(), start_t, units = "mins"), 2)
  print(total_t)
  return(pred_comp)
}
