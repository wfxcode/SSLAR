# data_name_list = c("mouse_brain","Young","Muraro","PBMC_Smart-Seq2","HNSCC", "melanoma",
#                    "PBMC_C1HT-medium", "PBMC_C1HT-small", "PBMC_CEL-Seq2", "PBMC_Chromium(sn)"
#                    , "PBMC_Chromium", "PBMC_ddSEQ", "PBMC_Drop-Seq", "PBMC_ICELL8"
#                    , "PBMC_inDrop", "PBMC_MARS-Seq", "PBMC_mcSCRB-Seq", "PBMC_Quartz-Seq2")
data_name = "PBMC_Smart-Seq2"
data_path = "/data/wangfx/data/my_deconvolution_spot/"
hvg=3000; cl_n=100;
numiters=100; lam=100
sle = 1e-7; sls=sle*2;
ajust_path = "/data/wangfx/test_result/SPOTlight/contrast/"

# depth_list = c(1000,2500,5000,7500,10000,15000,20000,30000)
depth_list = c(1000,5000,10000,20000,30000)
