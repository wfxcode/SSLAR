data_path = "/data/wangfx/data/my_deconvolution_spot/"
hvg=6000; cl_n="all"; min_cont=0; numiters=100; lam=100
nmf = "ssNMF"; regr = "lar"; regr_type = "nnls"
ajust_path = "/data/wangfx/test_result/SPOTlight/adjust_params/"
# data_name_list = c("mouse_brain","PBMC_Smart-Seq2","Young","Muraro")
# data_name_list = c("PBMC_Quartz-Seq2")
data_name_list = c("PBMC_Smart-Seq2")

# data_name_list = c("PBMC_Chromium","PBMC_CEL-Seq2","PBMC_inDrop","PBMC_ddSEQ")
# data_name_list = c("PBMC_C1HT-small","PBMC_Chromium(sn)","PBMC_Drop-Seq","PBMC_ddSEQ")
sle = 0; sls=0;

# params
{
  depth_list = c(1000,2500,5000,7500,10000,12500,15000,17500,20000,30000,40000)
  # 0 is only marker gene   "all" is all gene
  hvg_list = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,9999999)
  # hvg_list = c(21000,22000,23000,24000,25000,26000,27000,28000,29000,30000)
  # cl_n_list = c(1,50,100,150,200,250,300,350,400,450,500,9999999)
  cl_n_list = c(1,20,40,60,80,100,120,140,160,180,200,9999999)
  numiters_list = c(1)
  lam_list = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  # lam_list = c(1,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,220,240,260,280,300)
  # lam_list = c(150,200,250,300,350,400,450,500)
  # lam_list = c(1,100,200,300,400,500,600,700,800,900)
  # options(scipen=200, digits=3)
  # sle2_list = c(1e-10,1e-09,1e-08,1e-07,1e-06,1e-05,1e-04,0.001,0.01,0.1)
  # sle1.5_list = c(1e-10,1e-09,1e-08,1e-07,1e-06,1e-05,1e-04,0.001,0.01,0.1)
  # sle1.1_list = c(1e-10,1e-09,1e-08,1e-07,1e-06,1e-05,1e-04,0.001,0.01,0.1)
}
