data_path = "/data/wangfx/data/my_deconvolution_spot/"
# data_name_list = c("HNSCC", "PBMC_Smart-Seq2")
# data_name_list = c("mouse_brain","HNSCC", "melanoma", "PBMC_Smart-Seq2",
#                    "PBMC_C1HT-medium", "PBMC_C1HT-small", "PBMC_CEL-Seq2", "PBMC_Chromium(sn)"
#                    , "PBMC_Chromium", "PBMC_ddSEQ", "PBMC_Drop-Seq", "PBMC_ICELL8"
#                    , "PBMC_inDrop", "PBMC_MARS-Seq", "PBMC_mcSCRB-Seq", "PBMC_Quartz-Seq2")

data_name_list = c("Young", "Muraro","HNSCC","melanoma")
# data_name_list = c("PBMC_Smart-Seq2")
hvg=3000; cl_n=100; min_cont=0; numiters=100; lam=100
sle = 1e-7; sls=sle*1.1;

# params
{
  # min_cont_list = c(0.01,0.03,0.05,0.07,0.09)
  min_cont_list = c(0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10)
}
