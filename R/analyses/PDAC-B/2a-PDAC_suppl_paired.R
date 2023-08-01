library(dplyr)
library(NMF)
library(purrr)
library(tibble)
library(ggplot2)
library(Matrix)
library(Seurat)
library(Spaniel)
library(SPOTlight)
library(ggpubr)
# library(flextable)
source("R/analyses/utils/bin.r")
source("R/analyses/utils/spatial_plot_spaniel.R")

library(RColorBrewer)

indrop_pdac_b <- readRDS(file = "R/analyses/data/PDAC/PDAC-B_itai_processed.RDS")

dfa_col <- indrop_pdac_b@meta.data%>%
  dplyr::select(annotation,
                col_ct) %>%
  dplyr::distinct() %>%
  dplyr::arrange(annotation) %>%
  dplyr::mutate(annotation_mod = stringr::str_replace_all(
    string = annotation,
    pattern = "[[:punct:]]|[[:space:]]",
    replacement = ".")) %>%
  dplyr::rename(col_ct = col_ct,
                plt_names = annotation,
                df_names = annotation_mod)

st_se = readRDS("R/analyses/data/PDAC/PDAC-B_ST.RDS")
decon_mtrx = readRDS("R/analyses/data/PDAC/decon_mtrx_pdac_b.RDS")
st_se@meta.data = cbind(st_se@meta.data[1:6],decon_mtrx)
cell_types = dfa_col$df_names
# min(st_se@meta.data$x)
# max(st_se@meta.data$x)
# st_se@meta.data$x = st_se@meta.data$x-10
# source("/data/wangfx/R/ST/SSLAR/R/analyses/utils/spatial_plot_spaniel.R")
tmp_paired_ct_plt <- lapply(cell_types, function(ct){
  # print(ct)
  # browser()
  tmp_plt <- plot_spaniel_B(data_df = data.frame(st_se@meta.data),
                          grob = st_se@images[[1]],
                          x = "x",
                          y = "y",
                          point_colour = ct,
                          point_size = ct) +
    ggplot2::theme_void() +
    ggplot2::labs(title = dfa_col[dfa_col$df_name == ct, "plt_names"])+
    # theme(plot.title = element_text(vjust = 10, size = 10,color = "black"))
  theme(plot.title =
          element_text(hjust = 0.1, size = 14,color = "black"),
        plot.margin=unit(c(-0.1,-0.5,-0.1,-0.5), "cm"),
        legend.position="none")
  return(tmp_plt)
})
# tmp_paired_ct_plt[[3]]
paired_ct_plt = ggpubr::ggarrange(plotlist = tmp_paired_ct_plt,
                                  ncol = 5,
                                  nrow = 3,
                                  align = "hv")
paired_ct_plt
save_plot("R/analyses/img/PDAC-B/all_types/all_types.pdf",paired_ct_plt,bg="white",
          base_width = 21,
          base_height = 12)




i=1
for (ct in cell_types) {
  path = paste( "R/analyses/img/PDAC-B/all_types/",ct,".jpg", sep = "")
  print(path)
  save_plot(path,tmp_paired_ct_plt[[i]],bg="white",
            base_width = 4.3,
            base_height = 4)
  i = i+1
}

