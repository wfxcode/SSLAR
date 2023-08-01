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

indrop_pdac_a <- readRDS(file = "R/analyses/data/PDAC/PDAC-A_itai_processed.RDS")

dfa_col <- indrop_pdac_a@meta.data%>%
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

st_se = readRDS("R/analyses/data/PDAC/PDAC-A_ST.RDS")
decon_mtrx = readRDS("R/analyses/data/PDAC/decon_mtrx_pdac_a.RDS")
st_se@meta.data = cbind(st_se@meta.data[1:6],decon_mtrx)
cell_types = dfa_col$df_names
tmp_paired_ct_plt <- lapply(cell_types, function(ct){
  # print(ct)

  tmp_plt <- plot_spaniel(data_df = data.frame(st_se@meta.data),
                          grob = st_se@images[[1]],
                          x = "x",
                          y = "y",
                          point_colour = ct,
                          point_size = ct) +
    ggplot2::theme_void() +
    # ggplot2::labs(title = sprintf("Proportion of %s",
    #                               dfa_col[dfa_col$df_name == ct, "plt_names"])) +
    ggplot2::labs(title = dfa_col[dfa_col$df_name == ct, "plt_names"]) +
    theme(plot.title =
            element_text(hjust = 0.1, size = 18,color = "black"),
          plot.margin=unit(c(-0.1,-0.5,-0.1,-0.5), "cm"),
          legend.position="none")
  # print(dfa_col[dfa_col$df_name == ct, "plt_names"])
  return(tmp_plt)
})
tmp_paired_ct_plt[[19]]+labs(title = NULL)
paired_ct_plt = ggpubr::ggarrange(plotlist = tmp_paired_ct_plt,
                  ncol = 5,
                  nrow = 4,
                  align = "hv")
save_plot("R/analyses/img/PDAC/all_types/all_types.pdf",paired_ct_plt,bg="white",
          base_width = 24,
          base_height = 18)
paired_ct_plt

i=1
for (ct in cell_types) {
  path = paste( "R/analyses/img/PDAC/all_types/",ct,".jpg", sep = "")
  print(path)
  save_plot(path,tmp_paired_ct_plt[[i]]+labs(title = NULL),bg="white",
            base_width = 4,
            base_height = 4)
  i = i+1
}

