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

{
  decon_mtrx = readRDS("R/analyses/data/PDAC/decon_mtrx_ica_pdac_b.RDS")
  decon_mtrx_subs <- decon_mtrx[, colnames(decon_mtrx)[! colnames(decon_mtrx) %in% "res_ss"]]
  colnames(decon_mtrx_subs) <- gsub(pattern = "[[:punct:]]|[[:blank:]]",
                                    replacement = ".",
                                    x = colnames(decon_mtrx_subs),
                                    perl = TRUE)
  st_se = readRDS("R/analyses/data/PDAC/PDAC-B_ST.RDS")
  st_se@meta.data <- cbind(st_se@meta.data, decon_mtrx_subs)
  ## Preprocess data
  spatial_coord <- data.frame(st_se@meta.data) %>%
    tibble::rownames_to_column("ID")
  ica_se = readRDS("R/analyses/data/PDAC/ica_se_processed.RDS")
  ica_se_annot <- data.frame(
    df_names = sort(unique(as.character(ica_se$specific_cell_type))),
    plt_name = sort(unique(as.character(ica_se$cluster))),
    col_ct = c("#1f9747","#07b1c0","#0371ae","#f27d2f","#d91f2d","#8066a8","#814844","#787879","#d375ae","#a8b32e"))


  tmp_df <- ica_se_annot
  ind <- which(names(spatial_coord) %in% tmp_df$df_name)
  names(spatial_coord)[ind] <- tmp_df$plt_name
  ct_all <- names(spatial_coord)[names(spatial_coord) %in% tmp_df$plt_name]
  ind_rm <- which(colSums(spatial_coord[, ct_all] > 0) == 0)
  if (length(ind_rm) > 0) {
    ct_all <- ct_all[-ind_rm]
  }

  scatterpie_plt <- ggplot() +
    scatterpie::geom_scatterpie(data = spatial_coord,
                                aes(x = x,
                                    y = y),
                                cols = ct_all,
                                color = NA,
                                alpha = 1,
                                pie_scale = 1.4) +
    scale_y_reverse() +
    theme_half_open(11, rel_small = 1) +
    theme_void() +
    coord_fixed(ratio = 1) +
    scale_fill_manual(values = tmp_df[tmp_df$plt_name %in% ct_all, "col_ct"]) +
    labs(title = "PDAC-B Immune Spatial scatterpie",
         color = "Cell types") +
    theme(
      legend.position = 'none',
      # plot.background = element_rect(fill = "#FFFFFF"),
      # panel.background = element_blank(),
      # plot.margin = margin(20, 20, 20, 20),
      plot.title = element_text(hjust = 0.5, size = 20))

  save_plot("R/analyses/img/PDAC-B/PDAC_immune_all.jpg",scatterpie_plt,
            bg="white", base_width = 10,
            base_height = 9)

}
{
  cell_types = ica_se_annot$df_names

  ct_plt <- lapply(cell_types, function(ct){
    # print(ct)

    tmp_plt <- plot_spaniel_B(data_df = data.frame(st_se@meta.data),
                            grob = st_se@images[[1]],
                            x = "x",
                            y = "y",
                            point_colour = ct,
                            point_size = ct) +
      ggplot2::theme_void() +
      ggplot2::labs(title = tmp_df[tmp_df$df_name == ct, "plt_name"])+
      theme(plot.title =
              element_text(hjust = 0.1, size = 18,color = "black"),
            plot.margin=unit(c(-0.1,-0.5,-0.1,-0.5), "cm"),
            legend.position="none")
      # theme(plot.title = element_text(hjust = 0.5, size = 10,color = "black"))

    return(tmp_plt)
  })
  immune_ct_plt = ggpubr::ggarrange(plotlist = ct_plt,
                                    ncol = 5,
                                    nrow = 2,
                                    align = "hv")
  save_plot("R/analyses/img/PDAC-B/immune_all_types/all_types.pdf",immune_ct_plt,bg="white",
            base_width = 21,
            base_height = 8)
  i=1
  for (ct in cell_types) {
    path = paste( "R/analyses/img/PDAC-B/immune_all_types/",ct,".jpg", sep = "")
    print(path)
    save_plot(path,ct_plt[[i]],bg="white",
              base_width = 4.3,
              base_height = 4)
    i = i+1
  }
}
