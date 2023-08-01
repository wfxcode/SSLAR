## Libraries
source('R/import.R')
library(ggplot2)
library(purrr)
library(cowplot)
library(imager)
library(scatterpie)

## Paths
tech <- "sc"
tissue <- "allen_ref_70k"
dwn_smplng <- "both"
org <- "mm"
# source("R/analyses/misc/paths_vrs.R")
source("R/analyses/utils/bin.r")

## Set common parameters

data_dir <- "R/analyses/data/mouse_brain/"
options(stringsAsFactors = FALSE)


library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector  <-  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_low <- "green"; col_high <- "red"

## Load data
brain <- readRDS("R/analyses/data/mouse_brain/brain1_processed.RDS")
decon_mtrx <- readRDS(file = "R/analyses/data/mouse_brain/decon_mtrx.RDS")

# Join the data so its easier to work with
cell_types <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

# cell_types <- stringr::str_replace_all(
#   string = cell_types,
#   pattern = "[[:punct:]]|[[:space:]]",
#   replacement = ".")

colnames(decon_mtrx) <- c(cell_types, "res_ss")

# decon_mtrx_prop <- round(decon_mtrx[, cell_types] / rowSums(decon_mtrx[, cell_types]), 4)
# decon_mtrx_prop[is.na(decon_mtrx_prop)] <- 0
brain@meta.data <- cbind(brain@meta.data, decon_mtrx)
cell_types_metadata <- colnames(brain@meta.data)[colnames(brain@meta.data) %in% cell_types]

# To maintain consistent colors for each cell type between different plots
# we will create a dataframe of equivalencies which we will draw colors from

# Let us load the metadata which contains the subclass_label and the Allen Brain Institute associated color
allen_ref_70k = readRDS("R/analyses/data/mouse_brain/allen_ref_70k")
metadata <- allen_ref_70k@meta.data


col_df <- metadata %>%
  dplyr::select(subclass_color,
                subclass_label) %>%
  dplyr::distinct() %>%
  dplyr::arrange(subclass_label) %>%
  dplyr::mutate(subclass_mod = stringr::str_replace_all(
    string = subclass_label,
    pattern = "[[:punct:]]|[[:space:]]",
    replacement = ".")) %>%
  dplyr::bind_rows(data.frame(subclass_color = "#005e16",
                              subclass_label = "L2/3",
                              subclass_mod = "L2.3")) %>%
  dplyr::rename(col_vector = subclass_color,
                plt_names = subclass_label,
                ct_names = subclass_mod)
col_df = col_df[rownames(col_df) != "SM-GE8ZH_S041_E1-50", ]

## Plots
### Compositions
# In this sections we are going to build compositions of
# a cell's spatial location + marker gene for that cell type in the tissue + ISH image of that cell type

#### Cell types of interest
# Astro        CA1-ProS       CA2-IG-FC
# 969            1704              19
# CA3            Car3              CR
# 322            1193              38
# CT SUB              DG            Endo
# 158            2474             197
# L2 IT ENTl      L2 IT ENTm     L2/3 IT CTX
# 191              42            6003
# L2/3 IT ENTl     L2/3 IT PPP     L2/3 IT RHP
# 344            1497             134
# L3 IT ENT      L4 RSP-ACA     L4/5 IT CTX
# 593             256           10854
# L5 IT CTX          L5 PPP       L5 PT CTX
# 3750              49            1750
# L5/6 IT TPE-ENT     L5/6 NP CTX       L6 CT CTX
# 353            2305            6350
# L6 IT CTX      L6 IT ENTl         L6b CTX
# 4658              78            2153
# L6b/CT ENT           Lamp5           Meis2
# 662            4805             132
# Micro-PVM          NP PPP          NP SUB
# 178             142             253
# Oligo           Pvalb        SMC-Peri
# 231            4109             133
# Sncg             Sst       Sst Chodl
# 1555            5511             282
# SUB-ProS             Vip            VLMC
# 471            6329             120


#### Scatterpie spatial plots
##### All cell types
all_cell_types <- col_df$plt_names[col_df$plt_names != "L2/3" & !is.na(col_df$plt_names)]

cnames <- data.frame(ct_names = colnames(brain@meta.data)) %>%
  dplyr::left_join(col_df, by = "ct_names") %>%
  dplyr::mutate(ct_names = dplyr::if_else(is.na(plt_names), ct_names, plt_names)) %>%
  dplyr::pull(ct_names)

colnames(brain@meta.data) <- cnames
anterior_plt <- spatial_scatterpie(se_obj = brain,
                                   cell_types_all = all_cell_types,
                                   img_path = "R/analyses/data/mouse_brain/anterior-spatial/tissue_lowres_image.png",
                                   slice = "anterior1",
                                   cell_types_interest = NULL,
                                   pie_scale = 0.4
) +
  scale_fill_manual(values = col_df[col_df$ct_names %in% cell_types_metadata, "col_vector"],
                    breaks = all_cell_types) +
  labs(fill = "") +
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 5))

posterior_plt <- spatial_scatterpie(se_obj = brain,
                                    cell_types_all = all_cell_types,
                                    img_path = "R/analyses/data/mouse_brain/posterior-spatial/tissue_lowres_image.png",
                                    slice = "posterior1",
                                    cell_types_interest = NULL,
                                    pie_scale = 0.4
) +
  scale_fill_manual(values = col_df[col_df$ct_names %in% cell_types_metadata, "col_vector"],
                    breaks = all_cell_types) +
  labs(fill = "") +
  coord_fixed(ratio = 1) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 5))


# cowplot::save_plot(plot = anterior_plt,
#                    filename = "R/analyses/img/mouse/Figure_YYY_spatial_scatterpie_allen_anterior.svg",
#                    base_width = 16,
#                    base_height = 12)
#
# cowplot::save_plot(plot = posterior_plt,
#                    filename = "R/analyses/img/mouse/Figure_YYY_spatial_scatterpie_allen_posterior.svg",
#                    base_width = 16,
#                    base_height = 12)




##### Brain Layer cells
# Next we want to plot neuronal layers and subtypes in both anterior and posterior sections
scatterpie_tmp = joint_scatterpie_fun(se_obj = brain,
                                      cell_types_all = all_cell_types,
                                      img_path1 = "R/analyses/data/mouse_brain/anterior-spatial/tissue_lowres_image.png",
                                      img_path2 = "R/analyses/data/mouse_brain/posterior-spatial/tissue_lowres_image.png",
                                      slice1 = "anterior1",
                                      slice2 = "posterior1",
                                      cell_types_interest = NULL,
                                      scatterpie_alpha = 1,
                                      pie_scale = 0.4,
                                      col_df = col_df
)

save_plot("R/analyses/img/mouse/allen_all.jpg",scatterpie_tmp,bg="white",
          base_width = 10,
          base_height = 6)
# L2/3
# L23_layer <- c("L2/3 IT CTX","L2/3 IT ENTl","L2/3 IT PPP","L2/3 IT RHP")
L_layer <- c("L2 IT ENTl","L2 IT ENTm","L2/3 IT CTX",
             "L2/3 IT ENTl","L2/3 IT PPP","L2/3 IT RHP",
             "L3 IT ENT","L4 RSP-ACA","L4/5 IT CTX",
             "L5 IT CTX","L5 PPP","L5 PT CTX",
             "L5/6 IT TPE-ENT","L5/6 NP CTX","L6 CT CTX",
             "L6 IT CTX","L6 IT ENTl","L6b CTX",
             "L6b/CT ENT")
scatterpie_tmp = joint_scatterpie_fun(se_obj = brain,
                                      cell_types_all = all_cell_types,
                                      img_path1 = "R/analyses/data/mouse_brain/anterior-spatial/tissue_lowres_image.png",
                                      img_path2 = "R/analyses/data/mouse_brain/posterior-spatial/tissue_lowres_image.png",
                                      slice1 = "anterior1",
                                      slice2 = "posterior1",
                                      cell_types_interest = L_layer,
                                      scatterpie_alpha = 1,
                                      pie_scale = 0.4,
                                      col_df = col_df
)

save_plot("R/analyses/img/mouse/L_layer.jpg",scatterpie_tmp,
          base_width = 10,
          base_height = 6)



for (layer in all_cell_types) {
  scatterpie_tmp = joint_scatterpie_fun(se_obj = brain,
                                        cell_types_all = all_cell_types,
                                        img_path1 = "R/analyses/data/mouse_brain/anterior-spatial/tissue_lowres_image.png",
                                        img_path2 = "R/analyses/data/mouse_brain/posterior-spatial/tissue_lowres_image.png",
                                        slice1 = "anterior1",
                                        slice2 = "posterior1",
                                        cell_types_interest = layer,
                                        scatterpie_alpha = 1,
                                        pie_scale = 0.4,
                                        col_df = col_df
  )
  layer = gsub("/","-",layer)
  # print(layer)
  path = paste( "R/analyses/img/mouse/all-types/",layer,".jpg", sep = "")
  save_plot(path,scatterpie_tmp,
            base_width = 10,
            base_height = 6)
}

{
  brain_mod = brain
  arrange <- lapply(all_cell_types, function(ct) {
    # print(ct)
    ct_ls <- join_spatial_plots(spatial_obj = brain_mod, ct = ct)
    plt_tmp <- cowplot::plot_grid(plotlist = list(ct_ls[[1]] + theme(legend.position = "none"),
                                                  NULL,
                                                  ct_ls[[2]] + theme(legend.position = "none")),
                                                  nrow = 1,
                                                  rel_widths = c(1.1,-0.1, 1.1),
                                                  align = "hv",
                                                  labels = ct)

    leg_grobtable <- get_legend(ct_ls[[1]])
    plt_titleless <- cowplot::plot_grid(plotlist = list(plt_tmp,NULL, leg_grobtable),
                       ncol = 3,
                       nrow = 1, rel_widths = c(1, 0.01,0.12)) +
      # theme(plot.background = element_rect(fill = NA, colour ="black", size = 2))
      theme(plot.background = element_rect(fill = "white"))

    return(plt_tmp)
  }) %>%
    ggarrange(plotlist = ., ncol = 4, align = "hv")
  arrange

  cowplot::save_plot(plot = arrange, bg="white",
                     filename = "R/analyses/img/mouse/all-types-arrange/all_types.jpg",
                     base_width = 24,
                     base_height = 24)

  ggpubr::ggexport(plotlist = arrange,bg="white",
                   filename = "R/analyses/img/mouse/all-types-arrange/all_types.pdf",
                   width = 12,
                   height = 24,
                   res = 300)


  i=1
  for (layer in all_cell_types) {
    layer = col_df[col_df$plt_names == layer, "ct_names"]
    path = paste( "R/analyses/img/mouse/all-types-arrange/",layer,".jpg", sep = "")
    save_plot(path,arrange[[i]],bg="white",
              base_width = 10,
              base_height = 5)
    i = i+1
  }
}
# cowplot::draw_plot(scatterpie_tmp)

# cowplot::save_plot(plot = scatterpie_tmp,
#                    filename = "R/analyses/img/mouse/Figure_YYY_spatial_scatterpie_allen_all.svg",
#                    base_width = 16,
#                    base_height = 12)

###### Individual proportion plots
layer = c("L2.3")
# brain_mod <- brain
arrange <- lapply(layer, function(ct) {
  # print(ct)
  ct_ls <- join_spatial_plots(spatial_obj = brain_mod, ct = ct)
  plt_tmp <- cowplot::plot_grid(plotlist = list(ct_ls[[1]] + theme(legend.position = "none"),
                                                ct_ls[[2]] + theme(legend.position = "none")),
                                nrow = 1,
                                ncol = 2,
                                align = "hv",
                                labels = col_df[col_df$ct_names == ct, "plt_names"])

  # leg_grobtable <- get_legend(ct_ls[[1]])
  # plt_titleless <- cowplot::plot_grid(plotlist = list(plt_tmp, leg_grobtable),
  #                    ncol = 2,
  #                    nrow = 1, rel_widths = c(1, 0.2)) +
  #   theme(plot.background = element_rect(fill = NA, colour ="black", size = 2))

  return(plt_tmp)
})

save_plot("R/analyses/img/mouse/layer-arrange/L2.3.jpg",arrange[[1]],
          base_width = 10,
          base_height = 6)



tmp = rownames(brain_mod@assays[["Spatial"]]@data)
which(grepl("vl3",tmp,ignore.case = TRUE))
