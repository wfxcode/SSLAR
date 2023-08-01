## Libraries
{
  # source('R/import.R')
  library(ggplot2)
  library(purrr)
  library(cowplot)
  library(imager)
  library(scatterpie)
  library(RColorBrewer)
  # library(Cairo)
  # source("R/analyses/misc/paths_vrs.R")
  source("R/analyses/utils/bin.r")

  data_dir <- "R/analyses/data/mouse_brain/"
}

## Load data
brain <- readRDS("R/analyses/data/mouse_brain/brain1_processed.RDS")
decon_mtrx <- readRDS(file = "R/analyses/data/mouse_brain/decon_mtrx.RDS")

# Join the data so its easier to work with
cell_types <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

brain@meta.data <- cbind(brain@meta.data, decon_mtrx)
cell_types_metadata <- colnames(brain@meta.data)[colnames(brain@meta.data) %in% cell_types]

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
# Combine all L2 neurons
{
  brain_mod <- brain
  layer_cell_types = all_cell_types
  brain_mod@meta.data[["L2"]] <- rowSums(brain@meta.data[, c("L2.IT.ENTl", "L2.IT.ENTm")])
  brain_mod[["L2.IT.ENTl"]] <- NULL
  brain_mod[["L2.IT.ENTm"]] <- NULL


  # all_cell_types <- col_df$plt_names[col_df$plt_names != "L2" & !is.na(col_df$plt_names)]
  # cell_types_metadata_mod <- c(all_cell_types[!all_cell_types %in% c("L2 IT ENTl", "L2 IT ENTm")], "L2")

  # Combine all L2/3 neurons
  brain_mod@meta.data[["L2.3"]] <- rowSums(brain@meta.data[, c("L2.3.IT.CTX", "L2.3.IT.ENTl", "L2.3.IT.PPP", "L2.3.IT.RHP")])
  brain_mod[["L2.3.IT.CTX"]] <- NULL
  brain_mod[["L2.3.IT.ENTl"]] <- NULL
  brain_mod[["L2.3.IT.PPP"]] <- NULL
  brain_mod[["L2.3.IT.RHP"]] <- NULL

  # all_cell_types <- col_df$plt_names[col_df$plt_names != "L2/3" & !is.na(col_df$plt_names)]
  # cell_types_metadata_mod <- c(all_cell_types[!all_cell_types %in% c("L2/3 IT CTX", "L2/3 IT ENTl", "L2/3 IT PPP", "L2/3 IT RHP")], "L2/3")

  # Combine all L5 neurons
  brain_mod@meta.data[["L5"]] <- rowSums(brain@meta.data[, c("L5.IT.CTX", "L5.PPP", "L5.PT.CTX")])
  brain_mod[["L5.IT.CTX"]] <- NULL
  brain_mod[["L5.PPP"]] <- NULL
  brain_mod[["L5.PT.CTX"]] <- NULL

  # all_cell_types <- col_df$plt_names[col_df$plt_names != "L5" & !is.na(col_df$plt_names)]
  # cell_types_metadata_mod <- c(all_cell_types[!all_cell_types %in% c("L5 IT CTX", "L5 PPP", "L5 PT CTX")], "L5")

  # Combine all L5/6 neurons
  brain_mod@meta.data[["L5.6"]] <- rowSums(brain@meta.data[, c("L5.6.IT.TPE.ENT", "L5.6.NP.CTX")])
  brain_mod[["L5.6.IT.TPE.ENT"]] <- NULL
  brain_mod[["L5.6.NP.CTX"]] <- NULL

  # all_cell_types <- col_df$plt_names[col_df$plt_names != "L5" & !is.na(col_df$plt_names)]
  # cell_types_metadata_mod <- c(all_cell_types[!all_cell_types %in% c("L5/6 IT TPE-ENT", "L5/6 NP CTX")], "L5.6")

  # Combine all L6 neurons
  brain_mod@meta.data[["L6"]] <- rowSums(brain@meta.data[, c("L6.CT.CTX", "L6.IT.CTX", "L6.IT.ENTl","L6b.CTX", "L6b.CT.ENT")])
  brain_mod[["L6.CT.CTX"]] <- NULL
  brain_mod[["L6.IT.CTX"]] <- NULL
  brain_mod[["L6.IT.ENTl"]] <- NULL
  brain_mod[["L6b.CTX"]] <- NULL
  brain_mod[["L6b.CT.ENT"]] <- NULL


  # all_cell_types <- col_df$plt_names[col_df$plt_names != "L5" & !is.na(col_df$plt_names)]
  # cell_types_metadata_mod <- c(all_cell_types[!all_cell_types %in% c("L6 CT CTX", "L6 IT CTX", "L6 IT ENTl")], "L6")

  # Combine all L6b neurons
  # brain_mod@meta.data[["L6b"]] <- rowSums(brain@meta.data[, c("L6b.CTX", "L6b.CT.ENT")])
  # brain_mod[["L6b.CTX"]] <- NULL
  # brain_mod[["L6b.CT.ENT"]] <- NULL
}

{
  row = data.frame(col_vector=c( col_df[col_df$ct_names == "L2.IT.ENTl", "col_vector"],
                                 col_df[col_df$ct_names == "L5.IT.CTX", "col_vector"],
                                 col_df[col_df$ct_names == "L5.6.IT.TPE.ENT", "col_vector"],
                                 col_df[col_df$ct_names == "L6.CT.CTX", "col_vector"]),
                   plt_names=c("L2","L5","L5/6","L6"),
                   ct_names=c("L2","L5","L5.6","L6"))
  col_df_mod <- rbind(col_df, row)
  col_df_mod = col_df_mod %>%
    filter(!col_df_mod$ct_names %in% c("L2.IT.ENTl","L2.IT.ENTm","L2.3.IT.CTX","L2.3.IT.ENTl","L2.3.IT.PPP","L2.3.IT.RHP",
                               "L5.IT.CTX","L5.PPP","L5.PT.CTX","L5.6.IT.TPE.ENT","L5.6.NP.CTX","L6.CT.CTX",
                               "L6.IT.CTX","L6.IT.ENTl","L6b.CTX","L6b.CT.ENT"))

  layers = c("L2","L2.3","L3.IT.ENT","L4.RSP.ACA","L4.5.IT.CTX","L5","L5.6","L6")
  # brain_mod <- brain
  arrange <- lapply(layers, function(ct) {
    # print(ct)
    ct_ls <- join_spatial_plots(spatial_obj = brain_mod, ct = ct)
    plt_tmp <- cowplot::plot_grid(plotlist = list(ct_ls[[1]] + theme(legend.position = "none"),
                                                  NULL,
                                                  ct_ls[[2]] + theme(legend.position = "none")),
                                  nrow = 1,
                                  rel_widths = c(1.1,-0.1, 1.1),
                                  align = "hv"
                                  ,labels = col_df_mod[col_df_mod$ct_names == ct, "plt_names"])

    leg_grobtable <- get_legend(ct_ls[[1]])
    plt_titleless <- cowplot::plot_grid(plotlist = list(plt_tmp,NULL, leg_grobtable),
                       ncol = 3,
                       nrow = 1, rel_widths = c(1, 0.01,0.12)) +
      # theme(plot.background = element_rect(fill = NA, colour ="black", size = 2))
      theme(plot.background = element_rect(fill = "white"))

    return(plt_tmp)
  })

  i=1
  for (layer in layers) {
    path = paste( "R/analyses/img/mouse/layer-arrange/",layer,".jpg", sep = "")
    print(path)
    save_plot(path,arrange[[i]],bg="white",
              base_width = 10,
              base_height = 5)
    i = i+1
  }
}

{
  layers = c("L2","L2.3","L3.IT.ENT","L4.RSP.ACA","L4.5.IT.CTX","L5","L5.6","L6")
  row = data.frame(col_vector=c( col_df[col_df$ct_names == "L2.IT.ENTl", "col_vector"],
                                 col_df[col_df$ct_names == "L5.IT.CTX", "col_vector"],
                                 col_df[col_df$ct_names == "L5.6.IT.TPE.ENT", "col_vector"],
                                 col_df[col_df$ct_names == "L6.CT.CTX", "col_vector"]),
                   plt_names=c("L2","L5","L5/6","L6"),
                   ct_names=c("L2","L5","L5.6","L6"))
  col_df_mod <- rbind(col_df, row)
  col_df_mod = col_df_mod %>%
    filter(!col_df_mod$ct_names %in% c("L2.IT.ENTl","L2.IT.ENTm","L2.3.IT.CTX","L2.3.IT.ENTl","L2.3.IT.PPP","L2.3.IT.RHP",
                                       "L5.IT.CTX","L5.PPP","L5.PT.CTX","L5.6.IT.TPE.ENT","L5.6.NP.CTX","L6.CT.CTX",
                                       "L6.IT.CTX","L6.IT.ENTl","L6b.CTX","L6b.CT.ENT"))
  all_cell_types_mod = col_df$ct_names
  # layer <- c("L2","L2.3","L3.IT.ENT","L4.RSP.ACA","L4.5.IT.CTX","L5","L5.6","L6")
  layer <- c("L2.3")
  scatterpie_tmp = joint_scatterpie_fun(se_obj = brain_mod,
                                        cell_types_all = all_cell_types,
                                        img_path1 = "R/analyses/data/mouse_brain/anterior-spatial/tissue_lowres_image.png",
                                        img_path2 = "R/analyses/data/mouse_brain/posterior-spatial/tissue_lowres_image.png",
                                        slice1 = "anterior1",
                                        slice2 = "posterior1",
                                        cell_types_interest = layers,
                                        scatterpie_alpha = 1,
                                        pie_scale = 0.4,
                                        col_df = col_df
  )

  save_plot("R/analyses/img/mouse/L_layer.jpg",scatterpie_tmp,
            base_width = 10,
            base_height = 6)
}

{
  cortex_layers_mod <- c("L2.3")
  cell_types_metadata_mod <- c(all_cell_types[!all_cell_types %in% c("L2/3 IT.CTX","L2/3 IT.ENTl","L2/3 IT PPP","L2/3 IT RHP")], "L2.3")

  # general_cortex_plt <- joint_scatterpie_fun(se_obj = brain_mod,
  #                      cell_types_all = cell_types_metadata_mod,
  #                      img_path1 = "data/MusMusculus/sag_ant_1/spatial/tissue_lowres_image.png",
  #                      img_path2 = "data/MusMusculus/sag_post_1/spatial/tissue_lowres_image.png",
  #                      slice1 = "Anterior",
  #                      slice2 = "Posterior",
  #                      cell_types_interest = cortex_layers_mod,
  #                      return_legend = TRUE,
  #                      img_alpha = 1,
  #                      scatterpie_alpha = 1,
  #                      pie_scale = 0.2,
  #                      col_df = col_df)

  anterior_plt <- spatial_scatterpie(se_obj = brain_mod,
                                     cell_types_all = cell_types_metadata_mod,
                                     img_path = "R/analyses/data/mouse_brain/anterior-spatial/tissue_lowres_image.png",
                                     slice = "Anterior",
                                     cell_types_interest = cortex_layers_mod,
                                     pie_scale = 0.4) +
    scale_fill_manual(values = col_df[col_df$plt_names %in% cell_types_metadata_mod, "col_vector"],
                      breaks = cell_types_metadata_mod) +
    labs(fill = "") +
    coord_fixed(ratio = 1) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 5))

  posterior_plt <- spatial_scatterpie(se_obj = brain_mod,
                                      cell_types_all = cell_types_metadata_mod,
                                      img_path = "R/analyses/data/mouse_brain/posterior-spatial/tissue_lowres_image.png",
                                      slice = "Posterior",
                                      cell_types_interest = cortex_layers_mod,
                                      pie_scale = 0.4) +
    scale_fill_manual(values = col_df[col_df$plt_names %in% cell_types_metadata_mod, "col_vector"],
                      breaks = cell_types_metadata_mod) +
    labs(fill = "") +
    coord_fixed(ratio = 1) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 5))

  scatterpie_tmp <- cowplot::plot_grid(plotlist = list(anterior_plt[[1]] + theme(plot.margin = margin(0,0,0,0, "cm")),
                                                       posterior_plt[[1]] + theme(plot.margin = margin(0,0,0,0, "cm"))),
                                       ncol = 2,
                                       nrow = 1)


  general_cortex_plt <- cowplot::plot_grid(plotlist = list(anterior_plt,
                                                           posterior_plt + theme(legend.position = "none")),
                                           ncol = 2,
                                           nrow = 1,
                                           align = "hv",
                                           axis = "trbl") +
    coord_fixed()


  save_plot("R/analyses/img/mouse/L_layer.jpg",general_cortex_plt,
            base_width = 10,
            base_height = 6)
}

