source("R/analyses/utils/bin.r")


brain <- readRDS("R/analyses/data/mouse_brain/brain1_processed.RDS")
decon_mtrx <- readRDS(file = "R/analyses/data/mouse_brain/decon_mtrx.RDS")
cell_types <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
colnames(decon_mtrx) <- c(cell_types, "res_ss")
brain@meta.data <- cbind(brain@meta.data, decon_mtrx)
cell_types_metadata <- colnames(brain@meta.data)[colnames(brain@meta.data) %in% cell_types]
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


ct_dict <- list()
ct_dict[["CA1.ProS"]][["gene"]] <- "Fibcd1"
ct_dict[["CA1.ProS"]][["name"]] <- "CA1.ProS"
ct_dict[["CA1.ProS"]][["plot_name"]] <- "Cornu Ammonis 1"

ct_dict[["CA2.IG.FC"]][["gene"]] <- "Ccdc3"
ct_dict[["CA2.IG.FC"]][["name"]] <- "CA2.IG.FC"
ct_dict[["CA2.IG.FC"]][["plot_name"]] <- "Cornu Ammonis 2"

ct_dict[["CA3"]][["gene"]] <- "Pvrl3"
ct_dict[["CA3"]][["name"]] <- "CA3"
ct_dict[["CA3"]][["plot_name"]] <- "Cornu Ammonis 3"

ct_dict[["DG"]][["gene"]] <- "Prox1"
ct_dict[["DG"]][["name"]] <- "DG"
ct_dict[["DG"]][["plot_name"]] <- "Dentate Gyrus"

legend_theme <-  theme(legend.text = element_text(colour = "#3b3a39", size = 8),
                       legend.title = element_text(colour = "#3b3a39", vjust = 1))

ct_ls <- c("CA1.ProS", "CA2.IG.FC", "CA3", "DG")
plot_ls <- lapply(ct_ls, function(ct) {
  feat <- ct_dict[[ct]][["gene"]]
  nam <- ct_dict[[ct]][["name"]]

  # If a gene doesn't exist add its row as all 0
  if (!feat %in% rownames(brain@assays$SCT@data)) {
    feat_mtrx <- Matrix::Matrix(0, nrow = 1, ncol = ncol(brain@assays$SCT@data), sparse = TRUE)
    rownames(feat_mtrx) <- feat
    brain@assays$SCT@var.features <- c(brain@assays$SCT@var.features, feat)
    brain@assays$SCT@data <- rbind(brain@assays$SCT@data, feat_mtrx)
  }
  # Plot cell type
  # browser()
  suppressMessages(ct_arr <- join_spatial_plots_2(spatial_obj = brain, feat = ct, nam = nam))
  legend_ct <- get_legend(ct_arr[[1]] + legend_theme)

  # Plot gene marker for that cell type
  suppressMessages(feat_arr <- join_spatial_plots_2(spatial_obj = brain, feat = feat, nam = nam))
  legend_feat <- get_legend(feat_arr[[1]] + legend_theme)

  ### Load image to plot ###
  fn_all <- list.files("analysis/mouse_brain/img")
  fn_ls <- fn_all[grepl(pattern = paste0("^", nam), x = fn_all)]

  fn_anterior <- fn_ls[grepl(pattern = "anterior_crop", x = fn_ls)]
  anterior_ish <- plot_image(img_path = sprintf("R/analyses/data/mouse_brain/img/%s_left.jpg", feat)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    coord_fixed(ratio = 0.925, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

  fn_posterior <- fn_ls[grepl(pattern = "posterior_crop", x = fn_ls)]
  posterior_ish <- plot_image(img_path = sprintf("R/analyses/data/mouse_brain/img/%s_right.jpg", feat)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    coord_fixed(ratio = 0.925, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

  ct_plt <- cowplot::plot_grid(plotlist = list(ct_arr[[1]] + theme(legend.position = "none"),
                                               ct_arr[[2]] + theme(legend.position = "none")),
                               align = "vh",
                               # axis = "trbl",
                               ncol = 2,
                               nrow = 1)

  feat_plt <- cowplot::plot_grid(plotlist = list(feat_arr[[1]] + theme(legend.position = "none"),
                                                 feat_arr[[2]] + theme(legend.position = "none")),
                                 align = "vh",
                                 # axis = "trbl",
                                 ncol = 2,
                                 nrow = 1)
  ish_plt <- cowplot::plot_grid(plotlist = list(anterior_ish,
                                                posterior_ish),
                                align = "vh",
                                # axis = "trbl",
                                ncol = 2,
                                nrow = 1)

  plt_arr <- cowplot::plot_grid(plotlist = list(ct_plt,
                                                feat_plt,
                                                NULL,
                                                ish_plt),
                                nrow = 4,
                                ncol = 1,
                                # align = "vh",
                                axis = "tblr",
                                rel_heights = c(1, 1, 0.075, 1))

  leg_arr <- cowplot::plot_grid(plotlist = list(legend_ct,
                                                legend_feat,
                                                NULL),
                                nrow = 3,
                                ncol = 1)

  final_arr <- cowplot::plot_grid(plotlist = list(plt_arr, leg_arr),
                                  nrow = 1,
                                  ncol = 2,
                                  rel_widths = c(1 ,0.2),
                                  axis = "trbl") %>%
    ggpubr::annotate_figure(p = .,
                            top = ggpubr::text_grob(sprintf("%s", ct_dict[[ct]][["plot_name"]]),
                                                    face = "bold",
                                                    size = 20, vjust = 1))


  #### Save plot ####
  ggpubr::ggexport(plotlist = list(final_arr),
                   filename = sprintf("R/analyses/img/mouse/Fig_%s_arrange.jpeg", nam),
                   width = 3000,
                   height = 3700,
                   res = 300)

  return(final_arr + theme(plot.background = element_rect(fill = NA, colour = "black", size = 2)))
})
