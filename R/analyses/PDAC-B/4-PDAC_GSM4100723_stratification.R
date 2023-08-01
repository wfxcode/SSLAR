{
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
  library(cowplot)
  # library(flextable)
  source("R/analyses/utils/bin.r")
  source("R/analyses/utils/spatial_plot_spaniel.R")
}
st_se = readRDS("R/analyses/data/PDAC/PDAC-B_ST.RDS")
ica_se = readRDS("R/analyses/data/PDAC/ica_se_processed.RDS")

decon_mtrx_immune = readRDS("R/analyses/data/PDAC/decon_mtrx_ica_pdac_b.RDS")
decon_mtrx_paired = readRDS("R/analyses/data/PDAC/decon_mtrx_pdac_b.RDS")
colnames(decon_mtrx_paired)
cancer_prop <- decon_mtrx_paired[, c("Cancer.clone..TM4SF1.")]
healthy_prop <- rowSums(decon_mtrx_paired[, c("Acinar.cells", "Centroacinar.ductal.cells", "Terminal.ductal.cells", "Antigen.presenting.ductal.cells", "Tuft.cells" )])

# par(mfrow = c(2, 1))
st_se[["cancer_prop"]] <- cancer_prop
st_se[["healthy_prop"]] <- healthy_prop

{

  pt1 <- Spaniel::spanielPlot(object = st_se,
                              grob = st_se@images[[1]],
                              plotType = "Cluster",
                              clusterRes = "cancer_prop",
                              ptSize = 5) +
    theme_void() +
    coord_fixed(1) +
    scale_alpha(range = c(1, 1)) +
    scale_color_gradientn(
      colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous"))
  # colours = heat.colors(10, rev = TRUE))

  pt2 <- Spaniel::spanielPlot(object = st_se,
                              grob = st_se@images[[1]],
                              plotType = "Cluster",
                              clusterRes = "healthy_prop",
                              ptSize = 5) +
    theme_void() +
    coord_fixed(1) +
    scale_alpha(range = c(1, 1)) +
    scale_color_gradientn(
      colours = wesanderson::wes_palette("Zissou1", 100, type = "continuous"))
  # colours = heat.colors(10, rev = TRUE))

  # pt3 <- ggplot() +
  #   geom_histogram(aes(x = cancer_prop), binwidth = 0.01) +
  #   geom_vline(xintercept = 0.25, col = "red")

  ggpubr::ggarrange(pt1, pt2, ncol = 2)

}
{
  threshold = 0.4
  # 通过密度图进行可视化，将肿瘤细胞和非肿瘤细胞的准确鉴别阈值设定在0.4。
  st_se[["status_2_territories"]] <- if_else(st_se$cancer_prop > threshold, "Tumoral", "Non-Tumoral")

  data_df <- data.frame(st_se@meta.data)
  data_df$y_inv <- 36 - data_df$y

  tmp_plt_2 <- ggplot(data_df,
                      ggplot2::aes_string("x", "y_inv",
                                          color = "status_2_territories"
                                          # alpha = point_alpha
                      )) +
    ggplot2::xlim(14, 33) +
    ggplot2::ylim(14, 35) +
    # Layer 1 - Plot image
    ggplot2::annotation_custom(st_se@images[[1]],
                               xmin = 14,
                               xmax = 33,
                               ymin = 14,
                               ymax = 35) +
    # Layer 2 - Plot points
    geom_point(size = 5, alpha = 0.8) +
    labs(color = "Tissue stratification") +
    coord_fixed(1) +
    theme_void() +
    scale_color_manual(values = c("#00A087FF", "#E64B35FF"))

  tmp_plt_2
  save_plot("R/analyses/img/PDAC-B/PDAC_B tissue stratification.jpg",tmp_plt_2,
            bg="white", base_width = 10,
            base_height = 9)
}


# Joining immune cell type information
decon_mtrx_immune = readRDS("R/analyses/data/PDAC/decon_mtrx_ica_pdac_b.RDS")
decon_mtrx_immune[decon_mtrx_immune > 0] <- 1

decon_df_immune <- data.frame(decon_mtrx_immune)
decon_df_immune$status_2_territories <- st_se@meta.data[, "status_2_territories"]

{
  tmp_df <- decon_df_immune %>%
    dplyr::select(-c("res_ss")) %>%
    dplyr::mutate(total = 1) %>%
    dplyr::group_by(status_2_territories) %>%
    dplyr::summarise_if(is.numeric, ~sum(.)) %>%
    data.frame() %>%
    tibble::column_to_rownames("status_2_territories")

  prop_df <- tmp_df[, 1:(ncol(tmp_df) - 1)] / tmp_df[, ncol(tmp_df)]
  tmp_dif <- sapply(colnames(prop_df), function(i) prop_df[1, i] - prop_df[2, i])

  group_prop_differences <- function(df, grp_vr) {
    tmp_df <- df %>%
      # dplyr::select(-c("status_3_territories", "res_ss")) %>%
      dplyr::mutate(total = 1) %>%
      dplyr::group_by_at(grp_vr) %>%
      dplyr::summarise_if(is.numeric, ~sum(.)) %>%
      data.frame() %>%
      tibble::column_to_rownames(grp_vr)

    prop_df <- tmp_df[, 1:(ncol(tmp_df) - 1)] / tmp_df[, ncol(tmp_df)]
    tmp_dif <- sapply(colnames(prop_df), function(i) prop_df[1, i] - prop_df[2, i])

    return(tmp_dif)
  }
}
# Permutation test
{
  set.seed(123)

  decon_df_immune2 <- decon_df_immune %>%
    dplyr::select(-"res_ss")

  tmp <- as.data.frame(lapply(decon_df_immune2, sample)) %>%
    mutate(status_2_territories = decon_df_immune2$status_2_territories)
  # Get difference distributions
  perm_diff <- sapply(1:10000, function(iter){
    # browser()
    shuff_df <- as.data.frame(lapply(decon_df_immune2, sample)) %>%
      mutate(status_2_territories = decon_df_immune2$status_2_territories)

    tmp <- group_prop_differences(df = shuff_df, grp_vr = "status_2_territories")
    return(tmp)
  })

  # Get double tailed pvalue
  perm_pvals <- sapply(seq_len(length(tmp_dif)), function(i){
    # Get where on the distribution the observed value falls
    tmp_pval <- mean(perm_diff[i, ] >= tmp_dif[i])

    # Get the double tailed p value
    pval <- min(c(tmp_pval, 1 - tmp_pval))
    return(pval)
  })

  adj_perm_pvals <- p.adjust(p = perm_pvals, method = "bonferroni")

  perm_df <- data.frame(t(adj_perm_pvals))
  names(perm_df) <- names(prop_df)

  hm_df_2 <- dplyr::bind_rows(prop_df, perm_df)
  rownames(hm_df_2) <- c(rownames(prop_df), "pval")
}
# Plot heatmap
{
  cell_types_plt <- sort(unique(ica_se$cluster))
  pdac_plt_names <- data.frame(df_name = gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                              x = cell_types_plt,
                                              perl = TRUE),
                               plt_name = cell_types_plt,
                               col_ct = c("#1f9747","#07b1c0","#0371ae","#f27d2f","#d91f2d","#8066a8","#814844","#787879","#d375ae","#a8b32e"))

  hm_tmp <- t(hm_df_2) %>%
    data.frame() %>%
    tibble::rownames_to_column("cell_type") %>%
    tidyr::pivot_longer(cols = c("Non.Tumoral", "Tumoral", "pval"),
                        names_to = "status",
                        values_to = "proportion") %>%
    dplyr::mutate(
      status = factor(status, levels = c("Non.Tumoral", "Tumoral", "pval"), labels = c("Non", "Tumoral", "Pval"))
    ) %>%
    dplyr::left_join(pdac_plt_names, by = c("cell_type" = "df_name")) %>%
    dplyr::mutate(proportion_mod = dplyr::if_else(status == "Pval", 0, proportion),
                  plt_name = as.character(plt_name),
                  plt_name = dplyr::if_else(plt_name == "Macro_1", "Macrophage", plt_name),
                  plt_name = dplyr::if_else(plt_name == "Mono", "Monocyte", plt_name),
                  plt_name = stringr::str_wrap(string = plt_name, width = 20),
                  pval = dplyr::if_else(status == "Pval", proportion, NA_real_)) %>%
    ggplot(ggplot2::aes(x = status,
                        y = plt_name,
                        fill= proportion_mod)) +
    ggplot2::geom_tile(ggplot2::aes(alpha = proportion_mod)) +
    ggplot2::geom_text(ggplot2::aes(label = round(pval, 2)),
                       show_guide  = FALSE,
                       size = 7,
                       family = "Arial") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Cell type enrichment per region",
                  y = "",
                  x = "",
                  fill = "Proportion",
                  alpha = "Proportion") +
    ggplot2::scale_fill_gradient(low = "#FFF5F0",
                                 high = "#99000D",
                                 guide = "legend",
                                 limits = c(0, 1),
                                 breaks = seq(0, 1, 0.1)) +
    ggplot2::scale_alpha_continuous(guide = "legend",
                                    limits = c(0, 1),
                                    breaks = seq(0, 1, 0.1)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 20),
      axis.text.y = ggplot2::element_text(size = 20),
      legend.text = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_text(size = 22),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title = ggplot2::element_text(size = 22),
      axis.line = ggplot2::element_blank())

  hm_tmp
  save_plot("R/analyses/img/PDAC-B/Cell type emrichment per region.jpg",hm_tmp,
            bg="white", base_width = 8,
            base_height = 10)


  hm_tmp_uv <- t(hm_df_2[1:2,]) %>%
    data.frame() %>%
    tibble::rownames_to_column("cell_type") %>%
    tidyr::pivot_longer(cols = c("Non.Tumoral", "Tumoral"),
                        names_to = "status",
                        values_to = "proportion") %>%
    dplyr::mutate(
      status = factor(status, levels = c("Non.Tumoral", "Tumoral"), labels = c("Non", "Tumoral"))
    ) %>%
    dplyr::left_join(pdac_plt_names, by = c("cell_type" = "df_name")) %>%
    dplyr::mutate(proportion_mod = dplyr::if_else(status == "Adj. Pval", 0, proportion),
                  plt_name = as.character(plt_name),
                  plt_name = dplyr::if_else(plt_name == "Macro_1", "Macrophage", plt_name),
                  plt_name = dplyr::if_else(plt_name == "Mono", "Monocyte", plt_name),
                  plt_name = stringr::str_wrap(string = plt_name, width = 20),
    ) %>%
    ggplot(ggplot2::aes(x = status,
                        y = plt_name,
                        fill= proportion_mod)) +
    ggplot2::geom_tile(ggplot2::aes(alpha = proportion_mod)) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Cell type enrichment per region",
                  y = "",
                  x = "",
                  fill = "Proportion",
                  alpha = "Proportion") +
    ggplot2::scale_fill_gradient(low = "#FFF5F0",
                                 high = "#99000D",
                                 guide = "legend",
                                 limits = c(0, 1),
                                 breaks = seq(0, 1, 0.1)) +
    ggplot2::scale_alpha_continuous(guide = "legend",
                                    limits = c(0, 1),
                                    breaks = seq(0, 1, 0.1)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 20),
      axis.text.y = ggplot2::element_text(size = 20),
      legend.text = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_text(size = 22),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title = ggplot2::element_text(size = 22),
      axis.line = ggplot2::element_blank())
  hm_tmp_uv
  save_plot("R/analyses/img/PDAC-B/Cell type emrichment per region unpv.jpg",hm_tmp_uv,
            bg="white", base_width = 8.2,
            base_height = 6)
}

{
  decon_df_immune = data.frame(readRDS("R/analyses/data/PDAC/decon_mtrx_ica_pdac_b.RDS"))
  decon_df_immune$status_2_territories <- st_se@meta.data[, "status_2_territories"]
  decon_mtrx_long <- decon_df_immune %>%
    dplyr::select(-c("res_ss")) %>%
    dplyr::mutate(total = 1) %>%
    tidyr::pivot_longer(cols = -c(status_2_territories, total),
                        names_to = "cell_type",
                        values_to = "proportion") %>%
    dplyr::left_join(pdac_plt_names, by = c("cell_type" = "df_name")) %>%
    dplyr::mutate(plt_name = as.character(plt_name),
                  plt_name = dplyr::if_else(plt_name == "Macro_1", "Macrophage", plt_name),
                  plt_name = dplyr::if_else(plt_name == "Mono", "Monocyte", plt_name),
                  plt_name = dplyr::if_else(plt_name == "CD8 tumor-reactive (exhausted)", "CD8 tumor-reactive", plt_name))
  # Get Bonferroni adjusted P-values
  y_pos <- decon_mtrx_long %>%
    group_by(cell_type) %>%
    summarise(y_max = max(proportion)) %>%
    dplyr::pull(y_max)

  # annotation table with adjusted pvals and y-position of the labels
  anno_df <- ggpubr::compare_means(proportion ~ status_2_territories,
                                   method = "wilcox.test",
                                   group.by = "plt_name",
                                   data = decon_mtrx_long,
                                   p.adjust.method = "bonferroni") %>%
    mutate(y_pos = y_pos - 0.1 * y_pos,
           p.adj.txt = paste("p =", p.adj))
}
# Plot faceted boxplots All cell types
{
  ct_interest <- unique(ica_se@meta.data[["cluster"]])
  bplt <- decon_mtrx_long %>%
    left_join(anno_df, by = "plt_name") %>%
    filter(plt_name %in% ct_interest) %>%
    filter(p.adj < 0.05) %>%
    ggpubr::ggboxplot(.,
                      x = "status_2_territories",
                      y = "proportion",
                      color = "status_2_territories",
                      fill = "status_2_territories",
                      add = "jitter",
                      alpha = 0.6,
                      facet.by = "plt_name",
                      scales = "free",
                      # palette = "npg",
                      palette = c("#00A087FF", "#E64B35FF"),
                      c = 2,
                      outlier.shape = NA) +
    # ggplot2::facet_wrap(~plt_name) +
    ggpubr::geom_signif(
      data = anno_df %>% filter(p.adj < 0.05 & plt_name %in% ct_interest),
      aes(xmin = group1,
          xmax = group2,
          annotations = p.adj.txt,
          y_position = y_pos * 1,
          size = 20,
          family = "Arial"),
      manual = TRUE)
  bplt
  bplt <- bplt +
    ggplot2::labs(title = "Spot composition comparison",
                  x = "Tissue stratification",
                  y = "Capture location proportion",
                  fill = "Stratification",
                  color = "Stratification") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial"),
      strip.text = ggplot2::element_text(size = 23, face = "bold"),
      axis.text =  ggplot2::element_text(size = 20),
      axis.title = ggplot2::element_text(size = 25, face = "bold"),
      strip.background.x = ggplot2::element_rect(colour = "transparent",
                                                 fill = "transparent"),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 30),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 15),
      legend.key.size = ggplot2::unit(3, "line"))
  bplt
}

