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
  library(cowplot)
  # library(flextable)
  source("R/analyses/utils/bin.r")
  source("R/analyses/utils/spatial_plot_spaniel.R")
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
}
{
  pdac_A = readRDS("R/analyses/data/PDAC/PDAC-A_itai_processed.RDS")
  cell_types_plt <- sort(unique(pdac_A$annotation))
  decon_mtrx_pdac_a = readRDS("R/analyses/data/PDAC/decon_mtrx_pdac_a.RDS")

  pdac_plt_names <- pdac_A@meta.data%>%
    dplyr::select(annotation,
                  col_ct) %>%
    dplyr::distinct() %>%
    dplyr::arrange(annotation) %>%
    dplyr::mutate(annotation_mod = stringr::str_replace_all(
      string = annotation,
      pattern = "[[:punct:]]|[[:space:]]",
      replacement = ".")) %>%
    dplyr::rename(col_ct = col_ct,
                  plt_name = annotation,
                  df_name = annotation_mod)

  decon_mtrx <- decon_mtrx_pdac_a[, colnames(decon_mtrx_pdac_a)[!colnames(decon_mtrx_pdac_a) %in% "res_ss"]]

  colnames(decon_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                               x = colnames(decon_mtrx),
                               perl = TRUE)
  new_names <- data.frame(df_name = colnames(decon_mtrx)) %>%
    left_join(pdac_plt_names, by = "df_name") %>%
    pull(plt_name)

  colnames(decon_mtrx) <- new_names
  graph_ntw <- get_spatial_interaction_graph(decon_mtrx = decon_mtrx)
  # tmp = colSums(decon_mtrx > 0) / 10
  # deg <- (((tmp - min(tmp)) / (max(tmp) - min(tmp)))+1)*6
  deg <- colSums(decon_mtrx > 0) / 15

  # Get color palette for difusion
  # edge_importance <- E(graph_ntw)$importance
  tmp = E(graph_ntw)$importance
  edge_importance = ((tmp - min(tmp)) / (max(tmp) - min(tmp)))*2+0.4
  # Select a continuous palette
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
  # Create a color vetor
  getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[4:9])
  # Get how many values we need
  grad_edge <- seq(0, max(edge_importance), 0.1)
  # Generate extended gradient palette dataframe
  graph_col_df <- data.frame(value = as.character(grad_edge),
                             color = getPalette(length(grad_edge)),
                             stringsAsFactors = FALSE)
  # Assign color to each edge
  color_edge <- data.frame(value = as.character(round(edge_importance, 1)), stringsAsFactors = FALSE) %>%
    dplyr::left_join(graph_col_df, by = "value") %>%
    dplyr::pull(color)

  plot(graph_ntw,
       # Size of the edge
       edge.width = edge_importance,
       edge.color = color_edge,
       # Size of the buble
       vertex.size = deg,
       vertex.color = "#cde394",
       vertex.frame.color = "white",
       vertex.label.color = "black",
       vertex.label.family = "Helvetica", # Font family of the label (e.g.“Times”, “Helvetica”)
       layout = layout.circle,
       main = "PDAC-A spatial interaction network")
}

ica_se = readRDS("R/analyses/data/PDAC/ica_se_processed.RDS")
cell_types_plt <- sort(unique(ica_se$cluster))
decon_mtrx_immune = readRDS("R/analyses/data/PDAC/decon_mtrx_ica_pdac_a.RDS")
pdac_plt_names <- data.frame(df_name = gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                            x = cell_types_plt,
                                            perl = TRUE),
                             plt_name = cell_types_plt,
                             col_ct = c("#1f9747","#07b1c0","#0371ae","#f27d2f","#d91f2d","#8066a8","#814844","#787879","#d375ae","#a8b32e"))
decon_mtrx <- decon_mtrx_immune[, colnames(decon_mtrx_immune)[!colnames(decon_mtrx_immune) %in% "res_ss"]]
{
  colnames(decon_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                               x = colnames(decon_mtrx),
                               perl = TRUE)
  new_names <- data.frame(df_name = colnames(decon_mtrx)) %>%
    left_join(pdac_plt_names, by = "df_name") %>%
    pull(plt_name)

  colnames(decon_mtrx) <- new_names
  graph_ntw <- get_spatial_interaction_graph(decon_mtrx = decon_mtrx)
  tmp = colSums(decon_mtrx > 0) / 10
  deg <- (((tmp - min(tmp)) / (max(tmp) - min(tmp)))+1)*6
  # deg <- colSums(decon_mtrx > 0) / 10

  # Get color palette for difusion
  # edge_importance <- E(graph_ntw)$importance
  tmp = E(graph_ntw)$importance
  edge_importance = ((tmp - min(tmp)) / (max(tmp) - min(tmp)))*2+0.4
  # Select a continuous palette
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
  # Create a color vetor
  getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd")[5:9])
  # Get how many values we need
  grad_edge <- seq(0, max(edge_importance), 0.1)
  # Generate extended gradient palette dataframe
  graph_col_df <- data.frame(value = as.character(grad_edge),
                             color = getPalette(length(grad_edge)),
                             stringsAsFactors = FALSE)
  # Assign color to each edge
  color_edge <- data.frame(value = as.character(round(edge_importance, 1)), stringsAsFactors = FALSE) %>%
    dplyr::left_join(graph_col_df, by = "value") %>%
    dplyr::pull(color)

  plot(graph_ntw,
       # Size of the edge
       edge.width = edge_importance,
       edge.color = color_edge,
       # Size of the buble
       vertex.size = deg,
       vertex.color = "#cde394",
       vertex.frame.color = "white",
       vertex.label.color = "black",
       vertex.label.family = "Helvetica", # Font family of the label (e.g.“Times”, “Helvetica”)
       layout = layout.circle,
       main = "PDAC-A immune spatial interaction network")
}
