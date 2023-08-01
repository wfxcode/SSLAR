{
  source('R/import.R')
  library(cli)
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(scrattch.io)
  library(RColorBrewer)
  library(Spaniel)
  library(ggplot2)
  library(cowplot)
}

clust_vr <- "annotation"
cl_n = 99999; hvg = 99999
indrop_pdac_a <- readRDS(file = "R/analyses/data/PDAC/PDAC-A_itai_processed.RDS")
# Remove RBCs
# indrop_pdac_a <- indrop_pdac_a[, indrop_pdac_a$annotation != "RBCs"]
# ST data
{
  st_tibble = read.csv(file = "R/analyses/data/PDAC/PDAC-A-ST1-filtered.csv",header = TRUE,row.names=1)
  # st_position = read.csv(file = "R/analyses/data/PDAC/PDAC-A-ST1_positions.csv",header = TRUE,row.names=1)
  barcodes <- utils::read.csv("R/analyses/data/PDAC/PDAC-A-ST1_positions.csv", sep = "\t", header = FALSE)

  st_se = Spaniel::createSeurat(counts = st_tibble,
                                barcodeFile = ("R/analyses/data/PDAC/PDAC-A-ST1_positions.csv"),
                                projectName = "PDAC-A",
                                sectionNumber = "1")
  st_se <- Seurat::ScaleData(st_se)
  st_se = SCTransform(st_se)
  image <- Spaniel::parseImage("R/analyses/data/PDAC/PDAC-A-ST1.jpg")
  st_se@images <- list(image)
  saveRDS(object = st_se,
          file = sprintf("R/analyses/data/PDAC/PDAC-A_ST.RDS"))
  st_se = readRDS("R/analyses/data/PDAC/PDAC-A_ST.RDS")
}
## Spatial Deconvolution
Seurat::Idents(object = indrop_pdac_a) <- indrop_pdac_a@meta.data[, clust_vr]
cluster_markers_a <- Seurat::FindAllMarkers(object = indrop_pdac_a,
                                            verbose = TRUE,
                                            only.pos = TRUE,
                                            assay = "SCT",
                                            slot = "data")
# 只保留avg_log2FC>1且pct.1 > 0.75
cluster_markers_a <- cluster_markers_a %>%
  filter(avg_log2FC > 1 & pct.1 > 0.75)
saveRDS(object = cluster_markers_a,
        file = "R/analyses/data/PDAC/PDAC_A_markers.rds")
cluster_markers_a <- readRDS("R/analyses/data/PDAC/PDAC_A_markers.rds")

# cluster_markers_filt_a$cluster <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
#                                        x = cluster_markers_filt_a$cluster,
#                                        perl = TRUE)
indrop_pdac_a_down <- downsample_se_obj(se_obj = indrop_pdac_a, clust_vr = clust_vr,
                                cluster_markers = cluster_markers_a, cl_n = cl_n, hvg = 11000)

nmf_mod_ls <- train_nmf1(cluster_markers = cluster_markers_a,
                         se_sc = indrop_pdac_a_down,
                         mtrx_spatial = st_se@assays$SCT@counts,
                         ntop = NULL,
                         transf = "uv",
                         clust_vr = clust_vr,
                         nmf = "ssNMF")

nmf_mod <- nmf_mod_ls[[1]]
### Spot Deconvolution
# Extract count matrix
spot_counts <- st_se@assays$SCT@counts
# Subset to genes used to train the model
spot_counts <- spot_counts[rownames(spot_counts) %in% rownames(basis(nmf_mod)), ]

ct_topic_profiles <- topic_profile_per_cluster_nmf(h = coef(nmf_mod),
                                                   train_cell_clust = nmf_mod_ls[[2]])
decon_mtrx <- mixture_deconvolution_nmf1(nmf_mod = nmf_mod,
                                         transf =  "uv",
                                         mixture_transcriptome = spot_counts,
                                         reference_profiles = ct_topic_profiles,
                                         min_cont = 0.01)
saveRDS(decon_mtrx,"R/analyses/data/PDAC/decon_mtrx_pdac_a.RDS")
decon_mtrx = readRDS("R/analyses/data/PDAC/decon_mtrx_pdac_a.RDS")

# Plot the scatterplot

{
  decon_mtrx_subs <- decon_mtrx[, colnames(decon_mtrx)[! colnames(decon_mtrx) %in% "res_ss"]]
  colnames(decon_mtrx_subs) <- gsub(pattern = "[[:punct:]]|[[:blank:]]",
                                    replacement = ".",
                                    x = colnames(decon_mtrx_subs),
                                    perl = TRUE)
  st_se@meta.data <- cbind(st_se@meta.data, decon_mtrx_subs)
  ## Preprocess data
  spatial_coord <- data.frame(st_se@meta.data) %>%
    tibble::rownames_to_column("ID")

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

  cell_types <- dfa_col$df_names


  plt_df <- data.frame(plt_name = cell_types,
                       df_name = gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                      x = cell_types,
                                      perl = TRUE),
                       col_ct = dfa_col$col_ct)

  tmp_df <- data.frame(df_name = colnames(spatial_coord)[colnames(spatial_coord) %in% plt_df$df_name]) %>%
    left_join(plt_df)
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
                                pie_scale = 0.9) +
    scale_y_reverse() +
    theme_half_open(11, rel_small = 1) +
    theme_void() +
    coord_fixed(ratio = 1) +
    scale_fill_manual(values = tmp_df[tmp_df$plt_name %in% ct_all, "col_ct"]) +
    labs(title = "PDAC-A Spatial scatterpie",
         color = "Cell types") +
    theme(
      legend.position = 'none',
      # plot.background = element_rect(fill = "#FFFFFF"),
      # panel.background = element_blank(),
      # plot.margin = margin(20, 20, 20, 20),
      plot.title = element_text(hjust = 0.5, size = 20))

  save_plot("R/analyses/img/PDAC/PDAC_all.jpg",scatterpie_plt,
            bg="white", base_width = 10,
            base_height = 9)

}
{
  st_se = SCTransform(st_se)
  st_se@assays$RNA@counts = st_se@assays$SCT@counts

  spots = colnames(st_se)
  rownames(decon_mtrx) = colnames(st_se)
  tmp = as.matrix(st_se@assays$SCT@counts)
  rownames(decon_mtrx) = spots
  count = 1
  for (spot in spots) {
    s = sum(decon_mtrx[spot,])
    # print(sum(st_se@assays$SCT@counts[,spot]))
    # print(sum(decon_mtrx_subs[spot,]))
    if(s == 0){
      count = count+1
    }
    # print(s == 0)
    # print(sprintf("%d+%f",sum(st_se@assays$RNA@counts[,spot]),sum(decon_mtrx_subs[spot,])))
  }

  tmp = read.csv("R/analyses/data/PDAC/tmp.csv")
  count = 1
  for (spot in spots) {
    s = sum(tmp[,spot])
    # print(sum(st_se@assays$RNA@counts[,spot]))
    # print(sum(decon_mtrx_subs[spot,]))
    if(s == 0){
      count = count+1
    }
  }
}
