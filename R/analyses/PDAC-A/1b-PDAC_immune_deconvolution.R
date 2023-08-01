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
  library(Seurat)
}
# https://sci.bban.top/pdf/10.1038/s41422-019-0195-y.pdf#view=FitH
immune_count = as.data.frame(read.csv("R/analyses/data/PDAC/immune_count-matrix.txt",sep = " ",header=T,row.names = 1))
immune_count = t(immune_count)
update_annot = read.csv("R/analyses/data/PDAC/immune_all_celltype.txt"
                        ,sep = "\t") %>% dplyr::rename(barcode = cell.name)
unique(update_annot$cluster)
ica_se <- CreateSeuratObject(counts = immune_count,
                             project = "pdac_a",
                             assay = "RNA")

metadata <- ica_se@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::left_join(update_annot,
                   by = c("barcode")) %>%
  dplyr::mutate(specific_cell_type = cluster) %>%
  tibble::column_to_rownames("barcode")

table(update_annot$cluster)
ica_se@meta.data <- metadata
Seurat::Idents(object = ica_se) <- ica_se@meta.data[, "cluster"]



## Scale and normalize the data
ica_se <- Seurat::SCTransform(object = ica_se)
ica_se <- Seurat::RunPCA(ica_se, verbose = FALSE)
Seurat::ElbowPlot(ica_se, ndims = 50)

# From the elbow plot we can see that the elbow is around 12 so we will use the first 40 PC to proceed with the analysis.
ica_se <- Seurat::FindNeighbors(ica_se,
                                dims = 1:50)
ica_se <- Seurat::RunUMAP(ica_se,
                          dims = 1:50)
saveRDS(object = pdac_A,
        file = "R/analyses/data/PDAC/ica_se_joint.RDS")
ica_se <- readRDS(file = "R/analyses/data/PDAC/ica_se_joint.RDS")


ica_se_annot <- data.frame(
  specific_cell_type = sort(unique(as.character(ica_se$specific_cell_type))),
  plt_name = sort(unique(as.character(ica_se$cluster))),
  col_ct = c("#1f9747","#07b1c0","#0371ae","#f27d2f","#d91f2d","#8066a8","#814844","#787879","#d375ae","#a8b32e"))

metadata_ica_se = data.frame(ica_se@meta.data) %>%
  left_join(ica_se_annot, by = "specific_cell_type")
ica_se@meta.data[["col_ct"]] <- metadata_ica_se[["col_ct"]]

Idents(ica_se) <- ica_se$specific_cell_type

saveRDS(object = ica_se,
        file = "R/analyses/data/PDAC/ica_se_processed.RDS")
ica_se = readRDS("R/analyses/data/PDAC/ica_se_processed.RDS")

umap_ica_se <- DimPlot(ica_se, reduction = "umap", group.by = "cluster") +
  scale_color_manual(values = as.character(ica_se_annot[as.character(ica_se_annot$specific_cell_type) %in% ica_se$specific_cell_type, "col_ct"]))

save_plot("R/analyses/img/PDAC/Figure_UMAP_ice_se.jpg",umap_ica_se,
          base_width = 10,
          base_height = 6)

specific_markers <- Seurat::FindAllMarkers(
  object = ica_se,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE)
saveRDS(object = specific_markers,
        file = "R/analyses/data/PDAC/PDAC_specific_markers.rds")
specific_markers <- readRDS("R/analyses/data/PDAC/PDAC_specific_markers.rds")

sp_ct <- ica_se$specific_cell_type
ica_se$specific_cell_type <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                  x = ica_se$specific_cell_type,
                                  perl = TRUE)
st_se = readRDS("R/analyses/data/PDAC/PDAC-A_ST.RDS")

{
  ica_pdac_a_down <- downsample_se_obj(se_obj = ica_se, clust_vr = "cluster",
                                          cluster_markers = specific_markers, cl_n = 99999, hvg = 11000)
  nmf_mod_ls <- train_nmf1(cluster_markers = specific_markers,
                           se_sc = ica_pdac_a_down,
                           mtrx_spatial = st_se@assays$SCT@counts,
                           ntop = NULL,
                           transf = "uv",
                           clust_vr = "cluster",
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
  saveRDS(decon_mtrx,"R/analyses/data/PDAC/decon_mtrx_ica_pdac_a.RDS")
  decon_mtrx =   decon_mtrx = readRDS("R/analyses/data/PDAC/decon_mtrx_ica_pdac_a.RDS")
}
