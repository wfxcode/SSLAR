library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(Matrix)
library(Seurat)

library(wesanderson) # PDAC-A颜色盘
library(RColorBrewer)# PDAC-B颜色盘
library(scales)
source("R/analyses/utils/bin.r")



# PDAC-B
pdac_b_counts <- readr::read_tsv("R/analyses/data/PDAC/PDAC-B-indrop-filtered.txt", col_names = FALSE)
colnames(pdac_b_counts) <- c("Genes", paste("cell_a", 2:ncol(pdac_b_counts), sep = "_"))

# Find duplicate genes
pdac_b_counts$Genes[duplicated(pdac_b_counts$Genes)]

# Remove duplicates based on Genes columns
pdac_b_counts <- pdac_b_counts[!duplicated(pdac_b_counts$Genes), ] %>%
  tibble::column_to_rownames("Genes")

pdac_b_counts_mtrx <- Matrix::Matrix(data.matrix(pdac_b_counts[-1, ]), sparse = TRUE)

annotation_b <- data.frame(t(pdac_b_counts[1, ]))
colnames(annotation_b) <- "annotation"
pdac_B <- CreateSeuratObject(counts = pdac_b_counts_mtrx,
                             project = "pdac_b",
                             assay = "RNA",
                             meta.data = annotation_b)

pdac_B <- Seurat::SCTransform(object = pdac_B,min_cells=3,method="glmGamPoi")
pdac_B <- Seurat::RunPCA(pdac_B, verbose = FALSE)
Seurat::ElbowPlot(pdac_B, ndims = 50)
pdac_B <- Seurat::FindNeighbors(pdac_B,
                                dims = 1:40)
pdac_B <- Seurat::FindClusters(pdac_B,
                               verbose = FALSE,
                               resolution = c(1, 2, 3, 4, 5))
pdac_B <- Seurat::RunUMAP(pdac_B,
                          dims = 1:40)

saveRDS(object = pdac_B, file = "R/analyses/data/PDAC/PDAC-B_itai_joint.RDS")
pdac_B <- readRDS(file = "R/analyses/data/PDAC/PDAC-B_itai_joint.RDS")

pdac_b_annot <- data.frame(
  annotation = sort(unique(as.character(pdac_B$annotation))),
  plt_name = c("Acinar cells", "Cancer clone (TM4SF1)", "Centroacinar ductal cells",
               "Antigen-presenting ductal cells",
               "Terminal ductal cells", "Endocrine cells",
               "Endothelial cells", "Macrophages", "Mast cells", "mDCs", "Monocytes", "RBCs",
               "Tuft cells"),
  col_ct = c(brewer.pal(12, "Paired"),brewer.pal(1, "Dark2")[1]))

metadata_b_pbmc = data.frame(pdac_B@meta.data) %>%
  left_join(pdac_b_annot, by = "annotation")

pdac_B@meta.data[["annotation"]] <- metadata_b_pbmc[["plt_name"]]
pdac_B@meta.data[["col_ct"]] <- metadata_b_pbmc[["col_ct"]]
Idents(pdac_B) <- pdac_B$annotation
table(pdac_B$annotation)


saveRDS(object = pdac_B, file = "R/analyses/data/PDAC/PDAC-B_itai_processed.RDS")
pdac_B = readRDS("R/analyses/data/PDAC/PDAC-B_itai_processed.RDS")


Idents(pdac_B) <- pdac_B$annotation
dfb_col <- pdac_B@meta.data%>%
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
                ct_names = annotation_mod)
umap_pdac_b <- DimPlot(pdac_B, reduction = "umap", group.by = "annotation") +
  scale_color_manual(values = as.character(dfb_col[as.character(dfb_col$plt_name) %in% pdac_B$annotation, "col_ct"]))

save_plot("R/analyses/img/PDAC/Figure_UMAP_PDAC_B.jpg",umap_pdac_b,
          base_width = 10,
          base_height = 6)
# ggpubr::ggexport(plotlist = list(umap_pdac_b),
#                  filename = sprintf("%s/%s/Supplementary_Figure_QQQ_UMAP_PDAC-B.pdf",
#                                     an_pdac,plt_dir),
#                  width = 12,
#                  height = 9,
#                  res = 600)
