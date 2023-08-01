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


## Load data
### scRNAseq filtered matrices
#### PDAC_A
pdac_a_counts <- readr::read_tsv("R/analyses/data/PDAC/PDAC-A-indrop-filtered.txt", col_names = FALSE)
colnames(pdac_a_counts) <- c("Genes", paste("cell_a", 2:ncol(pdac_a_counts), sep = "_"))

# Find duplicate genes
pdac_a_counts$Genes[duplicated(pdac_a_counts$Genes)]

# Remove duplicates based on Genes columns
pdac_a_counts <- pdac_a_counts[!duplicated(pdac_a_counts$Genes), ] %>%
  tibble::column_to_rownames("Genes")

pdac_a_counts_mtrx <- Matrix::Matrix(data.matrix(pdac_a_counts[-1, ]), sparse = TRUE)
annotation_a <- data.frame(t(pdac_a_counts[1, ]))
colnames(annotation_a) <- "annotation"

# Create Seurat object
pdac_A <- CreateSeuratObject(counts = pdac_a_counts_mtrx,
                             project = "pdac_a",
                             assay = "RNA",
                             meta.data = annotation_a)
## Scale and normalize the data
pdac_A <- Seurat::SCTransform(object = pdac_A,min_cells=3,method="glmGamPoi")
pdac_A <- Seurat::RunPCA(pdac_A, verbose = FALSE)
Seurat::ElbowPlot(pdac_A, ndims = 50)

# From the elbow plot we can see that the elbow is around 12 so we will use the first 40 PC to proceed with the analysis.
pdac_A <- Seurat::FindNeighbors(pdac_A,
                                dims = 1:40)
pdac_A <- Seurat::FindClusters(pdac_A,
                               verbose = FALSE,
                               resolution = c(1, 2, 3, 4, 5))
pdac_A <- Seurat::RunUMAP(pdac_A,
                          dims = 1:40)



pdac_A <- RunTSNE(pdac_A)

saveRDS(object = pdac_A,
        file = "R/analyses/data/PDAC/PDAC-A_itai_joint.RDS")
pdac_A <- readRDS(file = "R/analyses/data/PDAC/PDAC-A_itai_joint.RDS")



# Change names to match the original manuscript
pdac_a_annot <- data.frame(
  annotation = sort(unique(as.character(pdac_A$annotation))),
  plt_name = c("Acinar cells", "Cancer clone (TM4SF1)", "Cancer clone (S100A4)",
               "High/Hypoxic ductal cells (APOL1)", "Centroacinar ductal cells",
               "Antigen-presenting ductal cells",
               "Terminal ductal cells", "Endocrine cells", "Endothelial cells",
               "Fibroblasts", "Macrophages M2", "Macrophages M1", "Mast cells",
               "mDCs", "mDCs", "Monocytes", "pDCs", "RBCs", "T cells & NK cells",
               "Tuft cells"),
  col_ct = c(wes_palette("Royal2")[1:5], wes_palette("Zissou1")[1:5],
             wes_palette("Darjeeling1")[1:4], wes_palette("Darjeeling1")[4:5],
             wes_palette("Darjeeling2")[1:4]))

pdac_a_annot[10,"col_ct"] = pdac_a_annot[12,"col_ct"]; pdac_a_annot[12,"col_ct"] = "#fc694d"

metadata_a_pbmc = data.frame(pdac_A@meta.data) %>%
  left_join(pdac_a_annot, by = "annotation")

pdac_A@meta.data[["annotation"]] <- metadata_a_pbmc[["plt_name"]]
pdac_A@meta.data[["col_ct"]] <- metadata_a_pbmc[["col_ct"]]
Idents(pdac_A) <- pdac_A$annotation
table(pdac_A$annotation)
# Save Seurat object
saveRDS(object = pdac_A,
        file = "R/analyses/data/PDAC/PDAC-A_itai_processed.RDS")
pdac_A = readRDS("R/analyses/data/PDAC/PDAC-A_itai_processed.RDS")

##### UMAP Author's Annotation
dfa_col <- pdac_A@meta.data%>%
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

Idents(pdac_A) <- pdac_A$annotation
umap_pdac_a <- DimPlot(pdac_A, reduction = "umap", group.by = "annotation") +
  scale_color_manual(values = as.character(dfa_col[as.character(dfa_col$plt_name) %in% pdac_A$annotation, "col_ct"]))

save_plot("R/analyses/img/PDAC/Figure_UMAP_PDAC_A.jpg",umap_pdac_a,
          base_width = 10,
          base_height = 6)

tsne_pdac_a <- DimPlot(pdac_A, reduction = "tsne", group.by = "annotation") +
  scale_color_manual(values = as.character(dfa_col[as.character(dfa_col$plt_name) %in% pdac_A$annotation, "col_ct"]))

save_plot("R/analyses/img/PDAC/Figure_TSNE_PDAC_A.jpg",tsne_pdac_a,
          base_width = 10,
          base_height = 6)
# ggpubr::ggexport(plotlist = list(umap_pdac_a),
#                  filename = "R/analyses/img/PDAC/Supplementary_Figure_UMAP_PDAC-A.pdf",
#                  width = 9,
#                  height = 6,
#                  res = 600)


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
