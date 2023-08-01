## Libraries
source('R/import.R')
library(cli)
library(Seurat)
library(ggplot2)
library(dplyr)
library(purrr)
library(SPOTlight)
library(NMF)
library(nnls)
library(scrattch.io)
library(xlsx)

## Paths
tech <- "sc"
tissue <- "allen_ref_70k"
dwn_smplng <- "both"
org <- "mm"
source("R/analyses/misc/paths_vrs.R")
seed_id <- 123
set.seed(seed_id)

## Set common parameters
clust_vr <- "subclass_label"
cl_n <- 999999
method <- "ssNMF"
transf <- "uv"
hvg <- 999999
FC <- 1
pct1 <- 0.9
date <- Sys.Date()

id_nmf <- sprintf("cln-%s_transf-%s_method-%s_hvg-%s_FC-%s_pct1-%s_seed-%s",
                  cl_n, transf, method, hvg, FC, pct1, seed_id)
options(stringsAsFactors = FALSE)

library(RColorBrewer)
n <- 60
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

## Load data
### scRNAseq-
#### Load data 70k and convert to Seurat object
allen_ref_70k = readRDS("R/analyses/data/mouse_brain/allen_ref_70k")

# Scale and normalize the data
allen_ref_70k <- Seurat::SCTransform(object = allen_ref_70k)
#
saveRDS(object = allen_ref_70k, file = "R/analyses/data/mouse_brain/allen_ref_70k_processed.RDS")
allen_ref_70k <- readRDS(file = "R/analyses/data/mouse_brain/allen_ref_70k_processed.RDS")



# Get cluster markers so we can downsample since it is a very large dataset.
Seurat::Idents(object = allen_ref_70k) <- allen_ref_70k@meta.data[, clust_vr]
cluster_markers_all <- Seurat::FindAllMarkers(object = allen_ref_70k,
                                              verbose = TRUE,
                                              only.pos = TRUE,
                                              assay = "SCT",
                                              slot = "data",
                                              min.pct = 0.9,
                                              max.cells.per.ident = 100)
#
saveRDS(object = cluster_markers_all, file = "R/analyses/data/mouse_brain/markers_allen_ref_70k.RDS")
cluster_markers_all <- readRDS(file = "R/analyses/data/mouse_brain/markers_allen_ref_70k.RDS")

### Downsampling + Data preprocessing
cluster_markers_filt <- cluster_markers_all %>%
  filter(avg_log2FC > FC & pct.1 > pct1)

cluster_markers_filt %>%
  count(cluster) %>%
  data.frame()

cluster_markers_filt %>% filter(cluster == "Astro")
table(unique(cluster_markers_all$cluster) %in% unique(cluster_markers_filt$cluster))

# As we can see all the cell types still have markers after filtering.
dim(allen_ref_70k)
allen_ref_70k_down <- downsample_se_obj(se_obj = allen_ref_70k,
                                        clust_vr = clust_vr,
                                        cluster_markers =cluster_markers_filt,
                                        cl_n = cl_n,
                                        hvg = hvg)
dim(allen_ref_70k_down)
saveRDS(object = allen_ref_70k_down,
        file = "R/analyses/data/mouse_brain/allen_ref_70k_downsampled_new.RDS")
allen_ref_70k_down <- readRDS(file = "R/analyses/data/mouse_brain/allen_ref_70k_downsampled_new.RDS")


### Spatial data
#### Anterior slice
anterior <- SeuratData::LoadData("stxBrain", type = "anterior1")
anterior$slice <- "Anterior"
anterior <- Seurat::SCTransform(anterior, assay = "Spatial", verbose = TRUE)

#### Posterior slice
posterior <- SeuratData::LoadData("stxBrain", type = "posterior1")
posterior$slice <- "Posterior"
posterior <- Seurat::SCTransform(posterior, assay = "Spatial", verbose = TRUE)

#### Merging slices
brain <- merge(anterior, posterior)

Seurat::DefaultAssay(brain) <- "SCT"
Seurat::VariableFeatures(brain) <- c(Seurat::VariableFeatures(anterior),
                                     Seurat::VariableFeatures(posterior))
brain <- Seurat::RunPCA(brain, verbose = FALSE)
Seurat::ElbowPlot(object = brain, ndims = 50)
brain <- Seurat::FindNeighbors(brain, dims = 1:40)
brain <- Seurat::FindClusters(brain, verbose = FALSE)
brain <- Seurat::RunUMAP(brain, dims = 1:40)

saveRDS(object = brain,
        file = "R/analyses/data/mouse_brain/brain1_processed.RDS")
brain <- readRDS("R/analyses/data/mouse_brain/brain1_processed.RDS")


# Visualization of the merged slices
Seurat::DimPlot(brain,
                reduction = "umap",
                group.by = c("ident", "slice"))

Seurat::SpatialDimPlot(brain)
Seurat::SpatialFeaturePlot(object = brain, features = "Cplx3")

## Spatial Deconvolution
### Train NMF model
# for (i in seq(2,10)) {

start_time <- Sys.time()
nmf_mod_ls <- train_nmf1(cluster_markers = cluster_markers_filt,
                         se_sc = allen_ref_70k_down,
                         mtrx_spatial = brain@assays$Spatial@counts,
                         ntop = NULL,
                         transf = transf,
                         clust_vr = clust_vr,
                         nmf = method)

tot_time <- difftime(time1 = Sys.time(), time2 =start_time , units = "mins")
tot_time

saveRDS(nmf_mod_ls,"R/analyses/data/mouse_brain/nmf_mod_ls_new.RDS")

nmf_mod_ls <- readRDS("R/analyses/data/mouse_brain/nmf_mod_ls_new.RDS")
nmf_mod <- nmf_mod_ls[[1]]

# Extract matrices form the model:
# get matrix W
w <- basis(nmf_mod)
dim(w)

# get matrix H
h <- coef(nmf_mod)
dim(h)

rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- dot_plot_profiles_fun(h = h,
                                            train_cell_clust = nmf_mod_ls[[2]])
topic_profile_plts[[2]] + theme(axis.text.x = element_text(angle = 90))

### Spot Deconvolution
# Extract count matrix
spot_counts <- brain@assays$Spatial@counts

# Subset to genes used to train the model
spot_counts <- spot_counts[rownames(spot_counts) %in% rownames(w), ]

#### Get mixture composition
# Run test spots through the basis to get the pertinent coefficients.
# To do this for every spot we are going to set up a system of linear equations
# where we need to find the coefficient, we will use non-negative least squares to determine the best coefficient fit.
ct_topic_profiles <- topic_profile_per_cluster_nmf(h = h,
                                                   train_cell_clust = nmf_mod_ls[[2]])

  # write.xlsx(ct_topic_profiles, "R/analyses/data/mouse_brain/celltype_topic/ct_topic.xls",
  #            sheetName = sprintf("Sheet %d",i), append = TRUE)
# }

decon_mtrx <- mixture_deconvolution_nmf1(nmf_mod = nmf_mod,
                                         transf =  "uv",
                                         mixture_transcriptome = spot_counts,
                                         reference_profiles = ct_topic_profiles,
                                         min_cont = 0.1)
saveRDS(decon_mtrx,"R/analyses/data/mouse_brain/decon_mtrx.RDS")
decon_mtrx <- readRDS("R/analyses/data/mouse_brain/decon_mtrx.RDS")
