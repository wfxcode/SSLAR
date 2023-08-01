source("R/method/import.R")
se_sc = readRDS("data/PBMC_Smart-Seq2.rds")
se_sc = Seurat::SCTransform(se_sc,verbose = FALSE)
se_sc@assays[["RNA"]] = se_sc@assays[["SCT"]]
Seurat::Idents(se_sc) = se_sc$subclass

marker_genes = Seurat::FindAllMarkers(object = se_sc,
                                       assay = "RNA",
                                       min.pct = 0,
                                       only.pos = TRUE,
                                       logfc.threshold = 0)

performance = test_synthetic_mixtures(se_sc = se_sc,
                                       clust_vr = "subclass",
                                       cluster_markers = marker_genes)
