library(Seurat)
library(ggplot2)
library(dplyr)

allen_ref_70k <- readRDS(file = "R/analyses/data/mouse_brain/allen_ref_70k_processed.RDS")
# matedata <- readr::read_csv(file = "R/analyses/data/mouse_brain/metadata.csv")

embed_tsne <- readr::read_csv(file = "R/analyses/data/mouse_brain/2d_coordinates.csv")

tsne_ds <- allen_ref_70k@meta.data %>%
  right_join(embed_tsne, by = "sample_name")

UMAP_allen_70k <- ggplot(tsne_ds, aes(x = Lim1,
                                      y = Lim2,
                                      colour = subclass_label)) +
  geom_point(alpha = 1, size = 0.7) +
  theme_classic() +
  labs(
    title = "Allen Institute 73k cells mouse brain",
    colour = "",
    x = "TSNE 1",
    y = "TSNE 2") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme(
    axis.title = element_text(size = 15),
    plot.title =  element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = col_vector) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 7)))

ggpubr::ggexport(UMAP_allen_70k,
                 filename = sprintf("R/analyses/img/mouse/Supplementary_Figure_allen_%s.pdf",
                                    "73k"),
                 width = 24,
                 height = 18)
