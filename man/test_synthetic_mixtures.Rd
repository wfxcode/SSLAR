% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_synthetic_mixtures.R
\name{test_synthetic_mixtures}
\alias{test_synthetic_mixtures}
\title{Test Synthetic Mixtures Performance}
\usage{
test_synthetic_mixtures(
  se_sc,
  clust_vr,
  cluster_markers,
  n_syn_mixt = 1000,
  ntop = NULL,
  transf = "uv",
  method = "ssNMF",
  min_cont = 0.07,
  assay = "RNA",
  slot = "counts",
  numiters = 100,
  lam = 100
)
}
\arguments{
\item{se_sc}{Object of class Seurat with the scRNAseq data.}

\item{clust_vr}{Object of class character. Name of the variable containing the cell clustering.}

\item{cluster_markers}{Object of class dataframe obtained from the function Seurat::FindAllMarkers().}

\item{n_syn_mixt}{Object of class integer specifying how many synthetic mixtures to generate.}

\item{ntop}{Object of class "numeric" or NULL; number of unique markers per cluster used to seed the model, by default NULL. If NULL it uses all of them.}

\item{transf}{Transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), sct (Seurat::SCTransform), raw (no transformation applied). By default CPM.}

\item{method}{Object of class character; Type of method to use to find W and H. Look at NMF package for the options and specifications, by default nsNMF.}

\item{min_cont}{Object of class numeric; Indicates the minimum contribution we expect from a cell in that spot. Since we're working with proportions by setting 0.01, by default, means that we will accept those cell types whose weight coefficient is at least 1\\% of the total.}

\item{assay}{Character string specifying the assay to use.}

\item{slot}{Character string specifying the slot to use.}

\item{numiters}{Integer specifying the number of iterations for NMF.}

\item{lam}{Numeric specifying the regularization parameter for NMF.}
}
\value{
This function returns a list where the first element is a list with the NMF model trained and the cell labels, the second is a list with the raw_statistics.
}
\description{
This function takes in a Seurat object with several tuning parameters and assesses its performance on synthetic test spots.
}
\examples{
# Load necessary libraries
library(Seurat)
library(NMF)


# Find all markers
cluster_markers <- FindAllMarkers(se_sc)

# Run the function
results <- test_synthetic_mixtures(
  se_sc = se_sc,
  clust_vr = "seurat_clusters",
  cluster_markers = cluster_markers,
  n_syn_mixt = 1000,
  ntop = 10,
  transf = "uv",
  method = "nsNMF",
  min_cont = 0.07,
  assay = "RNA",
  slot = "counts",
  numiters = 100,
  lam = 100
)

# View the results
print(results$model)
print(results$performance)

}
