# SSLAR

SSLAR conducted cell type annotation by leveraging single cell RNA sequencing and spatial transcriptomics data via machine learning.
![Fig1.jpg](pipeline.jpg)

## Data

All the data used in this article is shown below:

![Fig2.jpg](dataset.jpg)

We provide a dataset for PBMC(smart-seq) as an example.

## Intsall
The version under development is available on GitHub and can be installed as follows:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("wfxcode/SSLAR")
```
This package relies on the ssNMF library, requiring a functional Python interpreter with the library installed. The installPython function offers a convenient method to set up a conda environment named SSLAR-env, which includes all necessary dependencies. While we recommend using this method to install Python libraries, you also have the option to customize the installation process.
```r
library("SSLAR")
SSLAR::installPython()
```
## Run the SSLAR demo in the mouse visual cortex

```r
# Use the specified Conda environment for this R session
# The 'required' parameter ensures that the script will stop if the environment is not available
library(reticulate)
reticulate::use_condaenv("SSLAR-env", required = TRUE)
# Load the required datasets from the SSLAR package
mouse_vc_sc <- SSLAR::mouse_vc_sc
mouse_vc_st <- SSLAR::mouse_vc_st

# Run SSLAR analysis
res <- SSLAR::run.SSLAR(
  sc_count = mouse_vc_sc$sc_count,        # Single-cell count matrix
  sc_label = mouse_vc_sc$sc_label,        # Single-cell labels
  st_count = mouse_vc_st$st_count,        # Spatial transcriptomics count matrix
  st_position = mouse_vc_st$st_position,  # Spatial positions of spots
  min_cont = 0.03                         # Minimum contamination threshold
)

# Calculate overall performance for refine_comp
res1 <- test_performances(
  test_spots_metadata_mtrx = as.matrix(mouse_vc_st$st_composition),  # Test spots metadata matrix
  spot_composition_mtrx = res$refine_comp                            # Refined composition matrix
)

# Calculate performance for each class in refine_comp
res2 <- each_class_performances(
  test_spots_metadata_mtrx = as.matrix(mouse_vc_st$st_composition),  # Test spots metadata matrix
  spot_composition_mtrx = res$refine_comp                            # Refined composition matrix
)

# Calculate overall performance for pred_comp
res3 <- test_performances(
  test_spots_metadata_mtrx = as.matrix(mouse_vc_st$st_composition),  # Test spots metadata matrix
  spot_composition_mtrx = res$pred_comp                             # Predicted composition matrix
)

# Calculate performance for each class in pred_comp
res4 <- each_class_performances(  # Use a different variable name to avoid overwriting res2
  test_spots_metadata_mtrx = as.matrix(mouse_vc_st$st_composition),  # Test spots metadata matrix
  spot_composition_mtrx = res$pred_comp                             # Predicted composition matrix
)
```
