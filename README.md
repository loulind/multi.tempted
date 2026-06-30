multiTEMPTED Vignette
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Development setup (Docker)

If working in a Docker container, run the following:

``` zsh
# emulates amd64 if working on non-linux and builds container
docker build . --platform=linux/amd64 -t multi.tempted

# ports to 8787 for rocker container
docker run \
  --platform=linux/amd64 \
  -v $(pwd):/home/rstudio/work \
  -e PASSWORD=123 \
  -p 8787:8787 \
  multi.tempted
```

Then in the containerized RStudio session:

1.  USER: rstudio, PS: 123
2.  Open File \> Open Project and select `multi.tempted.Rproj`
3.  Once opened, run:

``` r
setwd('~/work')
library(devtools)
library(testthat)
library(roxygen2)
library(knitr)
```

------------------------------------------------------------------------

## Introduction of multiTEMPTED

This is a vignette for the R package `multi.tempted`, which implements
**multiTEMPTED** — a generalization of the TEMPoral TEnsor Decomposition
(TEMPTED) method to simultaneous multi-modality longitudinal data.

Where the original TEMPTED decomposes a single subject × feature × time
tensor, multiTEMPTED jointly decomposes M such tensors (one per
modality) into:

- A **shared subject loading matrix** A (subjects × r) — Per subject
  contributions to each component. This ties all modalities together.
- **Per-modality feature loading matrices** B_m (p_m × r) — which
  features drive each component within each modality.
- **Per-modality temporal loading functions** Zeta_m — smooth curves
  describing how each component evolves over time within each modality.
- **Per-modality scaling factors** Lambda (M × r) — the relative
  magnitude of each component within each modality.

Modalities may have different feature sets, different numbers of
measured features, and unaligned sampling times between subjects or
between modalities.

**Package dependencies:** R (\>= 4.5.0), np (\>= 0.60-17), ggplot2 (\>=
3.4.0), methods (\>= 4.2.1).

You can **cite this paper** for now using TEMPTED (the single-modality
precursor):

Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight
R, Shenhav L, Zhang AR. [Time-Informed Dimensionality Reduction for
Longitudinal Microbiome Studies.
bioRxiv.](https://doi.org/10.1101/2023.07.26.550749)

The statistical theory behind the functional tensor decomposition:

Han R, Shi P, Zhang AR. [Guaranteed Functional Tensor Singular Value
Decomposition. Journal of the American Statistical Association (2023):
1-13.](https://doi.org/10.1080/01621459.2022.2153689)

------------------------------------------------------------------------

## Installation

You can install the development version of `multi.tempted` from
[GitHub](https://github.com/loulind/multi.tempted) with:

``` r
# install.packages("pak")
pak::pak("loulind/multi.tempted")
```

------------------------------------------------------------------------

## Load packages for this vignette

``` r
library(tidyverse)
library(corrplot)
library(plotly)
library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)
library(magick)
library(pheatmap)

library(multi.tempted)
```

------------------------------------------------------------------------

## Read the example data

The example dataset is from a longitudinal exercise omics study (more
info can be found at <https://med.stanford.edu/snyderlab/ipop.html>).
Subjects were measured at five time intervals. Four modalities were
profiled at each visit:

- **cytokine** — plasma cytokine panel
- **metabolome** — serum metabolomics
- **lipid** — lipidomics
- **protein** — proteomics

All four modalities are provided pre-processed on a log₁₀ scale. Two
metadata objects accompany the data:

- `ipop[[1]]` (`meta_subj`): one row per subject with demographic
  variables.
- `ipop[[2]]` (`meta_visit`): one row per sample with `SubjectID` and
  `timepoint` (coded 1–5).

``` r
names(ipop)
# [1] "meta_subj" "meta_visit" "cytokine" "metabolome" "lipid" "protein"
```

``` r
# Stack the four modality tables into a named list
featuretables <- lapply(3:6, function(m) as.matrix(ipop[[m]]))
names(featuretables) <- names(ipop)[3:6]
M <- length(featuretables)

# Map visit codes to approximate minutes from baseline
timepoint <- as.vector(ipop[[2]]$timepoint)
timepoint <- case_when(
  timepoint == 1 ~ 0,
  timepoint == 2 ~ 12,
  timepoint == 3 ~ 25,
  timepoint == 4 ~ 40,
  timepoint == 5 ~ 70
)
timepoints <- rep(list(timepoint), M)
subjectID  <- rep(list(as.character(ipop[[2]]$SubjectID)), M)

# Group variables
group_sex <- ipop[[1]]$Sex # one column of sex for each subject
group_subID <- ipop[[1]][, c("subjectID", "Sex")] # two columns: subjectID and their corresponding sex
```

------------------------------------------------------------------------

## Running multiTEMPTED

### Run multiTEMPTED

A complete description of all parameters can be found in the function
documentation (`?multitempted_all`). Key parameters:

- `featuretables`: length-M named list of sample × feature matrices, one
  per modality.
- `timepoints`: length-M list of numeric vectors giving the sampling
  time of each row in the corresponding feature table.
- `subjectID`: length-M list of character vectors of subject IDs matched
  to rows of each feature table.
- `transforms`: transformation applied per modality before
  decomposition. `"clr"` (default) is recommended for raw microbiome
  counts; use `"none"` for data that are already on a log or otherwise
  pre-processed scale.
- `r`: number of components (rank of the CP-type decomposition). Default
  3.
- `smooth`: RKHS smoothing penalty for the temporal loading functions.
  Larger values produce smoother curves. Default 1e-8.
- `centralize`: if `TRUE` (default), removes the rank-`r_svd` mean
  structure from each modality before decomposition via
  `svd_centralize()`.
- `do_ratio`: if `TRUE` (default), computes log-ratio meta-features. Set
  to `FALSE` for data that are not raw counts.

**IMPORTANT NOTE:** As in matrix SVD, the signs of subject loadings,
feature loadings, and temporal loadings can be flipped in any consistent
combination. Interpretation should focus on the relative pattern of
values across subjects or features, not the absolute sign.

``` r
output <- multitempted_all(
  featuretables = featuretables,
  timepoints    = timepoints,
  subjectID     = subjectID,
  transforms    = "none", # data are already log10-transformed
  do_ratio      = FALSE,  # not raw counts
  r             = 3 # number computed components (will compute more later)
)
names(output)
```

------------------------------------------------------------------------

### Low-dimensional representation of subjects

The shared subject loading matrix `A_hat` (subjects × r) is the
multi-modality generalisation of a PCA score matrix: each row is one
subject’s position in the r-dimensional latent space learned jointly
from all modalities.

``` r
plot_subject_loading(output$A_hat, group = group_subID)
```

------------------------------------------------------------------------

### Plot the temporal loading functions

The temporal loading functions Zeta_m describe how each component
evolves over time within each modality. A peak in component l’s temporal
loading for modality m indicates that the contrast captured by component
l is strongest at that time point in that modality.

``` r
plot_time_loading(output) +
  labs(x = "Weeks from baseline", title = "Temporal loadings by modality")
```

------------------------------------------------------------------------

### Plot the feature loadings

Feature loadings B\_(m,l) rank the features of modality m by their
contribution to component l. The function below displays the top 1% of
features by absolute loading for each modality and component, with
negative loadings in red and positive loadings in blue.

``` r
plots_feat <- plot_feature_loading(output, pct = 0.01)

# display one modality at a time, e.g.:
plots_feat$cytokine
plots_feat$metabolome
```

------------------------------------------------------------------------

### Trajectories of individual features

The feature loading rankings can be used to identify biologically
relevant features to inspect directly. Below we plot the smoothed mean
trajectory of two cytokines (GLP-1 and Insulin) grouped by sex.

``` r
feat_mat <- ipop[["cytokine"]][, c("GLP1", "INSULIN")]

plot_feature_summary(
  feature_mat = feat_mat,
  time_vec    = timepoints[[1]],
  group_vec   = group_sex,
  bws         = 10
) + labs(x = "Weeks from baseline")
```

------------------------------------------------------------------------

### Subject trajectories (meta-features)

The feature loadings can be used as weights to aggregate all features of
a modality into a single “meta-feature” trajectory per subject per
component. This summarises each modality’s contribution to each
component over time and allows group-level trajectory comparisons.

`metafeature_aggregate` uses the observed data;
`metafeature_aggregate_est` uses the low-rank de-noised reconstruction.

``` r
# returns a named list of ggplots, one per modality
traj_plots <- plot_metafeature(output$metafeature_aggregate, group = group_subID)

# display one modality:
traj_plots$cytokine + labs(x = "Weeks from baseline")
```

------------------------------------------------------------------------

### Modality loading correlations

Within each modality, features that have similar loading profiles across
all r components will tend to respond together across the exercise
intervention. Computing pairwise Kendall correlations between feature
loading vectors reveals modules of co-responding features.

``` r
output2 <- multitempted_all(
  featuretables = featuretables,
  timepoints    = timepoints,
  subjectID     = subjectID,
  transforms    = "none", # data are already log10-transformed
  do_ratio      = FALSE,  # not raw counts
  r             = 10  # computing more components to compute correlations
)
```

``` r
# extract feature loadings per modality (transpose to r x p for cor())
cyto_loadings  <- t(as.matrix(output2$B_hat[["cytokine"]]))
metab_loadings <- t(as.matrix(output2$B_hat[["metabolome"]]))
lipid_loadings <- t(as.matrix(output2$B_hat[["lipid"]]))
prot_loadings  <- t(as.matrix(output2$B_hat[["protein"]]))

# within-modality feature correlation matrices
cyto_corr_mat  <- cor(cyto_loadings,  method = "kendall")
metab_corr_mat <- cor(metab_loadings, method = "kendall")
lipid_corr_mat <- cor(lipid_loadings, method = "kendall")
prot_corr_mat  <- cor(prot_loadings,  method = "kendall")

corrplot(cyto_corr_mat,  method = "color", type = "lower", tl.cex = 0.2, order = "hclust")
corrplot(metab_corr_mat, method = "color", type = "lower", tl.cex = 0.2, order = "hclust")
corrplot(lipid_corr_mat, method = "color", type = "lower", tl.cex = 0.2, order = "hclust")
corrplot(prot_corr_mat,  method = "color", type = "lower", tl.cex = 0.2, order = "hclust")
```

Cross-modality correlation matrices reveal features from different
modalities that share a loading profile — i.e., are driven by the same
latent components. (using pheatmap package to create rectangular ordered
heatmaps)

``` r
cyto_v_metab <- cor(cyto_loadings, metab_loadings, method = "kendall")
cyto_v_lipid <- cor(cyto_loadings, lipid_loadings, method = "kendall")
cyto_v_prot  <- cor(cyto_loadings, prot_loadings,  method = "kendall")
metab_v_lipid <- cor(metab_loadings, lipid_loadings, method = "kendall")
metab_v_prot  <- cor(metab_loadings, prot_loadings,  method = "kendall")
lipid_v_prot  <- cor(lipid_loadings, prot_loadings,  method = "kendall")

pheatmap(cyto_v_metab,
         clustering_method = "complete",
         color = colorRampPalette(c("red", "white", "blue"))(50),
         main = "Cytokine vs Metabolome Correlation Heatmap",
         fontsize_row = 5,
         fontsize_col = 1)
pheatmap(cyto_v_lipid,
         clustering_method = "complete",
         color = colorRampPalette(c("red", "white", "blue"))(50),
         main = "Cytokine vs Lipid Correlation Heatmap",
         fontsize_row = 5,
         fontsize_col = 4)
pheatmap(cyto_v_prot,
         clustering_method = "complete",
         color = colorRampPalette(c("red", "white", "blue"))(50),
         main = "Cytokine vs Protein Correlation Heatmap",
         fontsize_row = 5,
         fontsize_col = 4)
pheatmap(metab_v_lipid,
         clustering_method = "complete",
         color = colorRampPalette(c("red", "white", "blue"))(50),
         main = "Metablome vs Lipid Correlation Heatmap",
         fontsize_row = 1,
         fontsize_col = 4)
pheatmap(metab_v_prot,
         clustering_method = "complete",
         color = colorRampPalette(c("red", "white", "blue"))(50),
         main = "Metablome vs Protein Correlation Heatmap",
         fontsize_row = 1,
         fontsize_col = 4)
pheatmap(lipid_v_prot,
         clustering_method = "complete",
         color = colorRampPalette(c("red", "white", "blue"))(50),
         main = "Lipid vs Protein Correlation Heatmap",
         fontsize_row = 5,
         fontsize_col = 2)
```

------------------------------------------------------------------------

### Correlation network

Thresholding a within-modality correlation matrix yields a co-response
network where nodes are features and edges connect features whose
loading profiles exceed a correlation threshold. Edge colour indicates
the sign of the correlation (blue = positive, red = negative).

``` r
# choose a modality correlation matrix and threshold,

choose_corr_here <- metab_corr_mat # change "metab" to one of "cyto", "lipid", or "prot"

threshold       <- 0.70
adjacency_matrix <- ifelse(abs(choose_corr_here) >= threshold, choose_corr_here, 0)
diag(adjacency_matrix) <- 0

g <- graph_from_adjacency_matrix(
  adjacency_matrix,
  mode     = "undirected",
  weighted = TRUE,
  diag     = FALSE
)

plot(
  g,
  vertex.size  = 2,
  vertex.label = NA,
  edge.width   = abs(E(g)$weight) * 2,
  edge.color   = ifelse(E(g)$weight > 0, "steelblue", "tomato"),
  layout       = layout_with_fr(g, weights = abs(E(g)$weight))
)
```
