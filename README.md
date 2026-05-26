multiTEMPTED Vignette
================

If working in Docker container, run the following code…

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

Then in containerized R-studio session:

1.  USER: rstudio, PS: 123
2.  open File \> Open Project then select multi.tempted.Rproj (adds git
    tab)
3.  once opened, run…

``` r
setwd('~/work')
library(devtools)
library(testthat)
library(roxygen2)
library(knitr)
```

<!-- README.md is generated from README.Rmd. Please edit that file -->

# multiTEMPTED

<!-- badges: start -->

<!-- badges: end -->

The goal of multiTEMPTED is to to implement the statistical method
TEMPoral TEnsor Decomposition (TEMPTED) generalized to multiple
modalities.

Package dependencies: R (\>= 4.5.3), np (\>= 0.60-17), ggplot2 (\>=
3.4.0), methods (\>= 4.2.1)

Run time \_\_\_\_?\_\_\_\_

You can cite this paper for using TEMPTED: \_\_\_\_?\_\_\_\_

The statistical theories behind TEMPTED can be found in this paper:

\_\_\_\_?\_\_\_\_

\_\_\_\_?\_\_\_\_

## Installation

You can install the development version of multi.tempted from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("loulind/multi.tempted")
```

\_\_\_\_?\_\_\_\_\_

## Load packages for this vignette

## Read the example data

## Running TEMPTED for different formats of data

## Run TEMPTED for Microbiome Count Data (Straightforward Way)

### ***?steps?***\_

## Run TEMPTED for Microbiome Compositional Data (Straightforward Way)

### ***?steps?***\_

## Run TEMPTED for General Form of Multivariate Longitudinal Data (Straightforward Way)

### ***?steps?***\_

## Run TEMPTED in Customized Way

### ***?steps?***\_

## Transferring TEMPTED result from training to testing data

### ***?steps?***\_
