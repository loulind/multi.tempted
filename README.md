
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

Then in containerized R-studio session, (1) open File \> Open Project
then select multi.tempted.Rproj (adds git tab) (2) once opened, run…

``` r
setwd('~/work')
library(devtools)
library(testthat)
library(roxygen2)
library(knitr)
```

<!-- README.md is generated from README.Rmd. Please edit that file -->

# multi.tempted

<!-- badges: start -->

<!-- badges: end -->

The goal of multi.tempted is to …

## Installation

You can install the development version of multi.tempted from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("loulind/multi.tempted")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(multi.tempted)
#> Loading required package: np
#> np 0.70-2
#> Examples and guides at https://jeffreyracine.github.io/gallery/
#> See also vignette("np_getting_started", package = "np")
#> Loading required package: ggplot2
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
