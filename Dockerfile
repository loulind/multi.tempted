FROM rocker/verse:4.5.3

USER root

RUN R -e "install.packages(c('usethis', 'devtools', 'testthat', 'roxygen2', 'knitr', 'rmarkdown'), repos='https://cloud.r-project.org')"
