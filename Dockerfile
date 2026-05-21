FROM rocker/verse:4.6.0

USER root

RUN R -e "install.packages(c('usethis', 'devtools', 'testthat', 'roxygen2'), repos='https://cloud.r-project.org')"
