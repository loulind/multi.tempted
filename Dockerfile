FROM rocker/verse:4.6.0

USER root

RUN R -e "install.packages(c('usethis'), repos='https://cloud.r-project.org')"
