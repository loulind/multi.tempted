#' Integrated Personal Omics Profiling (iPOP)
#'
#' The example dataset is from a longitudinal exercise omics study. Subjects
#' were measured at five time intervals. Four omic modalities were profiled
#' at each visit: cytokine, metabolome, lipid, protein.
#'
#' @format A length-6 named list of data frames:
#' \describe{
#'   \item{meta_subj}{Subject metadata}
#'   \item{meta_visit}{Visit metadata}
#'   \item{cytokine}{Plasma cytokine panel}
#'   \item{metabolome}{Serum metabolomics}
#'   \item{lipid}{Lipidomics}
#'   \item{protein}{Proteomics}
#' }
#' @source \url{https://med.stanford.edu/snyderlab/ipop.html}
"ipop"
