#' Integrated Personal Omics Profiling (iPOP)
#'
#' The example dataset is from a longitudinal exercise omics study. Subjects
#' were measured at five time intervals. Four omic modalities were profiled
#' at each visit: Cytokine, Metablome, Lipid, Protein.
#'
#' @format A length-6 named list of data frames:
#' \describe{
#'   \item{meta_subj}{Subject metadata}
#'   \item{meta_visit}{Visit metadata}
#'   \item{cytokine}{Plasma cytokine panel}
#'   \item{metablome}{Serum metablomics}
#'   \item{lipid}{Lipidomics}
#'   \item{protein}{Proteomics}
#' }
#' @source \url{https://med.stanford.edu/snyderlab/ipop.html}
"ipop"
