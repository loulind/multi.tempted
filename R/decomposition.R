#' @title Decomposition of multiple temporal tensors
#' @description
#' Main function of multiTEMPTED
#'
#' @param datlists A length M named list of length n lists of matrices.
#' Each named list element represents a modality.
#' Each matrix represents a subject. Columns represent sample number.
#' The first rows represent sampling time points.
#' Row 2 through row (number features in modality m) + 1 represent feature values
#'
#' @param r Number components to decompose the M 3d tensors
#' (ie, rank of CP-type decomposition). Default is r=3.
#'
#' @returns
#' @export
#'
#' @examples
#'
multi_tempted <- function(datlists, r=3) {
  if (!(length(unique(lengths(datlists))) == 1)) {
    stop("All lists of matrices must be same length")
  }

  M <- length(datlists)  # number modalities
  n <- length(datlists[[1]])  # number subjects
  for (m in 1:M) {
    p <- sapply(datlists[[m]], nrow)  # number features per modality
    if (!(length(unique(p)) == 1)) {
      stop(paste("Modality", datlists[[m]], "has inconsistent feature counts across subjects"))
    }
  }

  # STEP 1: Initialize subject and feature loadings per modality
  ## (a) Subject loadings init'd with equal contribution
  a <- rep(1/sqrt(n), times = n)

  ## (b) Feature loadings init'd as matrix of leftmost singular vectors of SVD
  b <- vector("list", length = M)
  for (m in 1:M) {
    b[m] <- init_b(datlists[[m]], p[[m]])
  }

  # Calculate each component and remove contribution from feature values
  for (l in 1:r) {
    # STEP 2: Sequentially estimate loadings
    ## (a) Time loadings


    ## (b) Subject loadings


    ## (c) Feature loadings

    # STEP 3: Remove contribution of current component & repeat steps 1-2 for all components
  }

  # STEP 4: Estimate modality-specific scales
}

init_b <- function(datlist, p) {
  # 1. Matrication of datalist
  # 2. Perform SVD
  # 3. Extract leftmost singular vector
  return(leftmost_sv)
}

update_zeta <- function() {  # updates modality-specific time loadings
  # Kernel ridge regression code
}

update_a <- function() {  # updates cross-modality shared subject loading

}

update_b <- function() {  # updates feature loadings

}
