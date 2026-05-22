#' @title Decomposition of multiple temporal tensors
#' @description
#' Main function of multiTEMPTED
#'
#' @param datlists A length M named list of length n lists of matrices.
#' Each matrix represents a subject. Columns represent sample number.
#' The first rows represent sampling time points.
#' Row 2 through row (# features for modality m) + 1  represent feature values
#'
#' @param r Number components to decompose M 3d tensors (ie, rank of CP-type decomposition)
#' Default is r=3.
#'
#' @returns
#' @export
#'
#' @examples
#'
multi_tempted_all <- function(datlists, r=3, ) {
  if (!(length(unique(lengths(datlists))) == 1)) {
    stop("All lists of matrices must be same length", call. = FALSE)
  }

  M <- length(datlists)  # number modalities
  n <- length(datlists[1])  # number subjects

  for (l in 1:r) { # iterate over desired number of components
    # Step 1: Initialize subject and feature loadings per modality
    # (a) Subject loadings init'd with equal contribution
    a <- init_a(n)

    # (b) Feature loadings init'd as leftmost singular vector of SVD
    for (m in 1:M)
    b <- init_b(datlists)

    # Step 2: Sequentially estimate loadings
    # (a) Time loadings
    # (b) Subject loadings
    # (c) Feature loadings

    # Step 3: Remove contribution of current component & repeat steps 1-2 for all components

    # Step 4: Estimate modality-specific scales
  }
}

init_a <- function(n) {

}

init_b <- function(datlists) {

}
