#' @title Decomposition of temporal tensors
#' @description
#' CP-type decomposition of M 3d tensors (subject x feature x time).
#' Performed after formatting data!
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
multi_tempted_decomp <- function(datlists, r=3) {
  if (!(length(unique(lengths(datlists))) == 1)) {
    stop("All lists of matrices must be same length")
  }

  M <- length(datlists)  # number modalities
  n <- length(datlists[[1]])  # number subjects
  for (m in names(datlists)) {
    p <- sapply(datlists[[m]], nrow)  # number features per modality
    if (!(length(unique(p)) == 1)) {
      stop(paste("Modality", m, "has inconsistent feature counts across subjects"))
    }
  }

  # Calculate time interval and scale to 0,1

  # Calculate each component and remove contribution from feature values
  for (l in 1:r) {
    message(sprintf("Calculate the %dth Component", l))

    # STEP 1: Initialize subject and feature loadings per modality
    ## (i) Subject loadings, a, init'd with equal contribution
    a <- rep(1/sqrt(n), times = n)

    ## (ii) Feature loadings, b, init'd as matrix of leftmost singular vectors of SVD
    b <- vector("list", length = M)
    init_b(datlists, b, M, p)

    # STEP 2: Sequentially estimate loadings
    ## (i) Time loadings
    update_zeta()

    ## (ii) Subject loadings
    update_a()

    ## (iii) Feature loadings
    update_b()

    # STEP 3: Remove contribution of current component; repeat steps 1-2 for all r components
    centralize()
  }

  # STEP 4: Estimate modality-specific scales
}





# HELPER FUNCTIONS

#' Initialize feature loading's vector, b
#'
#' @param datlists Length M named list of length n lists of matrices.
#' @param M number total modalities
#' @param b list of M length p_m vectors
#' @param p length M vector containing number of features, p_m, for modality m
#' @returns list of M length p_m vectors
#'
#' @noRd No user-side documentation
init_b <- function(datlists, M, b, p) {
  b_hat <- vector(mode = "list", length = M)
  for (m in 1:M) {
    data_unfold <- NULL
    for (i in 1:n) {
      data_unfold = cbind(data_unfold, datlists[[m]][[i]][2:(p+1),])
    }
    b.intitials <- svd(data_unfold, nu=r, nv=r)$u
    b_hat[[m]] <- b.intitials[,1]
  }
  return(b_hat)
}

update_zeta <- function() {  # updates modality-specific time loadings
  # Kernel ridge regression code (RKHS penalty term)
}

update_a <- function() {  # updates cross-modality shared subject loading
  # 1. formula for updating a
  # 2. scale a by inverse norm of a
}

update_b <- function() {  # updates feature loadings
  # 1. formula for updating b
  # 2. scale b by inverse norm of b
}

svd_centralize <- function() {

}
