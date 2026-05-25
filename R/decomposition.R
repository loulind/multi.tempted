#' @title Decomposition of temporal tensors
#' @description
#' CP-type decomposition of multiple 3d tensors.
#' Each tensor represents a different modality (subject x feature x time).
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
#' @param smooth Smoothing parameter for RKHS norm.
#' Larger ==> smoother temporal loading functions. Default is 1e-8.
#' Check the smoothness of the estimated temporal loading function plot to adjust.
#'
#' @param interval The range of time points to run the decomposition for.
#' Default is set to be the range of all observed time points.
#' User can set it to be a shorter interval than the observed range.
#'
#' @param resolution Number of time points to evaluate the value of the temporal loading function.
#' Default is set to 101. It does not affect the subject or feature loadings.
#'
#' @param maxiter Maximum number of iteration. Default is 20.
#'
#' @param epsilon Convergence criteria for difference between iterations. Default is 1e-4.
#'
#' @return The estimations of the loadings for each modality.
#'
#' \describe{
#'   \item{A_hat}{Subject loading, a subject by r matrix.}
#'   \item{B_hat}{Feature loading, a feature by r matrix.}
#'   \item{Phi_hat}{Temporal loading function, a resolution by r matrix.}
#'   \item{time_Phi}{The time points where the temporal loading function is evaluated.}
#'   \item{Lambda}{Eigenvalue, a length r vector.}
#'   \item{r_square}{Variance explained by each component. This is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using individual components.}
#'   \item{accum_r_square}{Variance explained by the first few components accumulated. This is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using the first few components.}
#' }
#' @export
#'
#' @examples
#'
multi_tempted_decomp <- function(datlists, r=3, smooth=1e-8, interval=NULL,
                                 resolution = 101, maxiter=20, epsilon=1e-4) {
  if (!(length(unique(lengths(datlists))) == 1)) {
    stop("All lists of matrices must be same length")
  }

  # Initialize intermediate variables
  M <- length(datlists)  # number modalities
  n <- length(datlists[[1]])  # number subjects
  A <- matrix(0, nrow = n, ncol = r) # subject loadings
  B <- vector(mode = "list", length = M)
  p <- vector(mode = "numeric", length = M)
  for (m in 1:M) {
    p[m] <- sapply(datlists[[m]], nrow)  # number features per modality
    if (!(length(unique(p)) == 1)) {
      stop(paste("Modality", names(datlists)[m], "has inconsistent feature counts across subjects"))
    }
    B[[m]] <- matrix(0, p[m], r)
  }



  # Calculate each component and remove contribution from feature values
  for (l in 1:r) {
    message(sprintf("Calculating Component %d", l))

    # STEP 1: Initialize subject and feature loadings per modality
    ## (i) Subject loadings, a, set to l'th row of A matrix
    a_hat <- rep(1/sqrt(n), n)

    ## (ii) Feature loadings, b, init'd as list of matrices of leftmost singular vectors of SVD
    b_hat <- init_b(datlists, M, p)

    # STEP 2: Update loadings until max iterations reached or it converges
    t <- 0
    dif <- 1
    for (m in 1:M) {
      message(sprintf("...modality %d", m))
      while(t<=maxiter & dif>epsilon){
        ## (i) Time-varying function, Zeta
        Ly <- list()
        for (i in 1:n){
          Ly <- c(Ly, list(a_hat[i]*as.numeric(b_hat[[m]]%*%datlists[[m]][[i]][2:(p+1),])))
        }
        zeta_hat <- update_zeta(Ly, a_hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
        zeta_hat <- zeta_hat / sqrt(sum(zeta_hat^2))

        ## (ii) Subject loadings
        update_a()

        ## (iii) Feature loadings
        update_b()
      }

      # STEP 3: Remove contribution of current component; repeat steps 1-2 for all r components
      centralize()
    }

  }

  # STEP 4: Estimate modality-specific scales
}





# HELPER FUNCTIONS

#' Initialize feature loading's vector, b
#'
#' @param datlists Length M named list of length n lists of matrices.
#' @param M number total modalities
#' @param p length M vector containing number of features, p_m, for modality m
#' @returns list of M length p_m vectors
#'
#' @noRd No user-side documentation
init_b <- function(datlists, M, b_hat, p) {
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
