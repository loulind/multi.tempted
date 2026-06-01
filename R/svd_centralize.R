#' @title Remove the mean structure of the temporal tensors
#' @description
#' For each modality, this function averages each feature's value across all
#' observed time points for each subject, forming an n-by-p_m subject-by-feature
#' matrix. It then computes a rank-r SVD of that matrix and subtracts the
#' rank-r approximation from the full temporal tensor. The result is a tensor
#' whose constant-in-time (mean) component has been removed, so that the
#' subsequent decomposition captures time-varying structure rather than
#' dominant constant effects.
#' @param datlists A length-M named list of length-n lists of matrices.
#'   Each matrix is one subject: row 1 is sampling times, rows 2..p+1 are features.
#' @param r Number of ranks to remove from the mean structure. Default is 1.
#' @return A named list:
#' \describe{
#'   \item{datlists}{The centralised temporal tensors (same structure as input).}
#'   \item{A_tilde}{Length-M list of n x r subject loading matrices for the mean structure.}
#'   \item{B_tilde}{Length-M list of p_m x r feature loading matrices for the mean structure.}
#'   \item{lambda_tilde}{Length-M list of length-r singular value vectors for the mean structure.}
#' }
#' @references
#' Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight R,
#' Shenhav L, Zhang AR. (2023) \emph{Time-Informed Dimensionality Reduction for
#' Longitudinal Microbiome Studies}. bioRxiv. doi: 10.1101/550749.
#' \url{https://www.biorxiv.org/content/10.1101/550749}.
#' @export
#' @md
svd_centralize <- function(datlists, r = 1) {
  if (!(length(datlists) >= 1)) {
    stop("Must have a strictly positive number of modalities")
  }
  if (length(unique(lengths(datlists))) != 1) {
    stop("All modalities must have the same number of subjects.")
  }

  # Initialize data dimensions
  M <- length(datlists)  # number modalities
  n <- length(datlists[[1]])  # number subjects
  p <- sapply(1:M, function(m) { # list of number features per modality
    pm_vals <- sapply(datlists[[m]], nrow) - 1
    if (length(unique(pm_vals)) != 1) {
      stop(sprintf("Modality '%s' has inconsistent feature counts across subjects.",
                   names(datlists)[m]))
    }
    pm_vals[[1]]
  })

  # initialize output
  datlists_new <- datlists
  A_tilde <- vector("list", M)
  B_tilde <- vector("list", M)
  lambda_tilde <- vector("list", M)

  # Centralise each modality independently
  for (m in 1:M) {

    # Step 1: Average each feature over time for every subject -> n x p_m matrix
    mean_hat <- matrix(0, n, p[m])
    for (i in 1:n) {
      mean_hat[i, ] <- rowMeans(datlists[[m]][[i]][-1, , drop = FALSE])
    }

    # Step 2: Rank-r SVD of the subject-by-feature mean matrix
    mean_svd  <- svd(mean_hat, nu = r, nv = r)
    mean_rank_r <- mean_svd$u %*% t(mean_svd$v * mean_svd$d[1:r])

    # Step 3: Subtract the subject's mean profile from every time point
    for (i in seq_len(n)) {
      datlists_new[[m]][[i]][-1, ] <- datlists[[m]][[i]][-1, ] - mean_rank_r[i, ]
    }

    A_tilde[[m]] <- mean_svd$u
    B_tilde[[m]] <- mean_svd$v
    lambda_tilde[[m]] <- mean_svd$d[1:r]
  }

  names(A_tilde) <- names(datlists)
  names(B_tilde) <- names(datlists)
  names(lambda_tilde) <- names(datlists)

  return(list(
    datlists = datlists_new,
    A_tilde = A_tilde,
    B_tilde = B_tilde,
    lambda_tilde = lambda_tilde
  ))
}
