#' @title Run all major steps of multiTEMPTED
#' @description
#' Wrapper that runs the three core steps of multiTEMPTED in sequence:
#' \enumerate{
#'   \item \code{\link{format_tempted}} — transforms and formats each modality's
#'         feature table into the list-of-matrices input format.
#'   \item \code{\link{svd_centralize}} — optionally removes the rank-r mean
#'         structure from each modality so the decomposition captures
#'         time-varying rather than constant effects.
#'   \item \code{\link{multi_tempted_decomp}} — CP-type decomposition estimating
#'         the shared subject loading and modality-specific feature and temporal
#'         loadings.
#' }
#' @param featuretables A length-M **named** list of sample-by-feature matrices,
#'   one per modality. The list names become modality labels throughout all output.
#'   If unnamed, modalities are labelled \code{"modality1"}, \code{"modality2"}, etc.
#' @param timepoints Either a single numeric vector of time stamps shared by all
#'   modalities (rows of every feature table must be in the same sample order),
#'   or a length-M list of numeric vectors matched row-wise to each
#'   \code{featuretables[[m]]}.
#' @param subjectID Either a single character/integer vector of subject IDs shared
#'   by all modalities, or a length-M list of such vectors matched row-wise to
#'   each \code{featuretables[[m]]}. All modalities must ultimately contain the
#'   same set of subjects.
#' @param threshold Scalar (applied to all modalities) or length-M numeric vector.
#'   Features with zero-value percentage above this threshold are excluded.
#'   Default 0.95. Passed to \code{\link{format_tempted}}.
#' @param pseudo \code{NULL}, a scalar, or a length-M list. Small constant added
#'   before log-type transformations. \code{NULL} (default) uses half the
#'   minimum non-zero value per sample. Passed to \code{\link{format_tempted}}.
#' @param transforms Character string (applied to all modalities) or length-M
#'   character vector specifying the transformation for each modality.
#'   Options: \code{"clr"} (default), \code{"logcomp"}, \code{"comp"},
#'   \code{"ast"}, \code{"logit"}, \code{"lfb"}, \code{"none"}.
#'   Passed to \code{\link{format_tempted}}.
#' @param r Number of components (rank of CP decomposition). Default 3.
#'   Passed to \code{\link{multi_tempted_decomp}}.
#' @param smooth RKHS smoothing penalty for the temporal loading functions.
#'   Larger = smoother. Default 1e-8.
#'   Passed to \code{\link{multi_tempted_decomp}}.
#' @param interval Length-M list of length-2 vectors giving the time range to
#'   decompose for each modality, or \code{NULL} (default) to use each
#'   modality's full observed range.
#'   Passed to \code{\link{multi_tempted_decomp}}.
#' @param resolution Grid size for evaluating temporal loading functions.
#'   Default 101. Passed to \code{\link{multi_tempted_decomp}}.
#' @param maxiter Maximum iterations per component. Default 20.
#'   Passed to \code{\link{multi_tempted_decomp}}.
#' @param epsilon Convergence threshold on squared loading change. Default 1e-4.
#'   Passed to \code{\link{multi_tempted_decomp}}.
#' @param centralize Logical. If \code{TRUE} (default), runs
#'   \code{\link{svd_centralize}} before decomposition to remove the rank-\code{r_svd}
#'   mean structure from each modality.
#' @param r_svd Rank of mean structure to remove in \code{\link{svd_centralize}}.
#'   Only used when \code{centralize = TRUE}. Default 1.
#' @param weights Length-M non-negative numeric vector of modality weights for
#'   the decomposition objective. \code{NULL} (default) gives equal weight to
#'   all modalities. Passed to \code{\link{multi_tempted_decomp}}.
#' @return A named list with the following elements:
#' \describe{
#'   \item{datlists}{Length-M list of formatted (and transformed) data, output
#'     of \code{\link{format_tempted}} before centralisation.}
#'   \item{mean_svd}{Output of \code{\link{svd_centralize}}, or \code{NULL} if
#'     \code{centralize = FALSE}.}
#'   \item{A_hat}{Shared subject loading matrix (n x r).}
#'   \item{B_hat}{Length-M list of feature loading matrices (p_m x r).}
#'   \item{Zeta_hat}{Length-M list of temporal loading matrices (resolution x r).}
#'   \item{time_Zeta}{Length-M list of time grids for \code{Zeta_hat} on the
#'     original time scale.}
#'   \item{Lambda}{M x r matrix of modality-specific scales.}
#'   \item{r_square}{M x r matrix of per-modality R-squared values, one per
#'     (modality, component).}
#'   \item{accum_r_square}{M x r matrix of accumulated per-modality R-squared
#'     values across the first l components.}
#' }
#' @references
#' Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight R,
#' Shenhav L, Zhang AR. (2023) \emph{Time-Informed Dimensionality Reduction for
#' Longitudinal Microbiome Studies}. bioRxiv. doi: 10.1101/550749.
#' \url{https://www.biorxiv.org/content/10.1101/550749}.
#' @export
#' @md
multitempted_all <- function(featuretables, timepoints, subjectID,
                             threshold = 0.95, pseudo = NULL, transforms = "clr",
                             r = 3, smooth = 1e-8,
                             interval = NULL, resolution = 101,
                             maxiter = 20, epsilon = 1e-4,
                             centralize = TRUE, r_svd = 1,
                             weights = NULL) {

  M <- length(featuretables)
  if (M < 1) stop("'featuretables' must contain at least one modality.")
  if (is.null(names(featuretables))) { # auto names modalities if none specified
    names(featuretables) <- paste0("modality", seq_len(M))
  }

  # User can designate length 1 or length M list for parameters
  if (!is.list(timepoints)) timepoints <- replicate(M, timepoints, simplify = FALSE)
  if (!is.list(subjectID)) subjectID  <- replicate(M, subjectID, simplify = FALSE)
  if (length(transforms) == 1) transforms <- rep(transforms, M)
  if (length(threshold)  == 1) threshold  <- rep(threshold,  M)
  if (!is.list(pseudo)) pseudo <- replicate(M, pseudo, simplify = FALSE)

  if (length(timepoints) != M) stop("'timepoints' must be length 1 or length M.")
  if (length(subjectID) != M) stop("'subjectID' must be length 1 or length M.")
  if (length(transforms) != M) stop("'transforms' must be length 1 or length M.")
  if (length(threshold) != M) stop("'threshold' must be length 1 or length M.")

  # Format data
  datlists <- lapply(seq_len(M), function(m) {
    format_tempted(featuretable = featuretables[[m]],
                   timepoint = timepoints[[m]],
                   subjectID = subjectID[[m]],
                   threshold = threshold[m],
                   pseudo = pseudo[[m]],
                   transform = transforms[m])
  })
  names(datlists) <- names(featuretables)

  # Verify all modalities share the same subject set; reorder to a common order.
  subj_sets <- lapply(datlists, names)
  ref_subjs <- subj_sets[[1]]
  for (m in seq_len(M)) {
    if (!setequal(subj_sets[[m]], ref_subjs)) {
      stop(sprintf(
        "Modality '%s' has different subjects than modality '%s' after formatting.",
        names(datlists)[m], names(datlists)[1]))
    }
    datlists[[m]] <- datlists[[m]][ref_subjs]  # enforce consistent ordering
  }

  # Remove mean structure (if specified)
  if (centralize) {
    mean_svd       <- svd_centralize(datlists, r = r_svd)
    datlists_decomp <- mean_svd$datlists
  } else {
    mean_svd        <- NULL
    datlists_decomp <- datlists
  }

  # Decomposition
  res_decomp <- multi_tempted_decomp(datlists   = datlists_decomp,
                                     r          = r,
                                     smooth     = smooth,
                                     interval   = interval,
                                     resolution = resolution,
                                     maxiter    = maxiter,
                                     epsilon    = epsilon,
                                     weights    = weights)

  return(list(
    datlists       = datlists,
    mean_svd       = mean_svd,
    A_hat          = res_decomp$A_hat,
    B_hat          = res_decomp$B_hat,
    Zeta_hat       = res_decomp$Zeta_hat,
    time_Zeta      = res_decomp$time_Zeta,
    Lambda         = res_decomp$Lambda,
    r_square       = res_decomp$r_square,
    accum_r_square = res_decomp$accum_r_square
  ))
}
