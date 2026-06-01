#' @title Run all major functions of multiTEMPTED
#' @description This function wraps functions...
#' \code{\link{format_tempted}},
#' \code{\link{svd_centralize}},
#' \code{\link{tempted}},
#' @param r Number of components (rank). Default 3.
#' @param smooth RKHS smoothing penalty. Larger = smoother temporal functions.
#'    Adjust by checking the smoothness of the estimated temporal loading function plot.
#'    Default 1e-8.
#' @param featuretables A sample by feature matrix. It is an input for \code{\link{format_tempted}}.
#' @param modalities A list of strings (names of modalities)
#'    Equals length of featuretables
#' @param timepoint The time stamp of each sample, matched with the rows of \code{featuretable}.
#'    It is an input for \code{\link{format_tempted}}.
#' @param subjectID The subject ID of each sample, matched with the rows of \code{featuretable}.
#'    It is an input for \code{\link{format_tempted}}.
#' @param threshold A threshold for feature filtering for microbiome data.
#'    Features with zero value percentage >= threshold will be excluded. Default is 0.95.
#'    It is an input for \code{\link{format_tempted}}.
#' @param pseudo A small number to add to all the counts before
#'    normalizing into proportions and log transformation.
#'    Default is 1/2 of the smallest non-zero value that is specific for each sample.
#'    This pseudo count is added for \code{transform=c("logcomp", "clr", "logit")}.
#'    It is an input for \code{\link{format_tempted}}.
#' @param transforms A list of the transformations applied to each modality
#'    \code{"logcomp"} for log of compositions.
#'    \code{"comp"} for compositions.
#'    \code{"ast"} for arcsine squared transformation.
#'    \code{"clr"} for central log ratio transformation.
#'    \code{"logit"} for logit transformation.
#'    \code{"none"} for no transformation.
#'    Default \code{transform="clr"} is recommended for microbiome data.
#'    For data that are already transformed, use \code{transform="none"}.
#'    It is an input for \code{\link{format_tempted}}.
#'    Equals length of featuretables
#' @param r Number of components to decompose into, i.e. rank of the CP type decomposition.
#'    Default is set to 3.
#'    It is an input for \code{\link{tempted}}.
#' @param smooth Smoothing parameter for RKHS norm.
#'    Larger means smoother temporal loading functions. Default is set to be 1e-8.
#'    Value can be adjusted depending on the dataset by checking the smoothness of the estimated temporal loading function in plot.
#'    It is an input for \code{\link{tempted}}.
#' @param interval The range of time points to ran the decomposition for.
#'    Default is set to be the range of all observed time points.
#'    User can set it to be a shorter interval than the observed range.
#'    It is an input for \code{\link{tempted}}.
#' @param resolution Number of time points to evaluate the value of the temporal loading function.
#'    Default is set to 101. It does not affect the subject or feature loadings. It is an input for \code{\link{tempted}}.
#' @param maxiter Maximum number of iteration. Default is 20. It is an input for \code{\link{tempted}}.
#' @param epsilon Convergence criteria for difference between iterations. Default is 1e-4. It is an input for \code{\link{tempted}}.
#' @param r_svd The number of ranks in the mean structure. Default is 1. It is an input for \code{\link{svd_centralize}}.
#' @return A list including all the input and output of functions \code{\link{format_tempted}}, \code{\link{svd_centralize}}, \code{\link{tempted}},
#' #' \code{\link{ratio_feature}}, and \code{\link{aggregate_feature}}.
#' \describe{
#'   \item{input}{All the input options of function \code{\link{tempted_all}}.}
#'   \item{datalist_raw}{Output of \code{\link{format_tempted}} with option \code{transform="none"}.}
#'   \item{datlist}{Output of \code{\link{format_tempted}}.}
#'   \item{mean_svd}{Output of \code{\link{svd_centralize}}.}
#'   \item{A_hat}{Subject loading, a subject by r matrix.}
#'   \item{B_hat}{Feature loading, a feature by r matrix.}
#'   \item{Zeta_hat}{Temporal loading function, a resolution by r matrix.}
#'   \item{time_Zeta}{The time points where the temporal loading function is evaluated.}
#'   \item{Lambda}{Eigen value, a length r vector.}
#'   \item{r_square}{Variance explained by each component. This is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using individual components.}
#'   \item{accum_r_square}{Variance explained by the first few components accumulated. This is the R-squared of the linear regression of the vectorized temporal tensor against the vectorized low-rank reconstruction using the first few components.}
#' }
#' @references
#' Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight R, Shenhav L, Zhang AR. (2023) \emph{Time-Informed Dimensionality Reduction for Longitudinal Microbiome Studies}. bioRxiv. doi: 10.1101/550749. \url{https://www.biorxiv.org/content/10.1101/550749}.
#' @export
#' @md
tempted_all <- function(featuretable, modalities, timepoint, subjectID,
                        threshold=0.95, pseudo=NULL, transform="clr",
                        r = 3, smooth=1e-6,
                        interval = NULL, resolution = 51,
                        maxiter=20, epsilon=1e-4, r_svd=1,
                        do_ratio=TRUE, pct_ratio=0.05, absolute=FALSE,
                        pct_aggregate=1, contrast=NULL){
  datlists <- lapply(1:M, function(m) {
    format_tempted(featuretable=featuretable[[m]], timepoint=timepoints[[m]],
                   subjectID=subjectID, threshold=threshold, pseudo=pseudo,
                   transform=transform)
    }
  )

  datalists_raw <- lapply(1:M, function(m) {
    format_tempted(featuretable=featuretable[[m]], timepoint=timepoints[[m]],
                   subjectID=subjectID, threshold=threshold, pseudo=pseudo,
                   transform="none")
    }
  )

  names(datlists) <- names(datalist_raw) <- modalities # assign modality names
  mean_svd <- svd_centralize(datlists, r_svd)
  res_tempted <- multi_tempted_decomp(datlists=mean_svd$datlists, r = r,
                                      smooth=smooth, interval = interval,
                                      resolution = resolution,
                                      maxiter=maxiter, epsilon=epsilon)
  # later: add metafeatures calc
  return(res_tempted)

}
