#' De-noised temporal tensors for all modalities using estimated loadings.
#'
#' Reconstructs the rank-r approximation (plus mean structure if supplied)
#' for each modality evaluated on the resolution grid.
#'
#' @param res_decomp Output of \code{\link{multi_tempted_decomp}} or
#'   \code{\link{multitempted_all}}.
#' @param mean_svd Output of \code{\link{svd_centralize}}, or \code{NULL}.
#' @return Length-M named list of (n x p_m x resolution) arrays, with the
#'   third-dimension names set to the original-scale time grid.
#' @noRd
tdenoise <- function(res_decomp, mean_svd = NULL) {
  A <- res_decomp$A_hat
  n <- nrow(A)
  r <- ncol(A)
  M <- length(res_decomp$B_hat)
  resolution <- nrow(res_decomp$Zeta_hat[[1]])
  mod_names  <- names(res_decomp$B_hat)

  tensors <- lapply(seq_len(M), function(m) {
    p_m <- nrow(res_decomp$B_hat[[m]])
    tensor <- array(0, dim = c(n, p_m, resolution))

    if (!is.null(mean_svd)) {
      mean_m <- mean_svd$A_tilde[[m]] %*%
        t(sweep(mean_svd$B_tilde[[m]], 2, mean_svd$lambda_tilde[[m]], "*"))
      for (k in seq_len(resolution)) tensor[, , k] <- tensor[, , k] + mean_m
    }

    for (l in seq_len(r)) {
      tensor <- tensor +
        res_decomp$Lambda[m, l] *
        (A[, l] %o% res_decomp$B_hat[[m]][, l] %o% res_decomp$Zeta_hat[[m]][, l])
    }

    dimnames(tensor)[[3]] <- as.character(res_decomp$time_Zeta[[m]])
    tensor
  })

  return(setNames(tensors, mod_names))
}


#' @title Aggregate features using feature loadings (multi-modality)
#' @description For each modality, computes a weighted sum of features using the
#'   feature loading of each component (and any user-supplied contrasts) as
#'   weights. Returns both observed and de-noised estimated meta-features.
#' @param res_decomp Output of \code{\link{multi_tempted_decomp}} or
#'   \code{\link{multitempted_all}}.
#' @param mean_svd Output of \code{\link{svd_centralize}}, or \code{NULL}.
#' @param datlists Length-M named list of formatted (and transformed) data,
#'   one datlist per modality (output of \code{\link{format_tempted}} per modality).
#' @param pct Fraction of features to aggregate, ranked by absolute loading.
#'   Default 1 (all features). Setting \code{pct = 0.01} uses the top 1\%.
#' @param contrast An r x K matrix whose columns combine components via a linear
#'   contrast of the feature loadings, or \code{NULL}.
#' @return A list:
#' \describe{
#'   \item{metafeature_aggregate}{Data frame with columns \code{value},
#'     \code{subID}, \code{timepoint}, \code{PC}, and \code{modality}, giving
#'     the observed weighted-sum meta-feature at each subject's actual sample
#'     times.}
#'   \item{metafeature_aggregate_est}{Same structure as
#'     \code{metafeature_aggregate} but evaluated on the resolution grid from
#'     the de-noised tensor.}
#'   \item{contrast}{The contrast matrix from input.}
#'   \item{toppct}{Length-M named list of logical matrices (p_m x (r + K))
#'     indicating which features are included in each component (and contrast).}
#' }
#' @export
#' @md
aggregate_feature <- function(res_decomp, mean_svd = NULL, datlists,
                              pct = 1, contrast = NULL) {
  M <- length(res_decomp$B_hat)
  mod_names <- names(res_decomp$B_hat)
  r <- ncol(res_decomp$A_hat)
  subj_names <- rownames(res_decomp$A_hat)

  # Build augmented feature loading matrices (PC columns + contrast columns)
  B_data_list <- lapply(seq_len(M), function(m) {
    B_m <- as.data.frame(res_decomp$B_hat[[m]])
    if (!is.null(contrast)) {
      cont_m <- res_decomp$B_hat[[m]] %*% contrast
      colnames(cont_m) <- paste0("Contrast", seq_len(ncol(contrast)))
      B_m <- cbind(B_m, cont_m)
    }
    B_m
  })

  toppct_list <- lapply(seq_len(M), function(m) {
    apply(abs(B_data_list[[m]]), 2, function(x) x >= quantile(x, 1 - pct))
  })
  names(toppct_list) <- mod_names

  # Observed meta-features
  metafeature_aggregate <- NULL
  for (m in seq_len(M)) {
    B_m <- as.matrix(B_data_list[[m]])
    top_m <- toppct_list[[m]]
    datlist_m <- datlists[[m]]

    datlist_agg <- lapply(datlist_m, function(x) t(B_m * top_m) %*% x[-1, ])

    for (i in seq_along(datlist_agg)) {
      n_pc <- nrow(datlist_agg[[i]])
      n_ti <- ncol(datlist_agg[[i]])
      tmp <- data.frame(
        value = as.vector(datlist_agg[[i]]),
        subID = names(datlist_agg)[i],
        timepoint = rep(datlist_m[[i]][1, ], each = n_pc),
        PC = rep(rownames(datlist_agg[[i]]), n_ti),
        modality = mod_names[m]
      )
      metafeature_aggregate <- rbind(metafeature_aggregate, tmp)
    }
  }

  # Estimated meta-features from de-noised tensors
  tensor_list <- tdenoise(res_decomp, mean_svd)
  metafeature_aggregate_est <- NULL

  for (m in seq_len(M)) {
    B_m <- as.matrix(B_data_list[[m]])
    top_m <- toppct_list[[m]]
    tensor_m <- tensor_list[[m]]   # n x p_m x resolution
    time_grid <- res_decomp$time_Zeta[[m]]

    # apply over (subject, time): fn receives the p_m feature vector.
    # Force 3-D so the n_cols_B=1 (r=1) case doesn't get dimension-dropped.
    n_cols_B <- ncol(B_m)
    tensor_agg <- array(
      apply(tensor_m, c(1, 3), function(x) as.numeric(t(B_m * top_m) %*% x)),
      dim = c(n_cols_B, dim(tensor_m)[1], dim(tensor_m)[3])
    )
    # tensor_agg: n_cols_B x n x resolution

    for (l in seq_len(r)) {
      tmp <- data.frame(
        value = as.vector(tensor_agg[l, , ]),
        subID = rep(subj_names, length(time_grid)),
        timepoint = rep(time_grid, each = length(subj_names)),
        PC = colnames(B_data_list[[m]])[l],
        modality = mod_names[m]
      )
      metafeature_aggregate_est <- rbind(metafeature_aggregate_est, tmp)
    }
  }

  return(list(
    metafeature_aggregate = metafeature_aggregate,
    metafeature_aggregate_est = metafeature_aggregate_est,
    contrast = contrast,
    toppct = toppct_list
  ))
}


#' @title Log ratio of top vs bottom features (multi-modality)
#' @description For each modality, selects the top- and bottom-ranked features
#'   by component loading and returns the log ratio of their summed raw
#'   abundances. Designed for longitudinal microbiome count data; may not be
#'   meaningful for other data types.
#' @param res_decomp Output of \code{\link{multi_tempted_decomp}} or
#'   \code{\link{multitempted_all}}.
#' @param datlists_raw Length-M named list of raw (untransformed,
#'   \code{transform = "none"}) data, one datlist per modality.
#' @param pct Fraction of features to select. Default 0.05 (5\%).
#' @param absolute If \code{TRUE}, rank by |loading| and use sign to determine
#'   top/bottom. If \code{FALSE} (default), rank by loading value directly.
#' @param contrast An r x K contrast matrix, or \code{NULL}.
#' @return A list:
#' \describe{
#'   \item{metafeature_ratio}{Data frame with columns \code{value},
#'     \code{subID}, \code{timepoint}, \code{PC}, and \code{modality}.}
#'   \item{contrast}{The contrast matrix from input.}
#'   \item{toppct}{Length-M named list of logical matrices indicating the
#'     features used as the numerator (positive-loading features).}
#'   \item{bottompct}{Length-M named list of logical matrices indicating the
#'     features used as the denominator (negative-loading features).}
#' }
#' @export
#' @md
ratio_feature <- function(res_decomp, datlists_raw,
                          pct = 0.05, absolute = FALSE, contrast = NULL) {
  M <- length(res_decomp$B_hat)
  mod_names <- names(res_decomp$B_hat)

  B_data_list <- lapply(seq_len(M), function(m) {
    B_m <- as.data.frame(res_decomp$B_hat[[m]])
    if (!is.null(contrast)) {
      cont_m <- res_decomp$B_hat[[m]] %*% contrast
      colnames(cont_m) <- paste0("Contrast", seq_len(ncol(contrast)))
      B_m <- cbind(B_m, cont_m)
    }
    B_m
  })

  toppct_list <- vector("list", M)
  bottompct_list <- vector("list", M)
  for (m in seq_len(M)) {
    B_m <- B_data_list[[m]]
    if (!absolute) {
      toppct_list[[m]] <- apply(B_m,  2, function(x) x > quantile(x, 1 - pct) & x > 0)
      bottompct_list[[m]] <- apply(-B_m, 2, function(x) x > quantile(x, 1 - pct) & x > 0)
    } else {
      toppct_list[[m]] <- apply(B_m, 2, function(x) abs(x) > quantile(abs(x), 1 - pct) & x > 0)
      bottompct_list[[m]] <- apply(B_m, 2, function(x) abs(x) > quantile(abs(x), 1 - pct) & x < 0)
    }
  }
  names(toppct_list) <- mod_names
  names(bottompct_list) <- mod_names

  metafeature_ratio <- NULL

  for (m in seq_len(M)) {
    top_m <- toppct_list[[m]]
    bot_m <- bottompct_list[[m]]
    datlist_m <- datlists_raw[[m]]
    n_pc <- ncol(as.matrix(B_data_list[[m]]))

    pseudo <- min(sapply(datlist_m, function(x) {
      y <- x[-1, ]
      min(y[y != 0])
    })) / 2

    datlist_ratio <- lapply(datlist_m, function(x) {
      tt <- t(top_m) %*% x[-1, ]
      bb <- t(bot_m) %*% x[-1, ]
      log((tt + pseudo) / (bb + pseudo))
    })

    for (i in seq_along(datlist_ratio)) {
      n_ti <- ncol(datlist_ratio[[i]])
      tmp <- data.frame(
        value = as.vector(datlist_ratio[[i]]),
        subID = names(datlist_ratio)[i],
        timepoint = rep(datlist_m[[i]][1, ], each = n_pc),
        PC = rep(rownames(datlist_ratio[[i]]), n_ti),
        modality = mod_names[m]
      )
      metafeature_ratio <- rbind(metafeature_ratio, tmp)
    }
  }

  return(list(
    metafeature_ratio = metafeature_ratio,
    contrast = contrast,
    toppct = toppct_list,
    bottompct = bottompct_list
  ))
}
