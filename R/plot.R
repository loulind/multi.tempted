#' @title Plot nonparametric smoothed mean and error bands of features versus time
#' @description Plots the smoothed mean and error bands for multiple features,
#'   grouped by a factor variable. Each feature is shown as a separate facet.
#' @param feature_mat A sample by feature matrix. The features can be original
#'   features, meta-features, log ratios, or any variables of interest.
#' @param time_vec A vector of time points matched to the rows of
#'   \code{feature_mat}.
#' @param group_vec A factor variable indicating group membership of samples,
#'   matched to the rows of \code{feature_mat}.
#' @param coverage Coverage rate for the error band. Default 0.95.
#' @param bws Bandwidth for the smoothing lines and error bands. A larger value
#'   means a smoother line. Default \code{NULL} uses \code{np::npreg()} with AIC
#'   bandwidth selection.
#' @param nrow Number of rows for \code{ggplot2::facet_wrap()}. Default 1.
#' @return A ggplot2 object.
#' @importFrom np npreg
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon ylab facet_wrap
#' @importFrom stats qnorm setNames
#' @export
#' @md
plot_feature_summary <- function(feature_mat, time_vec, group_vec,
                                 coverage = 0.95, bws = NULL, nrow = 1) {
  nfeature <- ncol(feature_mat)
  if (!is(group_vec, "factor")) group_vec <- as.factor(group_vec)
  group_level <- levels(group_vec)
  CI_length <- -qnorm((1 - coverage) / 2)

  if (is.null(colnames(feature_mat))) stop("feature_mat needs to have column names!")

  time_all <- NULL
  mean_all <- NULL
  merr_all <- NULL
  feature_all <- NULL
  group_all <- NULL

  for (jj in 1:nfeature) {
    for (ii in 1:length(group_level)) {
      ind <- group_vec == group_level[ii]
      if (is.null(bws)) {
        model_np <- npreg(feature_mat[ind, jj] ~ time_vec[ind],
                          regtype = "ll", bwmethod = "cv.aic")
      } else {
        model_np <- npreg(feature_mat[ind, jj] ~ time_vec[ind], bws = bws,
                          regtype = "ll", bwmethod = "cv.aic")
      }
      time_eval <- as.vector(t(model_np$eval))
      ord <- order(time_eval)
      time_all <- c(time_all,  sort(time_eval))
      mean_all <- c(mean_all,  model_np$mean[ord])
      merr_all <- c(merr_all,  model_np$merr[ord])
      feature_all <- c(feature_all, rep(colnames(feature_mat)[jj], length(time_eval)))
      group_all <- c(group_all,   rep(group_level[ii],           length(time_eval)))
    }
  }

  group_all <- factor(group_all, levels = group_level)
  tab_summary <- data.frame(time_all = time_all, mean_all = mean_all,
                            merr_all = merr_all, group_all = group_all,
                            feature_all = feature_all)
  .data <- NULL

  ggplot(data = tab_summary,
         aes(x = .data$time_all, y = .data$mean_all,
             group = .data$group_all, color = .data$group_all)) +
    geom_line() +
    geom_ribbon(aes(ymin = .data$mean_all - CI_length * .data$merr_all,
                    ymax = .data$mean_all + CI_length * .data$merr_all,
                    color = .data$group_all, fill = .data$group_all),
                linetype = 2, alpha = 0.3) +
    ylab(paste0("mean +/- ", round(CI_length, 2), "*se")) +
    facet_wrap(~ .data$feature_all, scales = "free", nrow = nrow)
}


#' @title Plot top feature loadings as horizontal bar charts (multi-modality)
#' @description For each modality, selects the features whose absolute loading
#'   exceeds the \code{1 - pct} quantile for each component and displays them as
#'   horizontal bars: negative loadings in red pointing left, positive loadings
#'   in blue pointing right. Within each component panel, features are arranged
#'   from most-negative at the top down to least-negative, then most-positive
#'   down to least-positive, so the bars form two "wedges" opening from the
#'   centre. Each modality is returned as a separate ggplot2 object with one
#'   facet per component.
#' @param res Output of \code{\link{multi_tempted_decomp}} or
#'   \code{\link{multitempted_all}}.
#' @param pct Fraction of features to display per component, ranked by absolute
#'   loading. Default 0.05 (top 5 percent).
#' @param xlim Length-2 numeric vector giving the x-axis limits.
#'   Default \code{c(-0.5, 0.5)}.
#' @return A length-M named list of ggplot2 objects, one per modality. Each
#'   plot is faceted by component (one panel per PC).
#' @seealso \code{\link{multitempted_all}}, \code{\link{multi_tempted_decomp}}.
#' @importFrom ggplot2 ggplot aes geom_col geom_vline scale_fill_manual scale_y_discrete coord_cartesian facet_wrap labs theme_bw theme element_text
#' @export
#' @md
plot_feature_loading <- function(res, pct = 0.05, xlim = c(-0.5, 0.5)) {
  mod_names <- names(res$B_hat)
  PC_names <- colnames(res$B_hat[[1]])

  lapply(setNames(mod_names, mod_names), function(mod) {
    B_m <- res$B_hat[[mod]]  # p_m x r

    # Build a long data frame covering all components for this modality.
    all_data <- do.call(rbind, lapply(seq_along(PC_names), function(l) {
      pc <- PC_names[l]
      loading <- B_m[, pc]
      thresh <- quantile(abs(loading), 1 - pct)
      keep <- abs(loading) >= thresh

      vals <- loading[keep]
      feats <- rownames(B_m)[keep]

      neg_ord <- order(vals[vals < 0]) # most negative first
      pos_ord <- order(vals[vals >= 0], decreasing = TRUE) # most positive first

      ordered_feats <- c(feats[vals < 0][neg_ord], feats[vals >= 0][pos_ord])
      ordered_vals <- c(vals[vals < 0][neg_ord],  vals[vals >= 0][pos_ord])

      if (length(ordered_feats) == 0) return(NULL)

      data.frame(
        # Unique label per (feature, component) so free_y facets don't bleed.
        feat_label = paste0(ordered_feats, "__", l),
        feat_name = ordered_feats,
        value = ordered_vals,
        sign = ifelse(ordered_vals < 0, "negative", "positive"),
        component = pc,
        display_order = seq_along(ordered_feats), # 1 = top of panel
        stringsAsFactors = FALSE
      )
    }))

    if (is.null(all_data) || nrow(all_data) == 0) return(NULL)

    # Factor level order: within each component, ggplot y-axis runs bottom-to-top,
    # so sort by descending display_order so that display_order=1 (most negative)
    # ends up as the last (topmost) level in each component's subset.
    all_data_sorted <- all_data[order(all_data$component, -all_data$display_order), ]
    all_data$feat_factor <- factor(all_data$feat_label,
                                   levels = all_data_sorted$feat_label)

    .data <- NULL
    ggplot(all_data, aes(x = .data$value, y = .data$feat_factor,
                         fill = .data$sign)) +
      geom_col(width = 0.5, color="grey30", linewidth=0.25) +
      geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey30") +
      scale_fill_manual(values = c(negative = "red3", positive = "steelblue"),
                        guide = "none") +
      scale_y_discrete(labels = function(x) sub("__\\d+$", "", x)) +
      coord_cartesian(xlim = xlim) +
      facet_wrap(~ component, scales = "free_y", ncol = 1) +
      labs(title = mod, x = "Feature loading", y = NULL) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  })
}

#' @title Plot subject loadings in PC space
#' @description Scatter plot of the shared subject loading matrix (\code{A_hat})
#'   for two chosen components, optionally coloured by a grouping variable.
#' @param res Output of \code{\link{multi_tempted_decomp}} or
#'   \code{\link{multitempted_all}}.
#' @param group Optional subject x 2 data frame: first column is subject ID,
#'   second column is the grouping variable (e.g. treatment arm, sex). The
#'   second column's name is used as the legend title. \code{NULL} (default)
#'   plots all subjects in one colour.
#' @param pcs Length-2 integer vector selecting which two components to plot.
#'   Default \code{c(1, 2)}.
#' @return A ggplot2 object.
#' @seealso \code{\link{multitempted_all}}, \code{\link{multi_tempted_decomp}}.
#' @importFrom ggplot2 ggplot aes geom_point scale_color_brewer labs theme_minimal theme element_text
#' @export
#' @md
plot_subject_loading <- function(res, group = NULL, pcs = c(1, 2)) {
  A <- if (is.matrix(res) || is.data.frame(res)) as.matrix(res) else res$A_hat
  r <- ncol(A)

  if (length(pcs) != 2)  stop("'pcs' must be a length-2 integer vector.")
  if (any(pcs < 1 | pcs > r)) stop(sprintf("'pcs' values must be between 1 and r (%d).", r))

  pc_cols <- colnames(A)[pcs]   # e.g. c("PC1", "PC2")

  plot_df        <- as.data.frame(A[, pcs, drop = FALSE])
  plot_df$subID  <- rownames(A)

  group_label <- NULL
  if (!is.null(group)) {
    group_label    <- colnames(group)[2]
    colnames(group) <- c("subID", "group")
    plot_df        <- merge(plot_df, group, by = "subID")
  }

  .data <- NULL

  if (!is.null(group)) {
    p <- ggplot(plot_df,
                aes(x = .data[[pc_cols[1]]], y = .data[[pc_cols[2]]],
                    color = .data$group)) +
      scale_color_brewer(palette = "Set1") +
      labs(color = group_label)
  } else {
    p <- ggplot(plot_df,
                aes(x = .data[[pc_cols[1]]], y = .data[[pc_cols[2]]]))
  }

  p +
    geom_point(size = 4, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste0("Principal Components ", pcs[1], " and ", pcs[2]),
      x     = paste0("PC ", pcs[1]),
      y     = paste0("PC ", pcs[2])
    ) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}


#' @title Plot smoothed mean and error bands of meta-features versus time
#' @description Plots the smoothed mean and error bands of meta-features grouped
#'   by a factor variable. For multi-modality output, returns a named list of
#'   ggplot2 objects, one per modality.
#' @param metafeature \code{metafeature_ratio} from \code{\link{ratio_feature}}
#'   or \code{\link{multitempted_all}}, or \code{metafeature_aggregate} from
#'   \code{\link{aggregate_feature}} or \code{\link{multitempted_all}}. Must
#'   contain columns \code{value}, \code{subID}, \code{timepoint}, \code{PC},
#'   and \code{modality}.
#' @param group A subject x 2 data frame: first column is subject ID, second
#'   column is group membership.
#' @param coverage Coverage rate for the error band. Default 0.95.
#' @param bws Bandwidth for smoothing. Default \code{NULL} uses AIC selection.
#' @param nrow Number of rows for \code{ggplot2::facet_wrap()}. Default 1.
#' @return A named list of ggplot2 objects, one per modality.
#' @seealso \code{\link{aggregate_feature}}, \code{\link{ratio_feature}},
#'   \code{\link{multitempted_all}}.
#' @export
#' @md
plot_metafeature <- function(metafeature, group,
                             coverage = 0.95, bws = NULL, nrow = 1) {
  colnames(group) <- c("subID", "group")
  tab <- merge(metafeature, group, by = "subID")
  mod_names <- unique(tab$modality)

  plots <- lapply(setNames(mod_names, mod_names), function(mod) {
    tab_mod <- tab[tab$modality == mod, ]

    reshape_tab <- reshape(tab_mod[, c("subID", "timepoint", "group", "PC", "value")],
                           idvar = c("subID", "timepoint", "group"),
                           v.names = "value",
                           timevar = "PC",
                           direction = "wide")
    CC <- grep("^value\\.", colnames(reshape_tab))
    colnames(reshape_tab)[CC] <- paste0(
      sub("^value\\.", "", colnames(reshape_tab)[CC]), " (", mod, ")")

    feature_mat <- reshape_tab[, CC, drop = FALSE]
    time_vec <- reshape_tab$timepoint
    group_vec <- factor(reshape_tab$group)

    plot_feature_summary(feature_mat, time_vec, group_vec,
                         coverage = coverage, bws = bws, nrow = nrow)
  })

  return(plots)
}


#' @title Plot the temporal loading functions (multi-modality)
#' @description Plots the temporal loading functions (\code{Zeta_hat}) estimated
#'   by \code{\link{multi_tempted_decomp}}, faceted by modality and coloured by
#'   component.
#' @param res Output of \code{\link{multi_tempted_decomp}} or
#'   \code{\link{multitempted_all}}.
#' @param r Number of components to plot. Default: all components in \code{res}.
#' @param ... Additional aesthetics passed to \code{ggplot2::geom_line()}.
#' @return A ggplot2 object.
#' @seealso \code{\link{multitempted_all}}, \code{\link{multi_tempted_decomp}}.
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap
#' @export
#' @md
plot_time_loading <- function(res, r = NULL, ...) {
  M <- length(res$Zeta_hat)
  mod_names <- names(res$Zeta_hat)
  if (is.null(r)) r <- ncol(res$Zeta_hat[[1]])

  plot_data <- do.call(rbind, lapply(1:M, function(m) {
    Zeta_m <- res$Zeta_hat[[m]][, 1:r, drop = FALSE]
    ntime <- nrow(Zeta_m)
    data.frame(
      timepoint = rep(res$time_Zeta[[m]], r),
      value = as.vector(Zeta_m),
      component = factor(rep(1:r, each = ntime)),
      modality = mod_names[m],
      stringsAsFactors = FALSE
    )
  }))

  .data <- NULL

  ggplot(plot_data,
         aes(x = .data$timepoint, y = .data$value, color = .data$component)) +
    geom_line(aes(...)) +
    facet_wrap(~ .data$modality, scales = "free_x")
}
