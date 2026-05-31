#' @title Decomposition of temporal tensors across multiple modalities
#' @description
#' CP-type decomposition of M subject-by-feature-by-time tensors.
#' Each tensor ("modality") shares a subject loading with all others.
#' Run after formatting data with `format_multitempted()` &
#' after `svd_centralize` if centralize option is selected.
#' @param datlists A length-M named list of length-n lists of matrices.
#'   Each matrix is one subject: row 1 is sampling times, rows 2..p+1 are features.
#' @param r Number of components (rank). Default 3.
#' @param smooth RKHS smoothing penalty. Larger = smoother temporal functions.
#'   Default 1e-8.
#' @param interval Named list (length M) of length-2 vectors giving the time
#'   range to decompose for each modality. Default: full observed range per modality.
#' @param resolution Grid size for evaluating temporal loading functions. Default 101.
#' @param maxiter Maximum iterations per component. Default 20.
#' @param epsilon Convergence threshold on squared loading change. Default 1e-4.
#' @param weights Length-M numeric vector of non-negative modality weights (w_1,...,w_M)
#'   that scale the data-fit term for each modality in the objective function.
#'   Larger weight = greater emphasis on that modality when estimating the shared
#'   subject loading and each modality's temporal loading.
#'   Default: equal weights (all 1).
#' @return A list with:
#'   \describe{
#'     \item{A_hat}{Subject loading, n x r matrix (shared across modalities).}
#'     \item{B_hat}{Length-M list of feature loading matrices (p_m x r).}
#'     \item{Zeta_hat}{Length-M list of temporal loading matrices (resolution x r).}
#'     \item{time_Zeta}{Length-M list of time grids for Zeta (original time scale).}
#'     \item{Lambda}{M x r matrix of modality-specific scales.}
#'     \item{r_square}{M x r matrix. r_square[m, l] is the R-squared of component
#'       l's rank-1 reconstruction against the current residual for modality m
#'       (i.e. after deflating components 1..l-1).}
#'     \item{accum_r_square}{M x r matrix. accum_r_square[m, l] is the R-squared
#'       of the first l components' joint reconstruction against the original
#'       (pre-deflation) data for modality m.}
#'   }
#' @export
#' @md
multi_tempted_decomp <- function(datlists, r=3, smooth=1e-8, interval=NULL,
                                 resolution = 101, maxiter=20, epsilon=1e-4,
                                 weights=NULL) {
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

  # Default to equal weights
  if (is.null(weights)) weights <- rep(1, M)
  if (length(weights) != M) stop("'weights' must have length equal to the number of modalities (M).")
  if (any(weights < 0))    stop("'weights' must be non-negative.")

  # Initialize time intervals (rescale to [0,1], bin based on resolution, build kernel matrices)
  prep <- lapply(1:M, function(m) {
    interval_m <- if (!is.null(interval)) interval[[m]] else NULL
    init_time_intv(datlists[[m]], p[m], n, interval_m, resolution)
  })
  for (m in 1:M) datlists[[m]] <- prep[[m]]$datlist

  # Initialize output
  A <- matrix(0, nrow = n, ncol = r) # shared subject loading matrix
  B <- lapply(1:M, function(m) matrix(0, p[m], r)) # list of feature loading matrices
  Zeta <- lapply(1:M, function(m) matrix(0, resolution, r)) # list of time loading fns
  Lambda <- matrix(0, M, r)  # modality-specific scalings
  Rsq      <- matrix(0, M, r)  # per-modality R-squared, one value per (modality, component)
  accumRsq <- matrix(0, M, r)  # per-modality accumulated R-squared

  # Flatten original data for accumulated R-sqr tracking (one vector per modality)
  y0_per_modality <- lapply(1:M, function(m)
    flatten_features(datlists[[m]], p[m], prep[[m]]$tipos))
  X_accum <- vector("list", M)  # per-modality design matrix; grows one column per component

  # SEQUENTIAL ESTIMATION ALGORITHM
  for (l in 1:r) {
    message(sprintf("Estimating component %d of %d", l, r))

    # Current residual data (after subtracting components 1..l-1),
    # ...for component l's lambda regression and R-squared calc
    y_resid <- lapply(1:M, function(m)
      flatten_features(datlists[[m]], p[m], prep[[m]]$tipos))

    # STEP 1: Initialize subject and feature loadings
    a_hat  <- rep(1 / sqrt(n), n) # Subject loadings init'd to equal contribution
    b_hats <- lapply(1:M, function(m) # init'd as list of matrices of leftmost singular vectors of SVD
      init_b_hat(datlists[[m]], p[m], n))

    # STEP 2: Update a, b, and zeta until max iterations reached or it converges
    iter <- 0
    dif  <- 1
    while (iter <= maxiter && dif > epsilon) {
      # (a) Temporal loading: update zeta for each modality independently
      zeta_hats <- lapply(1:M, function(m)
        update_zeta(datlists[[m]], p[m], b_hats[[m]], a_hat, prep[[m]]$ind_vec,
                    prep[[m]]$Kmat, prep[[m]]$Kmat_output, smooth, weights[m]))

      # (b) Subject loading: update shared a_hat across modalities
      a_new <- update_a(datlists, p, b_hats, zeta_hats, prep, n, M, weights)
      dif <- sum((a_hat - a_new)^2)
      a_hat <- a_new

      # (c) Feature loading: update b for each modality independently
      for (m in 1:M) {
        b_new <- update_b(datlists[[m]], p[m], zeta_hats[[m]],
                          prep[[m]]$tipos, prep[[m]]$ti, a_hat, n)
        dif <- max(dif, sum((b_hats[[m]] - b_new)^2))
        b_hats[[m]] <- b_new
      }
      iter <- iter + 1
    }
    message(sprintf("  Converged: dif = %.2e after %d iterations", dif, iter))

    # STEP 3: Est'm modality-specific scales (lambda), R-squared, & remove l'th component
    A[, l] <- a_hat # record estimated a
    for (m in 1:M) {
      B[[m]][, l]    <- b_hats[[m]]   # record estimated b
      Zeta[[m]][, l] <- zeta_hats[[m]] # record estimated zeta

      lm_result    <- compute_lambda(y_resid[[m]], datlists[[m]], p[m],
                                     a_hat, b_hats[[m]], zeta_hats[[m]],
                                     prep[[m]]$tipos, prep[[m]]$ti, n)
      Lambda[m, l] <- lm_result$lambda

      # Per-modality R-squared: component l against current residual
      Rsq[m, l]      <- compute_rsq(y_resid[[m]], lm_result$x_m)
      # Per-modality accumulated R-squared: components 1..l against original data
      X_accum[[m]]   <- cbind(X_accum[[m]], lm_result$x_m)
      accumRsq[m, l] <- compute_rsq(y0_per_modality[[m]], X_accum[[m]])
    }

    # Remove l'th component contribution from datlists
    for (m in 1:M) {
      datlists[[m]] <- update_datlist(datlists[[m]], p[m], a_hat, b_hats[[m]],
                                      zeta_hats[[m]], Lambda[m, l],
                                      prep[[m]]$tipos, prep[[m]]$ti, n)
    }
  }
  # STEP 4: After est'm all r components, re-est'm modality scales
  Lambda <- reestimate_lambda(y0_per_modality, A, B, Zeta, prep, p, M, n, r)

  # Revise signs so results are comparable (no sign switching)
  signs  <- revise_signs(A, B, Zeta, Lambda, r, M)
  A <- signs$A
  B <- signs$B
  Zeta <- signs$Zeta
  Lambda <- signs$Lambda

  # Re-map time intervals back to original time scales
  time_Zeta <- lapply(1:M, function(m) {
    grid <- seq(prep[[m]]$interval[1], prep[[m]]$interval[2], length.out = resolution)
    grid * (prep[[m]]$input_time_range[2] - prep[[m]]$input_time_range[1]) +
      prep[[m]]$input_time_range[1]
  })
  names(time_Zeta) <- names(datlists)

  return(list(
    A_hat = A,
    B_hat = B,
    Zeta_hat = Zeta,
    time_Zeta = time_Zeta,
    Lambda = Lambda,
    r_square = Rsq,
    accum_r_square = accumRsq
  ))
}






# ----------------- HELPER FUNCTIONS ------------------

# -----(1) Preprocessing functions------

#' Rescale time to [0,1], bin samples to grid, build Bernoulli kernel matrices.
#'
#' @returns Named list: datlist (time-rescaled), ti, tipos, ind_vec, Kmat,
#'   Kmat_output, input_time_range, interval.
#' @noRd
init_time_intv <- function(datlist_m, p_m, n, interval_m, resolution) {

  timestamps_all <- sort(unique(unlist(lapply(datlist_m, function(s) s[1, ]))))
  input_time_range <- c(timestamps_all[1], timestamps_all[length(timestamps_all)])

  if (is.null(interval_m)) interval_m <- input_time_range

  # Rescale all time points (including interval bounds) to [0, 1].
  rescale <- function(t) {(t - input_time_range[1]) / (input_time_range[2] - input_time_range[1])}
  for (i in 1:n) datlist_m[[i]][1, ] <- rescale(datlist_m[[i]][1, ])
  interval_m <- rescale(interval_m)

  # Map each sample to a grid index (0 = outside interval).
  ti <- lapply(1:n, function(i) {
    idx <- 1 + round((resolution - 1) *
                       (datlist_m[[i]][1, ] - interval_m[1]) /
                       (interval_m[2] - interval_m[1]))
    idx[idx <= 0 | idx > resolution] <- 0
    idx
  })

  tipos <- lapply(1:n, function(i) ti[[i]] > 0)

  # ind_vec maps each sample (across all subjects) to its subject index.
  ind_vec <- unlist(lapply(1:n, function(i) rep(i, ncol(datlist_m[[i]]))))

  tm <- unlist(lapply(datlist_m, function(s) s[1, ]))
  grid <- seq(interval_m[1], interval_m[2], length.out = resolution)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(grid, tm)

  return(list(
    datlist = datlist_m,
    ti = ti,
    tipos = tipos,
    ind_vec = ind_vec,
    Kmat = Kmat,
    Kmat_output = Kmat_output,
    input_time_range = input_time_range,
    interval = interval_m
  ))
}


#' Concatenate in-interval feature values across subjects into a single vector.
#' @noRd
flatten_features <- function(datlist_m, p_m, tipos_m) {
  y <- NULL
  for (i in seq_along(datlist_m)) {
    y <- c(y, as.vector(t(datlist_m[[i]][2:(p_m + 1), tipos_m[[i]]])))
  }

  return(y)
}





# ------(2) Kernel functions------

#' Functional regression with RKHS penalty (Bernoulli kernel ridge regression).
#'
#' @param Ly Length-n list; Ly[[i]] is the projected time series for subject i.
#' @param a_hat Length-n subject loading vector.
#' @param ind_vec Integer vector mapping each sample to its subject (1..n).
#' @param Kmat Kernel matrix between all observed time points.
#' @param Kmat_output Kernel matrix between resolution grid and observed points.
#' @param smooth RKHS penalty weight (C_mK in the objective).
#' @param weight Modality weight (w_m in the objective). Scales the data term
#'   relative to the RKHS penalty, so the normal equation becomes
#'   (w_m * K_a + smooth * I) beta = w_m * cvec. Default 1.
#' @noRd
freg_rkhs <- function(Ly, a_hat, ind_vec, Kmat, Kmat_output, smooth = 1e-8,
                      weight = 1) {
  K <- Kmat
  for (i in seq_along(Ly)) {
    K[ind_vec == i, ] <- K[ind_vec == i, ] * a_hat[i]^2
  }
  cvec <- unlist(Ly)
  beta <- solve(weight * K + smooth * diag(ncol(K))) %*% (weight * cvec)

  return(Kmat_output %*% beta)
}


#' Bernoulli kernel between vectors x and y.
#'
#' @references
#' Han, R., Shi, P. and Zhang, A.R. (2023) Guaranteed functional tensor singular
#' value decomposition. JASA. doi:10.1080/01621459.2022.2153689.
#' @noRd
bernoulli_kernel <- function(x, y) {
  k1_x <- x - 0.5
  k1_y <- y - 0.5
  k2_x <- 0.5 * (k1_x^2 - 1/12)
  k2_y <- 0.5 * (k1_y^2 - 1/12)
  xy <- abs(x %*% t(rep(1, length(y))) - rep(1, length(x)) %*% t(y))
  k4_xy <- 1/24 * ((xy - 0.5)^4 - 0.5 * (xy - 0.5)^2 + 7/240)

  return(k1_x %*% t(k1_y) + k2_x %*% t(k2_y) - k4_xy + 1)
}





# -------(3) Initialization function--------

#' Initialize b as the first left singular vector of the mode-2 unfolding.
#' @noRd
init_b_hat <- function(datlist_m, p_m, n) {
  data_unfold <- NULL
  for (i in 1:n) {
    data_unfold <- cbind(data_unfold, datlist_m[[i]][2:(p_m + 1), ])
  }

  return(svd(data_unfold, nu = 1, nv = 0)$u[, 1])
}






# --------(4) Updating functions----------

#' Update temporal loading (zeta) for one modality via RKHS regression.
#'
#' Projects each subject's data onto b, then fits the result as a smooth
#' function of time (penalised by the RKHS norm). Returns the unit-norm
#' vector on the resolution grid.
#' @noRd
update_zeta <- function(datlist_m, p_m, b_hat, a_hat, ind_vec,
                        Kmat, Kmat_output, smooth, weight = 1) {
  Ly <- lapply(seq_along(datlist_m), function(i)
    a_hat[i] * as.numeric(b_hat %*% datlist_m[[i]][2:(p_m + 1), ]))

  zeta <- freg_rkhs(Ly, a_hat, ind_vec, Kmat, Kmat_output, smooth, weight)

  return(zeta / sqrt(sum(zeta^2)))
}


#' Update shared subject loading (a_hat) by pooling signal across all modalities.
#'
#' For each subject i the optimal (unnormalised) a is:
#'   numerator   = sum_m  w_m * b_m^T X_i^(m) zeta_i^(m)
#'   denominator = sum_m  w_m * || zeta_i^(m) ||^2
#' @noRd
update_a <- function(datlists, p, b_hats, zeta_hats, prep, n, M, weights) {
  a_tilde <- numeric(n)

  for (i in 1:n) {
    num <- 0
    den <- 0
    for (m in 1:M) {
      in_range <- prep[[m]]$tipos[[i]]
      grid_idx <- prep[[m]]$ti[[i]][in_range]
      zeta_i <- zeta_hats[[m]][grid_idx]
      num <- num + weights[m] * as.numeric(
        b_hats[[m]] %*% datlists[[m]][[i]][2:(p[m] + 1), in_range] %*% zeta_i)
      den <- den + weights[m] * sum(zeta_i^2)
    }
    a_tilde[i] <- num / den
  }

  return(a_tilde / sqrt(sum(a_tilde^2)))
}


#' Update feature loading (b_hat) for one modality.
#'
#' Solves the weighted least-squares problem for b given fixed a and zeta.
#' Returns the unit-norm solution.
#' @noRd
update_b <- function(datlist_m, p_m, zeta_hat, tipos_m, ti_m, a_hat, n) {
  num <- matrix(0, p_m, n)
  denom <- numeric(n)

  for (i in 1:n) {
    in_range <- tipos_m[[i]]
    grid_idx <- ti_m[[i]][in_range]
    zeta_i <- zeta_hat[grid_idx]
    num[, i] <- datlist_m[[i]][2:(p_m + 1), in_range] %*% zeta_i
    denom[i] <- sum(zeta_i^2)
  }

  b_tilde <- as.numeric(num %*% a_hat) / as.numeric(denom %*% (a_hat^2))

  return(b_tilde / sqrt(sum(b_tilde^2)))
}


#' Estimate modality-specific scale (lambda) via no-intercept regression.
#'
#' Regresses the vectorised residual data against the rank-1 reconstruction
#' a_hat[i] * outer(b_hat, zeta_hat[grid_idx]).
#'
#' @importFrom stats lm
#' @returns Named list: lambda (scalar), x_m (reconstruction vector).
#' @noRd
compute_lambda <- function(y_m, datlist_m, p_m, a_hat, b_hat, zeta_hat,
                           tipos_m, ti_m, n) {
  x_m <- NULL
  for (i in 1:n) {
    grid_idx <- ti_m[[i]][tipos_m[[i]]]
    x_m <- c(x_m, as.vector(t(a_hat[i] * outer(b_hat, zeta_hat[grid_idx]))))
  }
  lambda <- as.numeric(lm(y_m ~ x_m - 1)$coefficients)

  return(list(
    lambda = lambda,
    x_m = x_m
  ))
}





# --------(5) Post estimation algorithm functions---------

#' Re-estimate all component lambdas jointly from the original (pre-deflation) data.
#'
#' For each modality, regresses y0 against the rank-r reconstruction matrix
#' [x_1 | x_2 | ... | x_r] without intercept.
#' @importFrom stats lm
#' @noRd
reestimate_lambda <- function(y0_per_modality, A, B, Zeta, prep, p, M, n, r) {
  Lambda <- matrix(0, M, r, dimnames = list(NULL, colnames(A)))

  for (m in 1:M) {
    X_m <- NULL
    for (l in 1:r) {
      x_l <- NULL
      for (i in 1:n) {
        grid_idx <- prep[[m]]$ti[[i]][prep[[m]]$tipos[[i]]]
        x_l <- c(x_l, as.vector(t(A[i, l] * outer(B[[m]][, l], Zeta[[m]][grid_idx, l]))))
      }
      X_m <- cbind(X_m, x_l)
    }
    Lambda[m, ] <- as.numeric(lm(y0_per_modality[[m]] ~ X_m - 1)$coefficients)
  }
  return(Lambda)
}


#' Remove the contribution of one component from one modality's data.
#' @noRd
update_datlist <- function(datlist_m, p_m, a_hat, b_hat, zeta_hat,
                           lambda_ml, tipos_m, ti_m, n) {
  for (i in 1:n) {
    in_range <- which(tipos_m[[i]])
    grid_idx <- ti_m[[i]][tipos_m[[i]]]
    datlist_m[[i]][2:(p_m + 1), in_range] <-
      datlist_m[[i]][2:(p_m + 1), in_range] -
      lambda_ml * a_hat[i] * outer(b_hat, zeta_hat[grid_idx])
  }

  return(datlist_m)
}


#' R-squared of regressing y on X (no intercept).
#' @importFrom stats lm
#' @noRd
compute_rsq <- function(y, X) {
  return(summary(lm(y ~ X - 1))$r.squared)
}






# ---------(6) Revising Signs----------

#' Canonicalise signs so that loadings are interpretable.
#'
#' All sign corrections are absorbed into B[[m]][,l], which acts as the
#' per-modality "sign sink". This guarantees three compatible conventions
#' without circular conflicts:
#'
#'   1. Lambda[m,l] >= 0  (absorb sign into B).
#'   2. sum(Zeta[[m]][,l]) >= 0  (absorb sign into B).
#'   3. sum(A[,l]) >= 0  (shared; absorb sign into all B[[m]][,l]).
#'
#' Note: sum(B[[m]][,l]) is NOT guaranteed non-negative. Enforcing
#' sum(B) >= 0 simultaneously with sum(Zeta) >= 0 is impossible in
#' general for a shared-A multi-modality model, because doing so creates
#' a circular flip (fixing one undoes the other).
#'
#' @noRd
revise_signs <- function(A, B, Zeta, Lambda, r, M) {
  for (l in 1:r) {

    # 1. Ensure Lambda >= 0; absorb sign into B (preserves product).
    for (m in 1:M) {
      if (Lambda[m, l] < 0) {
        Lambda[m, l] <- -Lambda[m, l]
        B[[m]][, l]  <- -B[[m]][, l]
      }
    }

    # 2. Ensure sum(Zeta) >= 0 per modality; absorb sign into B (preserves product).
    for (m in 1:M) {
      sgn <- sign(sum(Zeta[[m]][, l]))
      if (sgn == 0) sgn <- 1
      if (sgn < 0) {
        Zeta[[m]][, l] <- -Zeta[[m]][, l]
        B[[m]][, l]    <- -B[[m]][, l]
      }
    }

    # 3. Ensure sum(A) >= 0 (shared loading); absorb sign into all B[[m]] (preserves product).
    sgn_A <- sign(sum(A[, l]))
    if (sgn_A == 0) sgn_A <- 1
    if (sgn_A < 0) {
      A[, l] <- -A[, l]
      for (m in 1:M) B[[m]][, l] <- -B[[m]][, l]
    }
  }

  return(list(
    A = A,
    B = B,
    Zeta = Zeta,
    Lambda = Lambda
  ))
}
