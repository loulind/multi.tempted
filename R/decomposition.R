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
#' @return A list with:
#'   \describe{
#'     \item{A_hat}{Subject loading, n x r matrix (shared across modalities).}
#'     \item{B_hat}{Length-M list of feature loading matrices (p_m x r).}
#'     \item{Zeta_hat}{Length-M list of temporal loading matrices (resolution x r).}
#'     \item{time_Zeta}{Length-M list of time grids for Zeta (original time scale).}
#'     \item{Lambda}{M x r matrix of modality-specific scales.}
#'     \item{r_square}{Variance explained per component (pooled across modalities).}
#'     \item{accum_r_square}{Accumulated variance explained by first l components.}
#'   }
#' @export
#' @md
multi_tempted_decomp <- function(datlists, r=3, smooth=1e-8, interval=NULL,
                                 resolution = 101, maxiter=20, epsilon=1e-4) {
    if (!(length(datlists) >= 1)) {
     stop("Must have a strictly positive number of modalities")
  }
    if (length(unique(lengths(datlists))) != 1) {
      stop("All modalities must have the same number of subjects.")
  }

  # Initialize data dimensions
  M <- length(datlists)  # number modalities
  n <- length(datlists[[1]])  # number subjects
  p <- sapply(1:M, function(m) { # list of features per modality
    pm_vals <- sapply(datlists[[m]], nrow) - 1
    if (length(unique(pm_vals)) != 1) {
      stop(sprintf("Modality '%s' has inconsistent feature counts across subjects.",
                   names(datlists)[m]))
    }
    pm_vals[[1]]
  })

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
  Rsq <- numeric(r)
  accumRsq <- numeric(r)

  # Flatten original data for accumulated R-squared tracking
  y0_per_modality <- lapply(1:M, function(m)
    flatten_features(datlists[[m]], p[m], prep[[m]]$tipos))
  y0_all  <- unlist(y0_per_modality)
  X_accum <- NULL  # grows one column per component

  # SEQUENTIAL ESTIMATION ALGORITHM
  for (l in 1:r) {
    message(sprintf("Estimating component %d of %d", l, r))

    # Current residual data (after subtracting components 1..l-1), for this
    # component l's lambda regression and R-squared calc
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
                    prep[[m]]$Kmat, prep[[m]]$Kmat_output, smooth))

      # (b) Subject loading: update shared a_hat across modalities
      a_new <- update_a(datlists, p, b_hats, zeta_hats, prep, n, M)
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

    # STEP 3: Est'm modality-specific scales (lambda) & Remove l'th component
    # Recording component values
    A[, l] <- a_hat # record estimated a
    x_comp <- NULL # rank-1 reconstruction vectorized across all modalities
    for (m in 1:M) {
      B[[m]][, l] <- b_hats[[m]]  # record estimated b
      Zeta[[m]][, l] <- zeta_hats[[m]]  # record estimated zeta

      lm_result <- compute_lambda(y_resid[[m]], datlists[[m]], p[m],
                                     a_hat, b_hats[[m]], zeta_hats[[m]],
                                     prep[[m]]$tipos, prep[[m]]$ti, n)
      Lambda[m, l] <- lm_result$lambda # update lambda
      x_comp <- c(x_comp, lm_result$x_m)
    }

    # R-squared for this component and accumulated total (pooled across modalities)
    Rsq[l] <- compute_rsq(unlist(y_resid), x_comp)
    X_accum <- cbind(X_accum, x_comp)
    accumRsq[l] <- compute_rsq(y0_all, X_accum)

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

# ------- (1) Preprocessing functions ----------
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

# --------- (5) Kernel functions -------------
#' Functional regression with RKHS penalty (Bernoulli kernel ridge regression).
#'
#' @param Ly Length-n list; Ly[[i]] is the projected time series for subject i.
#' @param a_hat Length-n subject loading vector.
#' @param ind_vec Integer vector mapping each sample to its subject (1..n).
#' @param Kmat Kernel matrix between all observed time points.
#' @param Kmat_output Kernel matrix between resolution grid and observed points.
#' @param smooth RKHS penalty weight.
#' @noRd
freg_rkhs <- function(Ly, a_hat, ind_vec, Kmat, Kmat_output, smooth = 1e-8) {
  K <- Kmat
  for (i in 1:Ly) {
    K[ind_vec == i, ] <- K[ind_vec == i, ] * a_hat[i]^2
  }
  cvec <- unlist(Ly)
  beta <- solve(K + smooth * diag(ncol(K))) %*% cvec
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
  xy   <- abs(x %*% t(rep(1, length(y))) - rep(1, length(x)) %*% t(y))
  k4_xy <- 1/24 * ((xy - 0.5)^4 - 0.5 * (xy - 0.5)^2 + 7/240)
  return(k1_x %*% t(k1_y) + k2_x %*% t(k2_y) - k4_xy + 1)
}



# ------- (3) Initialization function ----------

#' Initialize feature loading's vector using SVD of mode-2-matricized tensor
#'
#' @param datlists Length M named list of length n lists of matrices.
#' @param m modality
#' @param p_m number of features for modality m
#' @returns list of M length p_m vectors
#'
#' @noRd No user-side documentation
init_b <- function(datlists, m, b_hat, p) {
  data_unfold <- NULL
    for (i in 1:n) {
      data_unfold = cbind(data_unfold, datlists[[m]][[i]][2:(p_m+1),])
    }
  b.intitials <- svd(data_unfold, nu=r, nv=r)$u
  b_hat <- b.intitials[,1]
  return(b_hat)
}



# ------- (4) Updating functions ----------

update_zeta <- function() {  # updates modality-specific time loadings
  # Kernel ridge regression code (RKHS penalty term)
  Ly <- list()
  for (i in 1:n){
    Ly <- c(Ly, list(a_hat[i]*as.numeric(b_hat%*%datlist[[i]][2:(p+1),])))
  }
  phi_hat <- freg_rkhs(Ly, a_hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
  phi_hat <- phi_hat / sqrt(sum(phi_hat^2))

  # Normalize
  zeta_hat <- zeta_hat / sqrt(sum(zeta_hat^2))
  return(zeta_hat)
}

update_a <- function() {  # updates cross-modality shared subject loading
  # update b
  a_tilde <- rep(0,n)
  for (i in 1:n){
    t_temp <- tipos[[i]]
    a_tilde[i] <- b_hat %*% datlist[[i]][2:(p+1),t_temp] %*% phi_hat[ti[[i]][t_temp]]
    a_tilde[i] <- a_tilde[i] / sum((phi_hat[ti[[i]][t_temp]])^2)
  }

  # Normalize and calculate dif
  a.new <- a_tilde / sqrt(sum(a_tilde^2))
  dif <- sum((a_hat - a.new)^2)
  a_hat <- a.new
}

update_b <- function() {  # updates feature loadings
  temp_num <- matrix(0,p,n)
  temp_denom <- rep(0,n)
  for (i in 1:n){
    t_temp <- tipos[[i]]
    temp_num[,i] <- datlist[[i]][2:(p+1),t_temp] %*% phi_hat[ti[[i]][t_temp]]
    temp_denom[i] <-sum((phi_hat[ti[[i]][t_temp]])^2)
  }
  b_tilde <- as.numeric(temp_num%*%a_hat) / as.numeric(temp_denom%*%(a_hat^2))
  b.new <- b_tilde / sqrt(sum(b_tilde^2))
  dif <- max(dif, sum((b_hat - b.new)^2))
  b_hat <- b.new


}

update_datlists <- function() {
}

# ------- (5) Post algorithm functions ----------

calc_lambda <- function() {
  x <- NULL
  for (i in 1:n){
    t_temp <- ti[[i]]
    t_temp <- t_temp[t_temp>0]
    x <- c(x,as.vector(t(a_hat[i]*b_hat%o%phi_hat[t_temp])))
  }
  X <- cbind(X, x)
  lm_fit <- lm(y~x-1)
  lambda <- as.numeric(lm_fit$coefficients)
  A[,s] <- a_hat
  B[,s] <- b_hat
  Phi[,s] <- t(phi_hat)
  Lambda[s] <- lambda
  Rsq[s] <- summary(lm_fit)$r.squared
  accumRsq[s] <- summary(lm(y0~X-1))$r.squared
}

revise_signs <- function() {

}


