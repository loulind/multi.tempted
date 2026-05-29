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

  # for (m in 1:M) {
  #   p[m] <- sapply(datlists[[m]], nrow) - 1  # number of features per modality
  #   if (!(length(unique(p)) == 1)) {
  #     stop(paste0("Modality '", names(datlists)[m], "' has inconsistent feature counts across subjects"))
  #   }
  #   B[[m]] <- matrix(0, p[m], r)
  #
  #   # Calculate range for each subject in each modality
  #   timestamps_all <- vector(mode = "list", length = M)
  #   interval <- vector(mode = "list", length = M)
  #   timestamps_all[[m]] <- do.call(c,lapply(datlists[[m]], FUN=function(u){u[1,]}))  # extracts time point row into vector
  #   timestamps_all[[m]] <- sort(unique(timestamps_all[[m]]))  #
  #
  #   if (is.null(interval)){ # initializes interval as 1st timepoint to last if not specified
  #     interval[[m]] <- c(timestamps_all[[m]][1], timestamps_all[[m]][length(timestamps_all[[m]])])
  #   }
  #
  #   # Rescale time to 0-1
  #   input_time_range <- c(timestamps_all[[m]][1], timestamps_all[[m]][length(timestamps_all[[m]])])
  #   for (i in 1:n) {
  #     datlists[[m]][[i]][1,] <- (datlists[[m]][[i]][1,] - input_time_range[[m]][1]) / (input_time_range[[m]][2] - input_time_range[[m]][1])
  #   }
  #   interval[[m]] <- (interval[[m]] - input_time_range[[m]][1]) / (input_time_range[[m]][2] - input_time_range[[m]][1])
  #
  #   # Binning continuous time to a discrete grid (based on interval and resolution)
  #   for (i in 1:n) {
  #     temp <- 1 + round((resolution-1) * (datlists[[m]][[i]][1,] - interval[[m]][1]) / (interval[[m]][2] - interval[[m]][1]))
  #     temp[which(temp<=0 | temp>resolution)] <- 0
  #     ti[[m]][[i]] <- temp  # subject i's timepoint
  #   }
  #
  #   # Flattening out feature data into long vector
  #   for (i in 1:n){
  #     keep <- ti[[m]]][[i]]>0
  #     tipos[[m]][[i]] <- keep
  #     y0[[m]] <- c(y0[[m]], as.vector(t(datlists[[m]][[i]][2:(p+1),keep])))
  #   }
  #
  #   # creates long vector of timepoints for calculating bernoulli kernel from all timepoints
  #   Lt <- list()
  #   ind_vec <- NULL
  #   for (i in 1:n){
  #     Lt <- c(Lt, list(datlist[[i]][1,]))
  #     ind_vec <- c(ind_vec, rep(i,length(Lt[[i]])))
  #   }
  #
  #   tm <- unlist(Lt)
  #   Kmat[[m]] <- bernoulli_kernel(tm, tm)
  #   Kmat_output[[m]] <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)
  # }

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

    # STEP 3: Store loadings and estimate modality-specific scales (lambda)
    A[, l] <- a_hat # record estimated a for component l
    x_comp <- NULL  # rank-1 reconstruction vectorized across all modalities
    for (m in 1:M) {
      B[[m]][, l] <- b_hats[[m]]  # record estimated b for component l
      Zeta[[m]][, l] <- zeta_hats[[m]]  # record estimated zeta for component l

      lm_result <- compute_lambda(y_resid[[m]], datlists[[m]], p[m],
                                     a_hat, b_hats[[m]], zeta_hats[[m]],
                                     prep[[m]]$tipos, prep[[m]]$ti, n)
      Lambda[m, l] <- lm_result$lambda # update lambda
      x_comp <- c(x_comp, lm_result$x_m)
    }

    # R-squared for this component and accumulated total (pooled across modalities)
    Rsq[l]      <- compute_rsq(unlist(y_resid), x_comp)
    X_accum     <- cbind(X_accum, x_comp)
    accumRsq[l] <- compute_rsq(y0_all, X_accum)

    # STEP 4: Subtract l'th component from datlists
    for (m in 1:M) {
      datlists[[m]] <- update_datlist(datlists[[m]], p[m], a_hat, b_hats[[m]],
                                       zeta_hats[[m]], Lambda[m, l],
                                       prep[[m]]$tipos, prep[[m]]$ti, n)
    }
  }

  # for (r in 1:length(Lambda[[m]])){
  #   for (m in 1:M) {
  #     # revise the sign of Lambdas
  #     if (Lambda[[m]][r]<0){
  #       Lambda[[m]][r] <- -Lambda[[m]][r]
  #       A[,r] <- -A[,r]
  #     }
  #
  #     # revise the signs to make sure summation of zeta is nonnegative
  #     sgn.zeta <- sign(colSums(Zeta[[m]]))
  #     sgn.zeta[sgn.zeta==0] <- 1
  #     for (r in 1:ncol(Phi)){
  #       Zeta[[m]][,r] <- sgn.phi[r]*Zeta[[m]][,r]
  #       A[,r] <- sgn.zeta[r]*A[,r]
  #     }
  #
  #     # revise the signs to make sure summation of B is nonnegative
  #     sgn.B <- sign(colSums(B[[m]]))
  #     sgn.B[sgn.B==0] <- 1
  #     for (r in 1:ncol(Phi)){
  #       B[[m]][,r] <- sgn.B[r]*B[[m]][,r]
  #       A[,r] <- sgn.B[r]*A[,r]
  #     }
  #   }
  # }

  time_return <- seq(interval[1],interval[2],length.out = resolution)
  time_return <- time_return * (input_time_range[2] - input_time_range[1]) + input_time_range[1]
  results <- list("A_hat" = A, "B_hat" = B,
                  "Zeta_hat" = Zeta, "time_Phi" = time_return,
                  "Lambda" = Lambda, "r_square" = Rsq, "accum_r_square" = accumRsq)
  return(results)
}






# HELPER FUNCTIONS

init_time <- function(datlists) {

}


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

#' Functional Regression with RKHS penalty term
#'
#' @param Ly
#' @param a_hat
#' @param ind_vec
#' @param Kmat
#' @param Kmat_output
#' @param smooth
#'
#' @returns
#'
#' @noRd No user-side documentation
freg_rkhs <- function(Ly, a_hat, ind_vec, Kmat, Kmat_output, smooth=1e-8){
  A <- Kmat
  for (i in 1:length(Ly)){
    A[ind_vec==i,] <- A[ind_vec==i,]*a_hat[i]^2
  }
  cvec <- unlist(Ly)

  A_temp <- A + smooth*diag(ncol(A))
  beta <- solve(A_temp)%*%cvec

  zeta_est <- Kmat_output %*% beta
  return(zeta_est)
}

bernoulli_kernel <- function(x, y){
  k1_x <- x-0.5
  k1_y <- y-0.5
  k2_x <- 0.5*(k1_x^2-1/12)
  k2_y <- 0.5*(k1_y^2-1/12)
  xy <- abs(x %*% t(rep(1,length(y))) - rep(1,length(x)) %*% t(y))
  k4_xy <- 1/24 * ((xy-0.5)^4 - 0.5*(xy-0.5)^2 + 7/240)
  kern_xy <- k1_x %*% t(k1_y) + k2_x %*% t(k2_y) - k4_xy + 1
  return(kern_xy)
}

revise_signs <- function() {

}
