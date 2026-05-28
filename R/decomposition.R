#' @title Decomposition of temporal tensors
#' @description
#' CP-type decomposition of multiple 3d tensors.
#' Each tensor represents a different modality (subject x feature x time).
#' Performed after formatting data!
#' @param datlists A length M named list of length n lists of matrices.
#' Each named list element represents a modality.
#' Each matrix represents a subject. Columns represent sample number.
#' The first rows represent sampling time points.
#' Row 2 through row (number features in modality m) + 1 represent feature values
#' @param r Number components to decompose the M 3d tensors
#' (ie, rank of CP-type decomposition). Default is r=3.
#' @param smooth Smoothing parameter for RKHS norm.
#' Larger ==> smoother temporal loading functions. Default is 1e-8.
#' To adjust, check the smoothness of the estimated temporal loading function plot.
#' @param interval The range of time points to run the decomposition for.
#' Default is set to be the range of all observed time points.
#' User can set it to be a shorter interval than the observed range.
#' @param resolution Number of time points to evaluate the value of the temporal loading function.
#' Default is set to 101. It does not affect the subject or feature loadings.
#' @param maxiter Maximum number of iteration. Default is 20.
#' @param epsilon Convergence criteria for difference between iterations. Default is 1e-4.
#' @return The estimations of the loadings for each modality.
#' @export
#' @examples
multi_tempted_decomp <- function(datlists, r=3, smooth=1e-8, interval=NULL,
                                 resolution = 101, maxiter=20, epsilon=1e-4) {
  if (!(length(unique(lengths(datlists))) == 1)) {
    stop("Unequal number of subjects across modalities")
  }

  # Initialize intermediate variables
  M <- length(datlists)  # number modalities
  n <- length(datlists[[1]])  # number subjects
  A <- matrix(0, nrow = n, ncol = r) # subject loading components
  B <- vector(mode = "list", length = M) # feature loading components
  p <- vector(mode = "numeric", length = M) # number features per modality

  ti <- lapply(1:M, function(m) lapply(1:n, function(x) vector()))  # M lists containing subject timepoint indices
  tipos <- lapply(1:M, function(m) lapply(1:n, function(x) vector()))  # M lists of whether timepoints falls in range
  Kmat <- list() # list of matrices to calc Bernoulli kernel between all observed time points
  Kmat_output <- list() # list of matrices to calc Bernoulli kernel between resolution grid and actual observed time points

  Lambda <- replicate(M, numeric(r), simplify = FALSE)  # modality-specific scalings
  X <- NULL  # design matrix
  y0 <- vector(mode = "list", length=M)  #flattened feature data
  Rsq <- accumRsq <- rep(0, r)

  for (m in 1:M) {
    p[m] <- sapply(datlists[[m]], nrow) - 1  # number of features per modality
    if (!(length(unique(p)) == 1)) {
      stop(paste0("Modality '", names(datlists)[m], "' has inconsistent feature counts across subjects"))
    }
    B[[m]] <- matrix(0, p[m], r)

    # Calculate range for each subject in each modality
    timestamps_all <- vector(mode = "list", length = M)
    interval <- vector(mode = "list", length = M)
    timestamps_all[[m]] <- do.call(c,lapply(datlists[[m]], FUN=function(u){u[1,]}))  # extracts time point row into vector
    timestamps_all[[m]] <- sort(unique(timestamps_all[[m]]))  #

    if (is.null(interval)){ # initializes interval as 1st timepoint to last if not specified
      interval[[m]] <- c(timestamps_all[[m]][1], timestamps_all[[m]][length(timestamps_all[[m]])])
    }

    # Rescale time to 0-1
    input_time_range <- c(timestamps_all[[m]][1], timestamps_all[[m]][length(timestamps_all[[m]])])
    for (i in 1:n) {
      datlists[[m]][[i]][1,] <- (datlists[[m]][[i]][1,] - input_time_range[[m]][1]) / (input_time_range[[m]][2] - input_time_range[[m]][1])
    }
    interval[[m]] <- (interval[[m]] - input_time_range[[m]][1]) / (input_time_range[[m]][2] - input_time_range[[m]][1])

    # Binning continuous time to a discrete grid (based on interval and resolution)
    for (i in 1:n) {
      temp <- 1 + round((resolution-1) * (datlists[[m]][[i]][1,] - interval[[m]][1]) / (interval[[m]][2] - interval[[m]][1]))
      temp[which(temp<=0 | temp>resolution)] <- 0
      ti[[m]][[i]] <- temp  # subject i's timepoint
    }

    # Flattening out feature data into long vector
    for (i in 1:n){
      keep <- ti[[m]]][[i]]>0
      tipos[[m]][[i]] <- keep
      y0[[m]] <- c(y0[[m]], as.vector(t(datlists[[m]][[i]][2:(p+1),keep])))
    }

    # creates long vector of timepoints for calculating bernoulli kernel from all timepoints
    Lt <- list()
    ind_vec <- NULL
    for (i in 1:n){
      Lt <- c(Lt, list(datlist[[i]][1,]))
      ind_vec <- c(ind_vec, rep(i,length(Lt[[i]])))
    }

    tm <- unlist(Lt)
    Kmat[[m]] <- bernoulli_kernel(tm, tm)
    Kmat_output[[m]] <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)
  }

  # Calculate each component and remove contribution from feature values
  for (l in 1:r) {
    message(sprintf("Calculating Component %d", l))

    # SEQUENTIAL ESTIMATION ALGORITHM
    # STEP 1: Initialize subject and feature loadings per modality
    ## (i) Subject loadings, a, initially set to equal contribution
    a_hat <- rep(1/sqrt(n), n)

    for (m in 1:M) {
      message(sprintf("...modality '%s'", names(datlists)[m]))

      ## (ii) Feature loadings, b, init'd as list of matrices of leftmost singular vectors of SVD
      b_hat <- init_b(datlists, m, p[m])

      # STEP 2: Update a, b, and zeta until max iterations reached or it converges
      t <- 0
      dif <- 1
      while(t<=maxiter & dif>epsilon){
        ## (i) Update time-varying function, zeta
        zeta_hat <- update_zeta(Ly, a_hat, ind_vec, Kmat, Kmat_output, smooth=smooth)

        ## (ii) Update subject loadings
        a_hat <- update_a()

        ## (iii) Update feature loadings
        b_hat <- update_b()

        t <- t+1
      }

      # STEP 3: Remove contribution of current component; repeat steps 1-2 for all r components
      update_datlists()
    }

  # STEP 4: Estimate modality-specific scales, lambda
    calc_lambda()
  }


  for (r in 1:length(Lambda[[m]])){
    for (m in 1:M) {
      # revise the sign of Lambdas
      if (Lambda[[m]][r]<0){
        Lambda[[m]][r] <- -Lambda[[m]][r]
        A[,r] <- -A[,r]
      }

      # revise the signs to make sure summation of zeta is nonnegative
      sgn.zeta <- sign(colSums(Zeta[[m]]))
      sgn.zeta[sgn.zeta==0] <- 1
      for (r in 1:ncol(Phi)){
        Zeta[[m]][,r] <- sgn.phi[r]*Zeta[[m]][,r]
        A[,r] <- sgn.zeta[r]*A[,r]
      }

      # revise the signs to make sure summation of B is nonnegative
      sgn.B <- sign(colSums(B[[m]]))
      sgn.B[sgn.B==0] <- 1
      for (r in 1:ncol(Phi)){
        B[[m]][,r] <- sgn.B[r]*B[[m]][,r]
        A[,r] <- sgn.B[r]*A[,r]
      }
    }
  }

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
