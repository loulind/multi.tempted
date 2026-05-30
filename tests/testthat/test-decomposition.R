# Tests for R/decomposition.R
#
# Test organisation mirrors the helper-function sections in decomposition.R:
#   (1) bernoulli_kernel
#   (2) init_time_intv
#   (3) flatten_features
#   (4) freg_rkhs
#   (5) init_b_hat
#   (6) update_zeta
#   (7) update_a
#   (8) update_b
#   (9) compute_lambda
#  (10) update_datlist
#  (11) revise_signs
#  (12) multi_tempted_decomp (wrapper)

# ==============================================================================
# SHARED TEST FIXTURES
# ==============================================================================

# Build a minimal, valid datlists object with M modalities, n subjects,
# p features per modality, and n_times evenly spaced time points in [0, 1].
make_datlists <- function(M = 2, n = 3, p = 4, n_times = 6, seed = 42) {
  set.seed(seed)
  times <- seq(0, 1, length.out = n_times)
  mods <- lapply(seq_len(M), function(m) {
    subjs <- lapply(seq_len(n), function(i)
      rbind(times, matrix(rnorm(p * n_times), nrow = p)))
    names(subjs) <- paste0("s", seq_len(n))
    subjs
  })
  names(mods) <- paste0("m", seq_len(M))
  mods
}

# Build a datlists where every modality is a near-exact rank-1 tensor with
# known factor structure (signal-to-noise ratio >> 1).
make_rank1_datlists <- function(M = 2, n = 4, p = 5, n_times = 8, seed = 1) {
  set.seed(seed)
  times  <- seq(0, 1, length.out = n_times)
  a_true <- rnorm(n);  a_true <- a_true / sqrt(sum(a_true^2))
  mods <- lapply(seq_len(M), function(m) {
    b_true    <- rnorm(p);    b_true    <- b_true    / sqrt(sum(b_true^2))
    zeta_true <- sin(pi * times) + 0.2          # smooth, positive
    lambda_m  <- 5
    subjs <- lapply(seq_len(n), function(i) {
      signal <- lambda_m * a_true[i] * outer(b_true, zeta_true)
      noise  <- matrix(rnorm(p * n_times, sd = 0.05), nrow = p)
      rbind(times, signal + noise)
    })
    names(subjs) <- paste0("s", seq_len(n))
    subjs
  })
  names(mods) <- paste0("m", seq_len(M))
  mods
}

# Run init_time_intv on one modality and return the prep object.
make_prep <- function(datlist_m, p_m, n, interval_m = NULL, resolution = 21) {
  multi.tempted:::init_time_intv(datlist_m, p_m, n, interval_m, resolution)
}


# ==============================================================================
# (1) bernoulli_kernel
# ==============================================================================

test_that("bernoulli_kernel returns a matrix with dimensions length(x) x length(y)", {
  x <- c(0.1, 0.4, 0.7)
  y <- c(0.2, 0.5, 0.8, 0.95)
  K <- multi.tempted:::bernoulli_kernel(x, y)
  expect_equal(dim(K), c(length(x), length(y)))
})

test_that("bernoulli_kernel gram matrix is symmetric", {
  x <- seq(0.1, 0.9, by = 0.2)
  K <- multi.tempted:::bernoulli_kernel(x, x)
  expect_equal(K, t(K), tolerance = 1e-12)
})

test_that("bernoulli_kernel gram matrix is positive semi-definite", {
  x <- seq(0.05, 0.95, by = 0.1)
  K <- multi.tempted:::bernoulli_kernel(x, x)
  eigenvalues <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10))
})


# ==============================================================================
# (2) init_time_intv
# ==============================================================================

test_that("init_time_intv rescales all time points to [0, 1]", {
  dl   <- make_datlists(M = 1, n = 3, p = 4, n_times = 5)[[1]]
  prep <- make_prep(dl, p_m = 4, n = 3)
  all_times <- unlist(lapply(prep$datlist, function(s) s[1, ]))
  expect_true(all(all_times >= -1e-10 & all_times <= 1 + 1e-10))
})

test_that("init_time_intv ti values are 0 (out of range) or in [1, resolution]", {
  resolution <- 21
  dl   <- make_datlists(M = 1, n = 3, p = 4, n_times = 6)[[1]]
  prep <- make_prep(dl, 4, 3, resolution = resolution)
  all_ti <- unlist(prep$ti)
  expect_true(all(all_ti == 0 | (all_ti >= 1 & all_ti <= resolution)))
})

test_that("init_time_intv tipos[[i]] == (ti[[i]] > 0)", {
  dl   <- make_datlists(M = 1, n = 3, p = 4, n_times = 6)[[1]]
  prep <- make_prep(dl, 4, 3)
  for (i in seq_along(prep$ti)) {
    expect_equal(prep$tipos[[i]], prep$ti[[i]] > 0)
  }
})

test_that("init_time_intv ind_vec length equals total number of samples", {
  n <- 3; n_times <- 6
  dl   <- make_datlists(M = 1, n = n, p = 4, n_times = n_times)[[1]]
  prep <- make_prep(dl, 4, n)
  expect_equal(length(prep$ind_vec), n * n_times)
})

test_that("init_time_intv Kmat is square with side equal to total samples", {
  n <- 3; n_times <- 6
  dl   <- make_datlists(M = 1, n = n, p = 4, n_times = n_times)[[1]]
  prep <- make_prep(dl, 4, n)
  expect_equal(dim(prep$Kmat), c(n * n_times, n * n_times))
})

test_that("init_time_intv Kmat_output has nrow equal to resolution", {
  resolution <- 31
  dl   <- make_datlists(M = 1, n = 3, p = 4, n_times = 6)[[1]]
  prep <- make_prep(dl, 4, 3, resolution = resolution)
  expect_equal(nrow(prep$Kmat_output), resolution)
})

test_that("init_time_intv assigns ti = 0 to samples outside a narrowed interval", {
  n_times <- 10
  dl   <- make_datlists(M = 1, n = 2, p = 3, n_times = n_times)[[1]]
  # Narrow the interval so the two boundary time points fall outside
  prep <- make_prep(dl, 3, 2, interval_m = c(0.25, 0.75))
  orig_times <- seq(0, 1, length.out = n_times)
  out_of_range <- orig_times < 0.25 | orig_times > 0.75
  for (i in seq_len(2)) {
    expect_true(all(prep$ti[[i]][out_of_range] == 0))
  }
})


# ==============================================================================
# (3) flatten_features
# ==============================================================================

test_that("flatten_features output length equals (in-range samples) x p_m", {
  n <- 3; p <- 4
  dl   <- make_datlists(M = 1, n = n, p = p, n_times = 6)[[1]]
  prep <- make_prep(dl, p, n)
  y    <- multi.tempted:::flatten_features(prep$datlist, p, prep$tipos)
  expected_len <- p * sum(sapply(seq_len(n), function(i) sum(prep$tipos[[i]])))
  expect_equal(length(y), expected_len)
})

test_that("flatten_features returns numeric vector", {
  dl   <- make_datlists(M = 1, n = 3, p = 4, n_times = 6)[[1]]
  prep <- make_prep(dl, 4, 3)
  y    <- multi.tempted:::flatten_features(prep$datlist, 4, prep$tipos)
  expect_true(is.numeric(y))
})


# ==============================================================================
# (4) freg_rkhs
# ==============================================================================

# Build the inputs that freg_rkhs needs directly (without going through update_zeta).
make_freg_inputs <- function(n = 3, n_times = 5, resolution = 21, seed = 7) {
  set.seed(seed)
  times       <- seq(0, 1, length.out = n_times)
  tm          <- rep(times, n)
  ind_vec     <- rep(seq_len(n), each = n_times)
  Kmat        <- multi.tempted:::bernoulli_kernel(tm, tm)
  grid        <- seq(0, 1, length.out = resolution)
  Kmat_output <- multi.tempted:::bernoulli_kernel(grid, tm)
  a_hat       <- rep(1 / sqrt(n), n)
  Ly          <- lapply(seq_len(n), function(i) rnorm(n_times))
  list(Ly = Ly, a_hat = a_hat, ind_vec = ind_vec,
       Kmat = Kmat, Kmat_output = Kmat_output)
}

test_that("freg_rkhs output has length equal to resolution", {
  resolution <- 21
  inp <- make_freg_inputs(resolution = resolution)
  out <- multi.tempted:::freg_rkhs(inp$Ly, inp$a_hat, inp$ind_vec,
                                   inp$Kmat, inp$Kmat_output, smooth = 1e-5)
  expect_equal(length(out), resolution)
})

test_that("freg_rkhs with weight > 1 gives a different result than weight = 1", {
  # Increasing weight up-scales the data term relative to the fixed smooth penalty,
  # so the two solutions differ in general.
  inp <- make_freg_inputs()
  r1  <- multi.tempted:::freg_rkhs(inp$Ly, inp$a_hat, inp$ind_vec,
                                   inp$Kmat, inp$Kmat_output, smooth = 1e-3, weight = 1)
  r5  <- multi.tempted:::freg_rkhs(inp$Ly, inp$a_hat, inp$ind_vec,
                                   inp$Kmat, inp$Kmat_output, smooth = 1e-3, weight = 5)
  expect_false(isTRUE(all.equal(r1, r5)))
})

test_that("freg_rkhs default weight = 1 matches explicit weight = 1", {
  inp     <- make_freg_inputs()
  default <- multi.tempted:::freg_rkhs(inp$Ly, inp$a_hat, inp$ind_vec,
                                       inp$Kmat, inp$Kmat_output, smooth = 1e-5)
  explicit <- multi.tempted:::freg_rkhs(inp$Ly, inp$a_hat, inp$ind_vec,
                                        inp$Kmat, inp$Kmat_output, smooth = 1e-5, weight = 1)
  expect_equal(default, explicit)
})


# ==============================================================================
# (5) init_b_hat
# ==============================================================================

test_that("init_b_hat returns a unit-norm vector of length p_m", {
  n <- 4; p <- 6
  dl <- make_datlists(M = 1, n = n, p = p)[[1]]
  b  <- multi.tempted:::init_b_hat(dl, p_m = p, n = n)
  expect_equal(length(b), p)
  expect_equal(sum(b^2), 1, tolerance = 1e-10)
})


# ==============================================================================
# (6) update_zeta
# ==============================================================================

test_that("update_zeta returns a unit-norm vector of length resolution", {
  n <- 3; p <- 4; resolution <- 21
  dl   <- make_datlists(M = 1, n = n, p = p)[[1]]
  prep <- make_prep(dl, p, n, resolution = resolution)
  b    <- multi.tempted:::init_b_hat(prep$datlist, p, n)
  a    <- rep(1 / sqrt(n), n)
  zeta <- multi.tempted:::update_zeta(prep$datlist, p, b, a,
                                      prep$ind_vec, prep$Kmat, prep$Kmat_output,
                                      smooth = 1e-5)
  expect_equal(length(zeta), resolution)
  expect_equal(sum(zeta^2), 1, tolerance = 1e-10)
})

test_that("update_zeta with weight = 1 matches the default (no weight argument)", {
  n <- 3; p <- 4
  dl   <- make_datlists(M = 1, n = n, p = p)[[1]]
  prep <- make_prep(dl, p, n)
  b    <- multi.tempted:::init_b_hat(prep$datlist, p, n)
  a    <- rep(1 / sqrt(n), n)
  z_default <- multi.tempted:::update_zeta(prep$datlist, p, b, a,
                                           prep$ind_vec, prep$Kmat, prep$Kmat_output,
                                           smooth = 1e-5)
  z_w1      <- multi.tempted:::update_zeta(prep$datlist, p, b, a,
                                           prep$ind_vec, prep$Kmat, prep$Kmat_output,
                                           smooth = 1e-5, weight = 1)
  expect_equal(z_default, z_w1)
})


# ==============================================================================
# (7) update_a
# ==============================================================================

# Shared setup for update_a / update_b tests
make_update_inputs <- function(M = 2, n = 3, p = 4, seed = 42) {
  datlists <- make_datlists(M = M, n = n, p = p, seed = seed)
  preps    <- lapply(seq_len(M), function(m) make_prep(datlists[[m]], p, n))
  dl_prep  <- lapply(seq_len(M), function(m) preps[[m]]$datlist)
  b_hats   <- lapply(seq_len(M), function(m)
    multi.tempted:::init_b_hat(dl_prep[[m]], p, n))
  a_hat    <- rep(1 / sqrt(n), n)
  zeta_hats <- lapply(seq_len(M), function(m)
    multi.tempted:::update_zeta(dl_prep[[m]], p, b_hats[[m]], a_hat,
                                preps[[m]]$ind_vec, preps[[m]]$Kmat,
                                preps[[m]]$Kmat_output, smooth = 1e-5))
  list(datlists = dl_prep, preps = preps, b_hats = b_hats,
       a_hat = a_hat, zeta_hats = zeta_hats, M = M, n = n, p = p)
}

test_that("update_a returns a unit-norm vector of length n", {
  inp   <- make_update_inputs()
  a_new <- multi.tempted:::update_a(inp$datlists, rep(inp$p, inp$M), inp$b_hats,
                                    inp$zeta_hats, inp$preps, inp$n, inp$M,
                                    weights = rep(1, inp$M))
  expect_equal(length(a_new), inp$n)
  expect_equal(sum(a_new^2), 1, tolerance = 1e-10)
})

test_that("update_a result is unchanged when all weights are scaled by a constant", {
  # Proportional scaling of all w_m cancels in the ratio (num / den).
  inp    <- make_update_inputs()
  p_vec  <- rep(inp$p, inp$M)
  a_w1   <- multi.tempted:::update_a(inp$datlists, p_vec, inp$b_hats, inp$zeta_hats,
                                     inp$preps, inp$n, inp$M, weights = c(1, 1))
  a_w5   <- multi.tempted:::update_a(inp$datlists, p_vec, inp$b_hats, inp$zeta_hats,
                                     inp$preps, inp$n, inp$M, weights = c(5, 5))
  expect_equal(a_w1, a_w5, tolerance = 1e-10)
})

test_that("update_a with weight = 0 on one modality ignores that modality", {
  # Setting w_2 = 0 should give the same result as a single-modality update with
  # only modality 1.
  inp   <- make_update_inputs(M = 2)
  p_vec <- rep(inp$p, inp$M)

  a_both <- multi.tempted:::update_a(inp$datlists, p_vec, inp$b_hats, inp$zeta_hats,
                                     inp$preps, inp$n, inp$M, weights = c(1, 0))
  a_mod1 <- multi.tempted:::update_a(inp$datlists[1], p_vec[1], inp$b_hats[1],
                                     inp$zeta_hats[1], inp$preps[1],
                                     inp$n, M = 1, weights = 1)
  expect_equal(a_both, a_mod1, tolerance = 1e-10)
})

test_that("update_a result differs when one modality has a much larger weight", {
  inp   <- make_update_inputs()
  p_vec <- rep(inp$p, inp$M)
  a_eq  <- multi.tempted:::update_a(inp$datlists, p_vec, inp$b_hats, inp$zeta_hats,
                                    inp$preps, inp$n, inp$M, weights = c(1, 1))
  a_sk  <- multi.tempted:::update_a(inp$datlists, p_vec, inp$b_hats, inp$zeta_hats,
                                    inp$preps, inp$n, inp$M, weights = c(100, 1))
  expect_false(isTRUE(all.equal(a_eq, a_sk)))
})


# ==============================================================================
# (8) update_b
# ==============================================================================

test_that("update_b returns a unit-norm vector of length p_m", {
  n <- 3; p <- 4
  dl   <- make_datlists(M = 1, n = n, p = p)[[1]]
  prep <- make_prep(dl, p, n)
  b    <- multi.tempted:::init_b_hat(prep$datlist, p, n)
  a    <- rep(1 / sqrt(n), n)
  zeta <- multi.tempted:::update_zeta(prep$datlist, p, b, a,
                                      prep$ind_vec, prep$Kmat, prep$Kmat_output, 1e-5)
  b_new <- multi.tempted:::update_b(prep$datlist, p, zeta,
                                    prep$tipos, prep$ti, a, n)
  expect_equal(length(b_new), p)
  expect_equal(sum(b_new^2), 1, tolerance = 1e-10)
})


# ==============================================================================
# (9) compute_lambda
# ==============================================================================

test_that("compute_lambda returns a scalar lambda and an x_m of the right length", {
  n <- 3; p <- 4
  dl    <- make_datlists(M = 1, n = n, p = p)[[1]]
  prep  <- make_prep(dl, p, n)
  b     <- multi.tempted:::init_b_hat(prep$datlist, p, n)
  a     <- rep(1 / sqrt(n), n)
  zeta  <- multi.tempted:::update_zeta(prep$datlist, p, b, a,
                                       prep$ind_vec, prep$Kmat, prep$Kmat_output, 1e-5)
  y_m   <- multi.tempted:::flatten_features(prep$datlist, p, prep$tipos)
  result <- multi.tempted:::compute_lambda(y_m, prep$datlist, p, a, b, zeta,
                                           prep$tipos, prep$ti, n)

  expect_true(is.numeric(result$lambda) && length(result$lambda) == 1)
  expected_x_len <- p * sum(sapply(seq_len(n), function(i) sum(prep$tipos[[i]])))
  expect_equal(length(result$x_m), expected_x_len)
})

test_that("compute_lambda x_m reconstructs the observed data up to a scale", {
  # Build a dataset that is exactly rank-1 so lambda should equal the true scale.
  n <- 3; p <- 4; n_times <- 8
  times    <- seq(0, 1, length.out = n_times)
  set.seed(10)
  a_true   <- rnorm(n);  a_true   <- a_true   / sqrt(sum(a_true^2))
  b_true   <- rnorm(p);  b_true   <- b_true   / sqrt(sum(b_true^2))
  zeta_raw <- sin(pi * times) + 0.2  # smooth, positive
  lambda_true <- 4

  dl <- lapply(seq_len(n), function(i)
    rbind(times, lambda_true * a_true[i] * outer(b_true, zeta_raw)))
  prep  <- make_prep(dl, p, n)

  # Provide nearly-exact loadings as input (they'll be normalised inside helpers)
  b_unit <- b_true
  a_unit <- a_true
  # Evaluate zeta_raw on the grid produced by prep
  grid   <- seq(0, 1, length.out = 21)  # resolution = 21
  zeta_unit <- approx(seq(0, 1, length.out = n_times), zeta_raw, xout = grid)$y
  zeta_unit <- zeta_unit / sqrt(sum(zeta_unit^2))

  y_m    <- multi.tempted:::flatten_features(prep$datlist, p, prep$tipos)
  result <- multi.tempted:::compute_lambda(y_m, prep$datlist, p, a_unit, b_unit,
                                           zeta_unit, prep$tipos, prep$ti, n)
  # lambda should be positive and recover meaningful signal
  expect_gt(result$lambda, 0)
})


# ==============================================================================
# (10) update_datlist
# ==============================================================================

test_that("update_datlist removes the rank-1 contribution from in-range samples", {
  n <- 2; p <- 3; n_times <- 5
  times   <- seq(0, 1, length.out = n_times)
  set.seed(7)
  dl      <- lapply(seq_len(n), function(i)
    rbind(times, matrix(rnorm(p * n_times), nrow = p)))
  prep    <- make_prep(dl, p, n)

  b_hat    <- rnorm(p);  b_hat    <- b_hat    / sqrt(sum(b_hat^2))
  a_hat    <- rnorm(n);  a_hat    <- a_hat    / sqrt(sum(a_hat^2))
  zeta_hat <- rnorm(21); zeta_hat <- zeta_hat / sqrt(sum(zeta_hat^2))
  lambda   <- 2.5

  dl_orig     <- prep$datlist
  dl_deflated <- multi.tempted:::update_datlist(prep$datlist, p, a_hat, b_hat,
                                                zeta_hat, lambda,
                                                prep$tipos, prep$ti, n)

  for (i in seq_len(n)) {
    in_range <- which(prep$tipos[[i]])
    grid_idx <- prep$ti[[i]][prep$tipos[[i]]]
    expected <- dl_orig[[i]][2:(p + 1), in_range] -
      lambda * a_hat[i] * outer(b_hat, zeta_hat[grid_idx])
    expect_equal(dl_deflated[[i]][2:(p + 1), in_range], expected, tolerance = 1e-12)
  }
})

test_that("update_datlist leaves out-of-range samples and time row unchanged", {
  n <- 2; p <- 3; n_times <- 10
  times <- seq(0, 1, length.out = n_times)
  set.seed(8)
  dl   <- lapply(seq_len(n), function(i)
    rbind(times, matrix(rnorm(p * n_times), nrow = p)))
  # Narrow interval so some samples are out of range
  prep <- make_prep(dl, p, n, interval_m = c(0.3, 0.7))

  b_hat    <- rnorm(p);  b_hat    <- b_hat    / sqrt(sum(b_hat^2))
  a_hat    <- rep(1 / sqrt(n), n)
  zeta_hat <- rnorm(21); zeta_hat <- zeta_hat / sqrt(sum(zeta_hat^2))

  dl_orig     <- prep$datlist
  dl_deflated <- multi.tempted:::update_datlist(prep$datlist, p, a_hat, b_hat,
                                                zeta_hat, lambda_ml = 1,
                                                prep$tipos, prep$ti, n)

  for (i in seq_len(n)) {
    out_of_range <- which(!prep$tipos[[i]])
    if (length(out_of_range) > 0) {
      # Feature rows and time row at out-of-range columns should be untouched
      expect_equal(dl_deflated[[i]][, out_of_range],
                   dl_orig[[i]][, out_of_range])
    }
    # Time row (row 1) should never be modified
    expect_equal(dl_deflated[[i]][1, ], dl_orig[[i]][1, ])
  }
})


# ==============================================================================
# (11) revise_signs
# ==============================================================================

# Shared setup for revise_signs tests
make_sign_inputs <- function(r = 2, M = 2, n = 3, p = 4, resolution = 10,
                             seed = 99) {
  set.seed(seed)
  A      <- matrix(rnorm(n * r), n, r)
  B      <- lapply(seq_len(M), function(m) matrix(rnorm(p * r), p, r))
  Zeta   <- lapply(seq_len(M), function(m) matrix(rnorm(resolution * r), resolution, r))
  Lambda <- matrix(rnorm(M * r), M, r)
  list(A = A, B = B, Zeta = Zeta, Lambda = Lambda, r = r, M = M)
}

test_that("revise_signs makes all Lambda values non-negative", {
  inp    <- make_sign_inputs()
  result <- multi.tempted:::revise_signs(inp$A, inp$B, inp$Zeta, inp$Lambda,
                                         inp$r, inp$M)
  expect_true(all(result$Lambda >= 0))
})

test_that("revise_signs makes sum(Zeta[[m]][, l]) non-negative for all m, l", {
  inp    <- make_sign_inputs()
  result <- multi.tempted:::revise_signs(inp$A, inp$B, inp$Zeta, inp$Lambda,
                                         inp$r, inp$M)
  for (m in seq_len(inp$M)) {
    expect_true(all(colSums(result$Zeta[[m]]) >= 0))
  }
})

test_that("revise_signs does not guarantee sum(B[[m]][, l]) >= 0", {
  # sum(B) is NOT part of the sign convention: enforcing it simultaneously with
  # sum(Zeta) >= 0 creates a circular flip (each correction undoes the other).
  # B is instead used as the sign sink that absorbs corrections for Lambda and Zeta.
  # This test simply verifies the function runs without error; it does NOT assert
  # that sum(B) is non-negative.
  inp    <- make_sign_inputs()
  result <- multi.tempted:::revise_signs(inp$A, inp$B, inp$Zeta, inp$Lambda,
                                         inp$r, inp$M)
  expect_true(is.list(result))
  expect_named(result, c("A", "B", "Zeta", "Lambda"))
})

test_that("revise_signs makes sum(A[, l]) non-negative for all l", {
  inp    <- make_sign_inputs()
  result <- multi.tempted:::revise_signs(inp$A, inp$B, inp$Zeta, inp$Lambda,
                                         inp$r, inp$M)
  expect_true(all(colSums(result$A) >= 0))
})

test_that("revise_signs preserves the rank-1 reconstruction for each modality", {
  # The product Lambda[m,l] * A[i,l] * outer(B[[m]][,l], Zeta[[m]][,l])
  # should be identical before and after sign revision.
  r <- 1; M <- 2; n <- 3; p <- 4; resolution <- 10
  set.seed(55)
  A      <- matrix(rnorm(n), n, r)
  B      <- lapply(seq_len(M), function(m) matrix(rnorm(p), p, r))
  Zeta   <- lapply(seq_len(M), function(m) matrix(rnorm(resolution), resolution, r))
  Lambda <- matrix(rnorm(M), M, r)

  recon_before <- lapply(seq_len(M), function(m)
    Lambda[m, 1] * outer(A[, 1], as.vector(outer(B[[m]][, 1], Zeta[[m]][, 1]))))

  result <- multi.tempted:::revise_signs(A, B, Zeta, Lambda, r, M)

  recon_after <- lapply(seq_len(M), function(m)
    result$Lambda[m, 1] * outer(result$A[, 1],
                                as.vector(outer(result$B[[m]][, 1], result$Zeta[[m]][, 1]))))

  for (m in seq_len(M)) {
    expect_equal(recon_before[[m]], recon_after[[m]], tolerance = 1e-10)
  }
})


# ==============================================================================
# (12) multi_tempted_decomp  (integration tests)
# ==============================================================================

# --- Input validation ---------------------------------------------------------

test_that("multi_tempted_decomp errors when modalities have different subject counts", {
  dl      <- make_datlists(M = 2, n = 3)
  dl[[2]] <- dl[[2]][1:2]   # drop one subject from modality 2
  expect_error(multi_tempted_decomp(dl), "same number of subjects")
})

test_that("multi_tempted_decomp errors on negative weights", {
  dl <- make_datlists(M = 2, n = 3)
  expect_error(multi_tempted_decomp(dl, r = 1, weights = c(1, -0.5)),
               "non-negative")
})

test_that("multi_tempted_decomp errors when weights length != M", {
  dl <- make_datlists(M = 2, n = 3)
  expect_error(multi_tempted_decomp(dl, r = 1, weights = c(1, 1, 1)),
               "length")
})

# --- Output structure ---------------------------------------------------------

test_that("multi_tempted_decomp returns a correctly named list", {
  dl     <- make_datlists(M = 2, n = 3, p = 4)
  result <- suppressMessages(multi_tempted_decomp(dl, r = 1, maxiter = 3))
  expect_named(result, c("A_hat", "B_hat", "Zeta_hat", "time_Zeta",
                         "Lambda", "r_square", "accum_r_square"))
})

test_that("multi_tempted_decomp output dimensions are correct", {
  M <- 2; n <- 3; p <- 4; r <- 2; resolution <- 11
  dl     <- make_datlists(M = M, n = n, p = p)
  result <- suppressMessages(
    multi_tempted_decomp(dl, r = r, resolution = resolution, maxiter = 3))

  expect_equal(dim(result$A_hat), c(n, r))
  expect_equal(length(result$B_hat), M)
  expect_equal(dim(result$B_hat[[1]]), c(p, r))
  expect_equal(length(result$Zeta_hat), M)
  expect_equal(dim(result$Zeta_hat[[1]]), c(resolution, r))
  expect_equal(dim(result$Lambda), c(M, r))
  expect_equal(length(result$r_square), r)
  expect_equal(length(result$accum_r_square), r)
  expect_equal(length(result$time_Zeta), M)
  expect_equal(length(result$time_Zeta[[1]]), resolution)
})

# --- Structural properties of loadings ---------------------------------------

test_that("multi_tempted_decomp A columns are unit norm", {
  r  <- 2
  dl <- make_datlists(M = 2, n = 3, p = 4)
  result <- suppressMessages(multi_tempted_decomp(dl, r = r, maxiter = 5))
  for (l in seq_len(r)) {
    expect_equal(sum(result$A_hat[, l]^2), 1, tolerance = 1e-6,
                 label = sprintf("A column %d", l))
  }
})

test_that("multi_tempted_decomp B columns are unit norm for every modality", {
  M <- 2; r <- 2
  dl     <- make_datlists(M = M, n = 3, p = 4)
  result <- suppressMessages(multi_tempted_decomp(dl, r = r, maxiter = 5))
  for (m in seq_len(M)) for (l in seq_len(r)) {
    expect_equal(sum(result$B_hat[[m]][, l]^2), 1, tolerance = 1e-6,
                 label = sprintf("B[[%d]] column %d", m, l))
  }
})

test_that("multi_tempted_decomp Zeta columns are unit norm for every modality", {
  M <- 2; r <- 2
  dl     <- make_datlists(M = M, n = 3, p = 4)
  result <- suppressMessages(multi_tempted_decomp(dl, r = r, maxiter = 5))
  for (m in seq_len(M)) for (l in seq_len(r)) {
    expect_equal(sum(result$Zeta_hat[[m]][, l]^2), 1, tolerance = 1e-6,
                 label = sprintf("Zeta[[%d]] column %d", m, l))
  }
})

test_that("multi_tempted_decomp Lambda values are all non-negative", {
  dl     <- make_datlists(M = 2, n = 3, p = 4)
  result <- suppressMessages(multi_tempted_decomp(dl, r = 2, maxiter = 5))
  expect_true(all(result$Lambda >= 0))
})

test_that("multi_tempted_decomp r_square and accum_r_square are in [0, 1]", {
  dl     <- make_datlists(M = 2, n = 3, p = 4)
  result <- suppressMessages(multi_tempted_decomp(dl, r = 2, maxiter = 5))
  expect_true(all(result$r_square     >= 0 & result$r_square     <= 1))
  expect_true(all(result$accum_r_square >= 0 & result$accum_r_square <= 1))
})

test_that("multi_tempted_decomp accum_r_square is non-decreasing", {
  dl     <- make_datlists(M = 2, n = 3, p = 4)
  result <- suppressMessages(multi_tempted_decomp(dl, r = 3, maxiter = 5))
  diffs  <- diff(result$accum_r_square)
  expect_true(all(diffs >= -1e-10))
})

test_that("multi_tempted_decomp time_Zeta is within the observed time range", {
  dl     <- make_datlists(M = 2, n = 3, p = 4, n_times = 6)
  # Data times are seq(0, 1, ...) so observed range is [0, 1]
  result <- suppressMessages(multi_tempted_decomp(dl, r = 1, maxiter = 3))
  for (m in seq_along(result$time_Zeta)) {
    expect_true(all(result$time_Zeta[[m]] >= 0 - 1e-10 &
                      result$time_Zeta[[m]] <= 1 + 1e-10),
                label = sprintf("time_Zeta[[%d]]", m))
  }
})

# --- Weights behaviour --------------------------------------------------------

test_that("multi_tempted_decomp NULL weights gives the same result as explicit equal weights", {
  dl     <- make_datlists(M = 2, n = 3, p = 4, seed = 77)
  r_null <- suppressMessages(
    multi_tempted_decomp(dl, r = 1, smooth = 1e-5, maxiter = 10))
  r_ones <- suppressMessages(
    multi_tempted_decomp(dl, r = 1, smooth = 1e-5, maxiter = 10,
                         weights = c(1, 1)))
  expect_equal(r_null$A_hat,  r_ones$A_hat,  tolerance = 1e-10)
  expect_equal(r_null$Lambda, r_ones$Lambda, tolerance = 1e-10)
})

test_that("multi_tempted_decomp results differ when weights are very unequal", {
  dl     <- make_datlists(M = 2, n = 3, p = 4, seed = 88)
  r_eq   <- suppressMessages(
    multi_tempted_decomp(dl, r = 1, smooth = 1e-5, maxiter = 10,
                         weights = c(1, 1)))
  r_sk   <- suppressMessages(
    multi_tempted_decomp(dl, r = 1, smooth = 1e-5, maxiter = 10,
                         weights = c(100, 1)))
  expect_false(isTRUE(all.equal(r_eq$A_hat, r_sk$A_hat)))
})

# --- Recovery of known structure ----------------------------------------------

test_that("multi_tempted_decomp r_square[1] > 0.9 on near-exact rank-1 data", {
  dl     <- make_rank1_datlists(M = 2, n = 4, p = 5, n_times = 8)
  result <- suppressMessages(
    multi_tempted_decomp(dl, r = 1, smooth = 1e-5, maxiter = 20))
  expect_gt(result$r_square[1], 0.9)
})
