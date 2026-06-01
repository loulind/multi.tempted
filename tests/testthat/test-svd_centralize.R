# Tests for R/svd_centralize.R


# ==============================================================================
# SHARED TEST FIXTURES
# ==============================================================================

# Minimal valid datlists: M modalities, n subjects, p features, n_times time points.
make_datlists_svd <- function(M = 2, n = 4, p = 5, n_times = 6, seed = 42) {
  set.seed(seed)
  times <- seq(0, 1, length.out = n_times)
  mods <- lapply(seq_len(M), function(m) {
    subjs <- lapply(seq_len(n), function(i)
      rbind(times, matrix(rnorm(p * n_times), nrow = p)))
    names(subjs) <- paste0("s", seq_len(n))
    subjs
  })
  names(mods) <- paste0("mod", seq_len(M))
  mods
}

# Convenience: compute the time-averaged feature matrix for modality m.
# Returns an n x p_m matrix where entry [i, j] is the mean of feature j for subject i.
time_mean_matrix <- function(datlists, m) {
  n   <- length(datlists[[m]])
  p_m <- nrow(datlists[[m]][[1]]) - 1
  mat <- matrix(0, n, p_m)
  for (i in seq_len(n)) {
    mat[i, ] <- rowMeans(datlists[[m]][[i]][-1, , drop = FALSE])
  }
  mat
}


# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

test_that("svd_centralize errors when modalities have different subject counts", {
  dl      <- make_datlists_svd(M = 2, n = 4)
  dl[[2]] <- dl[[2]][1:3]   # drop one subject from modality 2
  expect_error(svd_centralize(dl), "same number of subjects")
})

test_that("svd_centralize errors when feature counts differ across subjects within a modality", {
  dl <- make_datlists_svd(M = 1, n = 3, p = 4)
  # Give subject 2 an extra feature row
  dl[[1]][[2]] <- rbind(dl[[1]][[2]], rnorm(6))
  expect_error(svd_centralize(dl), "inconsistent feature counts")
})


# ==============================================================================
# OUTPUT STRUCTURE
# ==============================================================================

test_that("svd_centralize returns a list with the four expected names", {
  dl     <- make_datlists_svd()
  result <- svd_centralize(dl, r = 1)
  expect_named(result, c("datlists", "A_tilde", "B_tilde", "lambda_tilde"))
})

test_that("svd_centralize output datlists preserves the input structure", {
  M <- 2; n <- 4; p <- 5; n_times <- 6
  dl     <- make_datlists_svd(M = M, n = n, p = p, n_times = n_times)
  result <- svd_centralize(dl, r = 1)
  expect_equal(length(result$datlists), M)
  for (m in seq_len(M)) {
    expect_equal(length(result$datlists[[m]]), n)
    for (i in seq_len(n)) {
      expect_equal(dim(result$datlists[[m]][[i]]), c(p + 1, n_times))
    }
  }
})

test_that("svd_centralize A_tilde dimensions are n x r for every modality", {
  M <- 2; n <- 4; r <- 2
  dl     <- make_datlists_svd(M = M, n = n)
  result <- svd_centralize(dl, r = r)
  expect_equal(length(result$A_tilde), M)
  for (m in seq_len(M)) {
    expect_equal(dim(result$A_tilde[[m]]), c(n, r))
  }
})

test_that("svd_centralize B_tilde dimensions are p_m x r for every modality", {
  M <- 2; p <- 5; r <- 2
  dl     <- make_datlists_svd(M = M, p = p)
  result <- svd_centralize(dl, r = r)
  expect_equal(length(result$B_tilde), M)
  for (m in seq_len(M)) {
    expect_equal(dim(result$B_tilde[[m]]), c(p, r))
  }
})

test_that("svd_centralize lambda_tilde is a length-r vector for every modality", {
  M <- 2; r <- 2
  dl     <- make_datlists_svd(M = M)
  result <- svd_centralize(dl, r = r)
  expect_equal(length(result$lambda_tilde), M)
  for (m in seq_len(M)) {
    expect_equal(length(result$lambda_tilde[[m]]), r)
  }
})

test_that("svd_centralize propagates modality names to all output lists", {
  dl     <- make_datlists_svd(M = 2)
  result <- svd_centralize(dl, r = 1)
  expect_equal(names(result$datlists),     names(dl))
  expect_equal(names(result$A_tilde),      names(dl))
  expect_equal(names(result$B_tilde),      names(dl))
  expect_equal(names(result$lambda_tilde), names(dl))
})


# ==============================================================================
# TIME ROW IS UNCHANGED
# ==============================================================================

test_that("svd_centralize does not alter the time row (row 1) of any subject matrix", {
  dl     <- make_datlists_svd(M = 2, n = 3, p = 4)
  result <- svd_centralize(dl, r = 1)
  for (m in seq_len(2)) {
    for (i in seq_len(3)) {
      expect_equal(result$datlists[[m]][[i]][1, ],
                   dl[[m]][[i]][1, ],
                   label = sprintf("time row: modality %d, subject %d", m, i))
    }
  }
})


# ==============================================================================
# MATHEMATICAL CORRECTNESS
# ==============================================================================

test_that("centralised data + rank-r mean reconstruction equals original feature rows", {
  # The rank-r mean removed from subject i is mean_rank_r[i, ], broadcast over
  # time. Adding it back to the centralised features must recover the original.
  M <- 2; n <- 3; p <- 4; r <- 1
  dl     <- make_datlists_svd(M = M, n = n, p = p)
  result <- svd_centralize(dl, r = r)

  for (m in seq_len(M)) {
    # Reconstruct the rank-r mean matrix: U %*% diag(d) %*% t(V)
    mean_rank_r <- result$A_tilde[[m]] %*%
      t(result$B_tilde[[m]] * result$lambda_tilde[[m]])  # n x p_m

    for (i in seq_len(n)) {
      recovered <- result$datlists[[m]][[i]][-1, ] + mean_rank_r[i, ]
      expect_equal(recovered, dl[[m]][[i]][-1, ], tolerance = 1e-12,
                   label = sprintf("reconstruction: modality %d, subject %d", m, i))
    }
  }
})

test_that("when r >= min(n, p_m) the time-mean of centralised features is zero", {
  # A full-rank SVD removes the entire mean matrix, so every subject's
  # time-averaged feature vector should be zero after centralisation.
  n <- 3; p <- 5
  r_full <- n   # r = n removes all variance in the n x p mean matrix
  dl     <- make_datlists_svd(M = 1, n = n, p = p)
  result <- svd_centralize(dl, r = r_full)

  residual_mean <- time_mean_matrix(result$datlists, m = 1)
  expect_equal(residual_mean, matrix(0, n, p), tolerance = 1e-10)
})

test_that("with r = 1 the residual mean matrix is orthogonal to the retained singular vectors", {
  # The residual mean (after removing rank-1 approximation) must be orthogonal
  # to A_tilde[[m]] and B_tilde[[m]], since those span the removed subspace.
  M <- 2; n <- 4; p <- 5; r <- 1
  dl     <- make_datlists_svd(M = M, n = n, p = p)
  result <- svd_centralize(dl, r = r)

  for (m in seq_len(M)) {
    residual_mean <- time_mean_matrix(result$datlists, m)   # n x p_m
    # Projection of residual onto left singular vectors should be ~0
    proj_left  <- t(result$A_tilde[[m]]) %*% residual_mean        # r x p_m
    expect_equal(proj_left, matrix(0, r, p), tolerance = 1e-10,
                 label = sprintf("left projection: modality %d", m))
    # Projection of residual onto right singular vectors should be ~0
    proj_right <- residual_mean %*% result$B_tilde[[m]]            # n x r
    expect_equal(proj_right, matrix(0, n, r), tolerance = 1e-10,
                 label = sprintf("right projection: modality %d", m))
  }
})

test_that("removing more ranks reduces the norm of the residual mean", {
  # Higher r removes more of the mean structure, so the Frobenius norm of the
  # residual mean matrix should be non-increasing as r grows.
  n <- 5; p <- 6
  dl  <- make_datlists_svd(M = 1, n = n, p = p, seed = 7)

  norms <- sapply(seq_len(min(n, p)), function(r) {
    res <- svd_centralize(dl, r = r)
    norm(time_mean_matrix(res$datlists, m = 1), type = "F")
  })
  expect_true(all(diff(norms) <= 1e-10))
})


# ==============================================================================
# SVD COMPONENT PROPERTIES
# ==============================================================================

test_that("A_tilde columns are orthonormal for every modality", {
  M <- 2; r <- 2
  dl     <- make_datlists_svd(M = M, n = 5, p = 6)
  result <- svd_centralize(dl, r = r)
  for (m in seq_len(M)) {
    gram <- t(result$A_tilde[[m]]) %*% result$A_tilde[[m]]
    expect_equal(gram, diag(r), tolerance = 1e-10,
                 label = sprintf("A_tilde orthonormality: modality %d", m))
  }
})

test_that("B_tilde columns are orthonormal for every modality", {
  M <- 2; r <- 2
  dl     <- make_datlists_svd(M = M, n = 5, p = 6)
  result <- svd_centralize(dl, r = r)
  for (m in seq_len(M)) {
    gram <- t(result$B_tilde[[m]]) %*% result$B_tilde[[m]]
    expect_equal(gram, diag(r), tolerance = 1e-10,
                 label = sprintf("B_tilde orthonormality: modality %d", m))
  }
})

test_that("lambda_tilde values are positive and non-increasing for every modality", {
  M <- 2; r <- 3
  dl     <- make_datlists_svd(M = M, n = 5, p = 6)
  result <- svd_centralize(dl, r = r)
  for (m in seq_len(M)) {
    lam <- result$lambda_tilde[[m]]
    expect_true(all(lam > 0),
                label = sprintf("lambda positive: modality %d", m))
    expect_true(all(diff(lam) <= 0),
                label = sprintf("lambda non-increasing: modality %d", m))
  }
})


# ==============================================================================
# BEHAVIOUR WITH r = 1 (DEFAULT)
# ==============================================================================

test_that("svd_centralize with default r = 1 produces the same result as explicit r = 1", {
  dl       <- make_datlists_svd(M = 2, n = 4, p = 5, seed = 99)
  r_default  <- svd_centralize(dl)
  r_explicit <- svd_centralize(dl, r = 1)
  expect_equal(r_default$datlists,     r_explicit$datlists)
  expect_equal(r_default$lambda_tilde, r_explicit$lambda_tilde)
})
