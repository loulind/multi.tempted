# Tests for R/multiTEMPTED_all.R
#
# multitempted_all is a wrapper around format_tempted, svd_centralize, and
# multi_tempted_decomp.  These tests focus on the wrapper's own logic:
#   - input validation
#   - subject ordering alignment across modalities
#   - the centralize toggle
#   - output structure and naming
#   - weights passthrough
#
# The correctness of the individual steps (transformation, SVD, decomposition)
# is already covered by test-decomposition.R and test-svd_centralize.R.
# Tests here use transform = "none" and tiny data (maxiter = 2) to stay fast.


# ==============================================================================
# SHARED TEST FIXTURE
# ==============================================================================

# Build a minimal, valid multitempted_all input.
# timepoints and subjectID are returned as length-M lists, matching the
# required API where each modality supplies its own sample metadata.
make_wrapper_input <- function(M = 2, n = 4, p = 5, n_times = 3, seed = 42) {
  set.seed(seed)
  n_samples <- n * n_times
  subj_id   <- rep(paste0("s", seq_len(n)), each  = n_times)
  times     <- rep(seq_len(n_times),        times = n)

  tables <- lapply(seq_len(M), function(m) {
    mat <- matrix(abs(rnorm(n_samples * p)) + 1, nrow = n_samples, ncol = p)
    colnames(mat) <- paste0("feat", seq_len(p))
    mat
  })
  names(tables) <- paste0("mod", seq_len(M))

  list(
    featuretables = tables,
    timepoints    = replicate(M, times,   simplify = FALSE),
    subjectID     = replicate(M, subj_id, simplify = FALSE)
  )
}

# Run multitempted_all with quiet messages and fast settings.
# r, resolution, and transforms are explicit parameters so tests can override
# them without triggering "matched by multiple actual arguments".
# do_ratio = FALSE keeps the output predictable for tests using random data.
run_wrapper <- function(dat, r = 1, resolution = 11, transforms = "none",
                        do_ratio = FALSE, ...) {
  suppressMessages(
    multitempted_all(dat$featuretables, dat$timepoints, dat$subjectID,
                     transforms = transforms, r = r, resolution = resolution,
                     maxiter = 2, do_ratio = do_ratio, ...))
}


# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

test_that("multitempted_all errors when timepoints list length != M", {
  dat <- make_wrapper_input(M = 2)
  # Wrap the full list in another list to make length 1 instead of 2
  expect_error(
    suppressMessages(
      multitempted_all(dat$featuretables,
                       timepoints = list(dat$timepoints[[1]]),   # length 1, M = 2
                       dat$subjectID, transforms = "none", r = 1, maxiter = 2)),
    "length M"
  )
})

test_that("multitempted_all errors when subjectID list length != M", {
  dat <- make_wrapper_input(M = 2)
  expect_error(
    suppressMessages(
      multitempted_all(dat$featuretables, dat$timepoints,
                       subjectID = c(dat$subjectID, dat$subjectID[1]),  # length 3, M = 2
                       transforms = "none", r = 1, maxiter = 2)),
    "length M"
  )
})

test_that("multitempted_all errors when transforms length != 1 and != M", {
  dat <- make_wrapper_input(M = 2)
  expect_error(
    suppressMessages(
      multitempted_all(dat$featuretables, dat$timepoints, dat$subjectID,
                       transforms = c("none", "none", "none"),   # length 3, M = 2
                       r = 1, maxiter = 2)),
    "length 1 or length M"
  )
})

test_that("multitempted_all errors when threshold length != 1 and != M", {
  dat <- make_wrapper_input(M = 2)
  expect_error(
    suppressMessages(
      multitempted_all(dat$featuretables, dat$timepoints, dat$subjectID,
                       threshold = c(0.9, 0.9, 0.9),   # length 3, M = 2
                       transforms = "none", r = 1, maxiter = 2)),
    "length 1 or length M"
  )
})

test_that("multitempted_all errors when modalities end up with different subjects", {
  dat <- make_wrapper_input(M = 2, n = 4)
  # Give modality 2 a different fourth subject
  subj_mod2 <- sub("s4", "s99", dat$subjectID[[2]])
  expect_error(
    suppressMessages(
      multitempted_all(dat$featuretables, dat$timepoints,
                       subjectID  = list(dat$subjectID[[1]], subj_mod2),
                       transforms = "none", r = 1, maxiter = 2)),
    "different subjects"
  )
})


# ==============================================================================
# NAMING
# ==============================================================================

test_that("unnamed featuretables get default names 'modality1', 'modality2', ...", {
  dat <- make_wrapper_input(M = 2)
  names(dat$featuretables) <- NULL
  result <- run_wrapper(dat)
  expect_equal(names(result$datlists), c("modality1", "modality2"))
})

test_that("featuretable names propagate to datlists, B_hat, Zeta_hat, time_Zeta", {
  dat <- make_wrapper_input(M = 2)
  names(dat$featuretables) <- c("stool", "sputum")
  result <- run_wrapper(dat)
  expect_equal(names(result$datlists),  c("stool", "sputum"))
  expect_equal(names(result$B_hat),     c("stool", "sputum"))
  expect_equal(names(result$Zeta_hat),  c("stool", "sputum"))
  expect_equal(names(result$time_Zeta), c("stool", "sputum"))
})

test_that("modalities with different numbers of time points per subject are accepted", {
  # The per-modality list API is specifically designed to support this case:
  # different subjects × time structure in each modality.
  set.seed(99)
  n <- 3; p <- 4

  make_mod <- function(n_tp) {
    n_samp <- n * n_tp
    list(
      table    = {m <- matrix(abs(rnorm(n_samp * p)) + 1, n_samp, p);
      colnames(m) <- paste0("f", seq_len(p)); m},
      times    = rep(seq_len(n_tp), times = n),
      subjects = rep(paste0("s", seq_len(n)), each = n_tp)
    )
  }
  m1 <- make_mod(3)
  m2 <- make_mod(4)   # different number of time points

  result <- suppressMessages(
    multitempted_all(
      featuretables = list(mod1 = m1$table, mod2 = m2$table),
      timepoints    = list(m1$times,    m2$times),
      subjectID     = list(m1$subjects, m2$subjects),
      transforms = "none", r = 1, resolution = 11, maxiter = 2))

  expect_equal(dim(result$A_hat), c(n, 1))
})


# ==============================================================================
# SUBJECT ORDERING
# ==============================================================================

test_that("subjects in different row orders across modalities are aligned consistently", {
  dat <- make_wrapper_input(M = 2, n = 4, n_times = 3)

  # Shuffle all rows of modality 2 (keep subjectID and timepoints matched)
  set.seed(7)
  shuf <- sample(nrow(dat$featuretables[[2]]))
  dat$featuretables[[2]] <- dat$featuretables[[2]][shuf, ]
  subj_mod2 <- dat$subjectID[[2]][shuf]
  time_mod2 <- dat$timepoints[[2]][shuf]

  result <- suppressMessages(
    multitempted_all(dat$featuretables,
                     timepoints = list(dat$timepoints[[1]], time_mod2),
                     subjectID  = list(dat$subjectID[[1]],  subj_mod2),
                     transforms = "none", r = 1, resolution = 11, maxiter = 2))

  # After alignment, every modality's datlists must have subjects in the same order
  for (m in seq_along(result$datlists)) {
    expect_equal(names(result$datlists[[m]]), names(result$datlists[[1]]))
  }
})


# ==============================================================================
# OUTPUT STRUCTURE
# ==============================================================================

test_that("multitempted_all without ratio returns the core decomposition elements", {
  dat    <- make_wrapper_input()
  result <- run_wrapper(dat, do_ratio = FALSE)
  core <- c("datlists", "mean_svd", "A_hat", "B_hat", "Zeta_hat",
            "time_Zeta", "Lambda", "r_square", "accum_r_square",
            "metafeature_aggregate", "toppct_aggregate", "contrast")
  expect_true(all(core %in% names(result)))
  expect_false("metafeature_ratio" %in% names(result))
})

test_that("multitempted_all with do_ratio = TRUE adds ratio output elements", {
  dat    <- make_wrapper_input()
  result <- run_wrapper(dat, do_ratio = TRUE)
  expect_true("metafeature_ratio"  %in% names(result))
  expect_true("toppct_ratio"       %in% names(result))
  expect_true("bottompct_ratio"    %in% names(result))
})

test_that("A_hat is an n x r matrix", {
  M <- 2; n <- 4; r <- 1
  dat    <- make_wrapper_input(M = M, n = n)
  result <- run_wrapper(dat, r = r)
  expect_equal(dim(result$A_hat), c(n, r))
})

test_that("B_hat is a length-M list of p_m x r matrices", {
  M <- 2; p <- 5; r <- 1
  dat    <- make_wrapper_input(M = M, p = p)
  result <- run_wrapper(dat, r = r)
  expect_equal(length(result$B_hat), M)
  for (m in seq_len(M)) {
    expect_equal(dim(result$B_hat[[m]]), c(p, r),
                 label = sprintf("B_hat[[%d]] dimensions", m))
  }
})

test_that("Zeta_hat is a length-M list of resolution x r matrices", {
  M <- 2; resolution <- 11; r <- 1
  dat    <- make_wrapper_input(M = M)
  result <- run_wrapper(dat, r = r, resolution = resolution)
  expect_equal(length(result$Zeta_hat), M)
  for (m in seq_len(M)) {
    expect_equal(dim(result$Zeta_hat[[m]]), c(resolution, r),
                 label = sprintf("Zeta_hat[[%d]] dimensions", m))
  }
})

test_that("Lambda is an M x r matrix", {
  M <- 2; r <- 1
  dat    <- make_wrapper_input(M = M)
  result <- run_wrapper(dat, r = r)
  expect_equal(dim(result$Lambda), c(M, r))
})

test_that("r_square and accum_r_square are M x r matrices", {
  M <- 2; r <- 1
  dat    <- make_wrapper_input(M = M)
  result <- run_wrapper(dat, r = r)
  expect_equal(dim(result$r_square),       c(M, r))
  expect_equal(dim(result$accum_r_square), c(M, r))
})

test_that("time_Zeta is a length-M list of length-resolution vectors", {
  M <- 2; resolution <- 11
  dat    <- make_wrapper_input(M = M)
  result <- run_wrapper(dat, resolution = resolution)
  expect_equal(length(result$time_Zeta), M)
  for (m in seq_len(M)) {
    expect_equal(length(result$time_Zeta[[m]]), resolution)
  }
})


# ==============================================================================
# CENTRALIZE TOGGLE
# ==============================================================================

test_that("centralize = TRUE produces a non-NULL mean_svd with the expected structure", {
  dat    <- make_wrapper_input()
  result <- run_wrapper(dat, centralize = TRUE, r_svd = 1)
  expect_false(is.null(result$mean_svd))
  expect_named(result$mean_svd, c("datlists", "A_tilde", "B_tilde", "lambda_tilde"))
})

test_that("centralize = FALSE leaves mean_svd as NULL", {
  dat    <- make_wrapper_input()
  result <- run_wrapper(dat, centralize = FALSE)
  expect_null(result$mean_svd)
})

test_that("datlists in output is the pre-centralisation data, not the centralised data", {
  # result$datlists should match format_tempted output directly.
  # result$mean_svd$datlists (centralised) should differ in feature rows.
  dat    <- make_wrapper_input(M = 1)
  result <- run_wrapper(dat, centralize = TRUE, r_svd = 1)

  # Direct format_tempted call for reference — use [[1]] since dat$timepoints
  # and dat$subjectID are now length-M lists.
  expected_datlist <- format_tempted(featuretable = dat$featuretables[[1]],
                                     timepoint    = dat$timepoints[[1]],
                                     subjectID    = dat$subjectID[[1]],
                                     transform    = "none")

  for (i in seq_along(expected_datlist)) {
    expect_equal(result$datlists[[1]][[i]][-1, ],
                 expected_datlist[[i]][-1, ],
                 label = sprintf("datlists subject %d features", i))
  }

  # Centralised feature rows should differ from the raw formatted values
  expect_false(isTRUE(all.equal(result$mean_svd$datlists[[1]][[1]][-1, ],
                                result$datlists[[1]][[1]][-1, ])))
})


# ==============================================================================
# WEIGHTS PASSTHROUGH
# ==============================================================================

test_that("NULL weights gives the same result as explicit equal weights", {
  dat     <- make_wrapper_input(seed = 55)
  r_null  <- run_wrapper(dat, weights = NULL)
  r_equal <- run_wrapper(dat, weights = c(1, 1))
  expect_equal(r_null$A_hat,  r_equal$A_hat,  tolerance = 1e-10)
  expect_equal(r_null$Lambda, r_equal$Lambda, tolerance = 1e-10)
})

test_that("very unequal weights give a different result from equal weights", {
  dat     <- make_wrapper_input(seed = 66)
  r_equal <- run_wrapper(dat, weights = c(1, 1))
  r_skew  <- run_wrapper(dat, weights = c(100, 1))
  expect_false(isTRUE(all.equal(r_equal$A_hat, r_skew$A_hat)))
})


# ==============================================================================
# PER-MODALITY TRANSFORMS
# ==============================================================================

test_that("per-modality transforms (vector) give same result as scalar when equal", {
  dat      <- make_wrapper_input()
  r_scalar <- run_wrapper(dat, transforms = "none")
  r_vector <- suppressMessages(
    multitempted_all(dat$featuretables, dat$timepoints, dat$subjectID,
                     transforms = c("none", "none"),
                     r = 1, resolution = 11, maxiter = 2))
  expect_equal(r_scalar$A_hat,    r_vector$A_hat)
  expect_equal(r_scalar$datlists, r_vector$datlists)
})
