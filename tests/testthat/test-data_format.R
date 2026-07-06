# Tests for R/data_format.R (format_tempted)
#
# format_tempted transforms a sample-by-feature table plus sample metadata into
# the per-subject list-of-matrices format consumed by the decomposition. These
# tests focus on:
#   - input validation
#   - each transform producing finite output with the expected structure
#   - the "lfb" (log-2 fold change over baseline) path specifically, which used
#     to be shadowed by the "logcomp" branch


# ==============================================================================
# SHARED FIXTURE
# ==============================================================================

# A small positive count table with 2 subjects x 3 time points and p features.
make_counts <- function(p = 4, seed = 1) {
  set.seed(seed)
  n_samples <- 6
  counts <- matrix(rpois(n_samples * p, lambda = 20) + 1, nrow = n_samples)
  colnames(counts) <- paste0("f", seq_len(p))
  list(
    counts    = counts,
    timepoint = c(0, 1, 2, 0, 1, 2),
    subjectID = c("a", "a", "a", "b", "b", "b")
  )
}


# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

test_that("format_tempted errors on subjects with a single time point", {
  d <- make_counts()
  sid <- d$subjectID
  sid[4:6] <- "c"; sid[4] <- "d"   # subject 'd' has one sample
  expect_error(format_tempted(d$counts, d$timepoint, sid, transform = "none"),
               "only one time point")
})

test_that("format_tempted errors on subjectID / timepoint length mismatch", {
  d <- make_counts()
  expect_error(format_tempted(d$counts, d$timepoint, d$subjectID[-1]),
               "subjectID does not match")
  expect_error(format_tempted(d$counts, d$timepoint[-1], d$subjectID),
               "timepoint does not match")
})

test_that("format_tempted errors on an unknown transform", {
  d <- make_counts()
  expect_error(format_tempted(d$counts, d$timepoint, d$subjectID,
                              transform = "not_a_transform"),
               "Unknown transform")
})

test_that("format_tempted errors when a log transform is given negative data", {
  d <- make_counts()
  neg <- d$counts; neg[1, 1] <- -5
  expect_error(format_tempted(neg, d$timepoint, d$subjectID, transform = "clr"),
               "non-negative")
})


# ==============================================================================
# TRANSFORMS PRODUCE FINITE, WELL-STRUCTURED OUTPUT
# ==============================================================================

test_that("every non-lfb transform returns finite per-subject matrices", {
  d <- make_counts(p = 4)
  for (tf in c("logcomp", "comp", "ast", "clr", "logit", "none")) {
    dl <- format_tempted(d$counts, d$timepoint, d$subjectID, transform = tf)
    expect_equal(length(dl), 2, label = sprintf("n subjects (%s)", tf))
    for (i in seq_along(dl)) {
      expect_true(all(is.finite(dl[[i]])),
                  label = sprintf("finite values (%s, subject %d)", tf, i))
      # (p features + 1 time row) x 3 time points
      expect_equal(dim(dl[[i]]), c(5, 3),
                   label = sprintf("dims (%s, subject %d)", tf, i))
    }
  }
})


# ==============================================================================
# lfb: log-2 fold change over baseline
# ==============================================================================

test_that("lfb runs, drops the baseline column, and returns finite values", {
  d  <- make_counts(p = 4)
  dl <- format_tempted(d$counts, d$timepoint, d$subjectID, transform = "lfb")
  for (i in seq_along(dl)) {
    # baseline (first) time point is subtracted out and dropped: 3 -> 2 columns
    expect_equal(ncol(dl[[i]]), 2, label = sprintf("lfb n columns (subject %d)", i))
    expect_true(all(is.finite(dl[[i]])),
                label = sprintf("lfb finite (subject %d)", i))
  }
})

test_that("lfb feature values equal log2 composition differenced from baseline", {
  d <- make_counts(p = 3)
  # Reconstruct the expected lfb values for subject 'a' by hand.
  pseudo <- apply(d$counts, 1, function(x) min(x[x != 0]) / 2)
  comp_log2 <- log2((d$counts + pseudo) / rowSums(d$counts + pseudo))
  a_rows <- which(d$subjectID == "a")            # already time-ordered 0,1,2
  expected <- t(comp_log2[a_rows, ])             # features x time
  expected <- expected[, -1, drop = FALSE] - expected[, 1]

  dl <- format_tempted(d$counts, d$timepoint, d$subjectID, transform = "lfb")
  expect_equal(unname(dl[["a"]][-1, ]), unname(expected), tolerance = 1e-10)
})
