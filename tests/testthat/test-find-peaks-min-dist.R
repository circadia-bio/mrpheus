# tests/testthat/test-find-peaks-min-dist.R
#
# Regression tests for `.find_peaks_min_dist()`, the internal helper behind
# detect_qrs() Stage 5 peak selection. Fixed in mrpheus 0.1.4.9000 to replace
# an O(n*k) candidate-vs-candidate scan with a binary-search insertion into a
# pre-allocated sorted buffer (see NEWS.md).
#
# All fixtures here are fully synthetic (no participant-derived data).

# ── Brute-force reference (pre-fix semantics; used only to check
#    correctness, deliberately a different code path) ───────────────────────

.find_peaks_min_dist_ref <- function(x, min_dist = 1L) {
  n <- length(x)
  if (n < 3L) return(list(pks = double(0), locs = integer(0)))

  is_peak  <- x[2:(n - 1L)] >= x[1:(n - 2L)] & x[2:(n - 1L)] >= x[3:n]
  locs_all <- which(is_peak) + 1L
  if (length(locs_all) == 0L) return(list(pks = double(0), locs = integer(0)))

  pks_all <- x[locs_all]
  ord     <- order(pks_all, decreasing = TRUE)

  accepted <- integer(0)
  for (k in ord) {
    loc <- locs_all[k]
    if (length(accepted) == 0L || all(abs(accepted - loc) >= min_dist)) {
      accepted <- c(accepted, loc)
    }
  }
  final <- sort(accepted)
  list(pks = x[final], locs = final)
}

# ── Correctness vs brute-force reference ─────────────────────────────────────

test_that(".find_peaks_min_dist matches brute-force reference on random signal", {
  set.seed(101)
  x <- cumsum(rnorm(2000))  # smooth-ish random walk, plenty of local maxima
  fast <- .find_peaks_min_dist(x, min_dist = 10L)
  ref  <- .find_peaks_min_dist_ref(x, min_dist = 10L)

  expect_equal(fast$locs, ref$locs)
  expect_equal(fast$pks, ref$pks)
})

test_that(".find_peaks_min_dist matches brute-force reference with small min_dist", {
  set.seed(102)
  x <- sin(seq(0, 40 * pi, length.out = 3000)) + rnorm(3000, sd = 0.05)
  fast <- .find_peaks_min_dist(x, min_dist = 3L)
  ref  <- .find_peaks_min_dist_ref(x, min_dist = 3L)

  expect_equal(fast$locs, ref$locs)
  expect_equal(fast$pks, ref$pks)
})

test_that(".find_peaks_min_dist enforces the minimum distance constraint", {
  set.seed(103)
  x <- cumsum(rnorm(5000))
  out <- .find_peaks_min_dist(x, min_dist = 20L)
  expect_true(all(diff(out$locs) >= 20L))
})

test_that(".find_peaks_min_dist keeps the tallest peak in a close cluster", {
  x <- c(0, 1, 0, 5, 0, 2, 0, 1, 0)  # local maxima at 2,4,6,8 (values 1,5,2,1)
  out <- .find_peaks_min_dist(x, min_dist = 5L)
  expect_true(4L %in% out$locs)
  expect_false(2L %in% out$locs)
  expect_false(6L %in% out$locs)
})

test_that(".find_peaks_min_dist returns empty result for short signals", {
  expect_equal(.find_peaks_min_dist(c(1, 2), min_dist = 1L)$locs, integer(0))
})

test_that(".find_peaks_min_dist handles a flat (all-tied) signal via stable tie-break", {
  # Candidate detection uses `>=` on both neighbours, so on a perfectly flat
  # signal every interior point ties as a "peak" (pre-existing plateau
  # handling, not touched by this fix). Ties are broken by original index
  # order, so the greedy min-distance selection should pick evenly spaced
  # points starting at the first candidate (index 2).
  out <- .find_peaks_min_dist(rep(1, 100), min_dist = 5L)
  expect_true(length(out$locs) > 0L)
  expect_equal(out$locs[1], 2L)
  expect_true(all(diff(out$locs) == 5L))
  expect_true(all(out$pks == 1))
})

test_that(".find_peaks_min_dist output is ordered by location", {
  set.seed(104)
  x <- cumsum(rnorm(1000))
  out <- .find_peaks_min_dist(x, min_dist = 5L)
  expect_true(all(diff(out$locs) > 0))
})

# ── Scale check: regression test for the original bug. A pure O(n^2) or
#    O(n*k) implementation would be unusably slow or hang here; the fixed
#    version should return well within a few seconds. ───────────────────────

test_that(".find_peaks_min_dist completes quickly on a large, realistically-smoothed candidate set", {
  set.seed(105)
  # A raw random walk (cumsum(rnorm(n))) has far more local extrema (~1/3 of
  # points) than the signal this helper actually sees in production: ecg_m,
  # which has already been through 150ms boxcar smoothing in detect_qrs()
  # Stage 4. Smoothing here first keeps the candidate density realistic;
  # otherwise this becomes an adversarial worst case dominated by pure R
  # for-loop/findInterval() overhead rather than a meaningful regression
  # check against the original O(n*k) blow-up.
  n   <- 470000  # comparable to the 15-min/500Hz recording that originally hung
  win <- 75L
  raw <- cumsum(rnorm(n))
  x   <- as.numeric(stats::filter(raw, rep(1 / win, win), sides = 2))
  x   <- x[!is.na(x)]

  t0  <- Sys.time()
  out <- .find_peaks_min_dist(x, min_dist = 100L)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  expect_true(all(diff(out$locs) >= 100L))
  expect_lt(elapsed, 20)  # generous ceiling; pre-fix this would not return at all
})
