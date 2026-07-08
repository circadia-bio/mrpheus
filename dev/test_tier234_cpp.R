# dev/test_tier234_cpp.R
#
# Parity validation for Tier 2-4 Rcpp functions:
#   rowmedian_cpp, roll_right_mean_cpp, robust_scale_cpp  (staging pipeline)
#   roll_rms_cpp, detect_so_candidates_cpp                (event detection)
#
# Run after:
#   Rcpp::compileAttributes()
#   devtools::load_all()

library(mrpheus)

cat("Tier 2-4 Rcpp parity validation\n")
cat("================================\n\n")

check <- function(label, cpp_val, r_val, tol = 1e-10) {
  ok <- isTRUE(all.equal(cpp_val, r_val, tolerance = tol, check.names = FALSE))
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", label))
  if (!ok) { cat("  C++:", head(cpp_val, 5), "\n"); cat("  R:  ", head(r_val,  5), "\n") }
  invisible(ok)
}

# ── rowmedian_cpp ─────────────────────────────────────────────────────────────
cat("--- rowmedian_cpp ---\n")
set.seed(1)
m1 <- matrix(rnorm(251 * 9), nrow = 251)  # typical pgrams size
check("251x9 matrix", rowmedian_cpp(m1), apply(m1, 1L, median))

m2 <- matrix(rnorm(251 * 1), nrow = 251)  # single segment
check("251x1 matrix (n_segs=1)", rowmedian_cpp(m2), apply(m2, 1L, median))

m3 <- matrix(c(1,3,5,2,4,6), nrow = 2)
check("2x3 matrix: row medians c(4, 3)", rowmedian_cpp(m3), c(4, 3))

# ── roll_right_mean_cpp ───────────────────────────────────────────────────────
cat("\n--- roll_right_mean_cpp ---\n")

.roll_right_r <- function(x, k) {
  as.vector(zoo::rollapply(x, k, mean, fill = NA, partial = TRUE, align = "right"))
}

set.seed(2)
x <- rnorm(200)
check("k=4 random (n=200)", roll_right_mean_cpp(x, 4L), .roll_right_r(x, 4L))
check("k=1 is identity",    roll_right_mean_cpp(x, 1L), x)
check("k=4 partial window at i=1: mean(x[1])",
      roll_right_mean_cpp(x, 4L)[1],  x[1])
check("k=4 partial window at i=2: mean(x[1:2])",
      roll_right_mean_cpp(x, 4L)[2],  mean(x[1:2]))
check("k=4 full window at i=4: mean(x[1:4])",
      roll_right_mean_cpp(x, 4L)[4],  mean(x[1:4]))
check("constant input",
      roll_right_mean_cpp(rep(7, 50), 4L), rep(7, 50))

# ── robust_scale_cpp ──────────────────────────────────────────────────────────
cat("\n--- robust_scale_cpp ---\n")

.robust_r <- function(x, q_low = 0.05, q_high = 0.95) {
  eps <- 1e-10
  med <- median(x, na.rm = TRUE)
  q   <- quantile(x, c(q_low, q_high), na.rm = TRUE)
  (x - med) / (q[2L] - q[1L] + eps)
}

set.seed(3)
x <- rnorm(1000)
check("random n=1000 defaults",
      robust_scale_cpp(x),        .robust_r(x))
check("custom q_low=0.1, q_high=0.9",
      robust_scale_cpp(x, 0.1, 0.9), .robust_r(x, 0.1, 0.9))
check("median of result is 0",
      median(robust_scale_cpp(x)), 0, tol = 1e-6)

x_na <- c(1, 2, NA, 4, 5)
r_na <- robust_scale_cpp(x_na)
check("NA preserved at position 3", is.na(r_na[3]), TRUE)
check("non-NA values scaled correctly (NA excluded from stats)",
      r_na[!is.na(r_na)], .robust_r(x_na[!is.na(x_na)]), tol = 1e-9)

check("constant input -> all zeros",
      unique(robust_scale_cpp(rep(5, 50))), 0, tol = 1e-5)

# ── roll_rms_cpp ──────────────────────────────────────────────────────────────
cat("\n--- roll_rms_cpp ---\n")

.roll_rms_r <- function(x, k) {
  as.vector(zoo::rollapply(x^2, k, function(w) sqrt(mean(w)), fill = NA, align = "center"))
}

set.seed(4)
x <- rnorm(300)
for (k in c(3L, 5L, 15L, 30L)) {
  check(sprintf("k=%d random (n=300)", k),
        roll_rms_cpp(x, k), .roll_rms_r(x, k), tol = 1e-10)
}
check("k=5 edge positions are NA",
      is.na(roll_rms_cpp(x, 5L)[1:2]), c(TRUE, TRUE))
check("constant input: rms = |constant|",
      roll_rms_cpp(rep(3, 100), 5L)[3:98],
      rep(3, 96), tol = 1e-10)

# ── detect_so_candidates_cpp ──────────────────────────────────────────────────
cat("\n--- detect_so_candidates_cpp ---\n")

.so_r <- function(sig, sr, dur_neg, dur_pos, amp) {
  zc <- which(diff(sign(sig)) != 0)
  if (length(zc) < 4) return(matrix(numeric(0), ncol = 6))
  rows <- lapply(seq(1, length(zc) - 3, by = 2), function(k) {
    ns <- zc[k]; ne <- zc[k+1]; pe <- zc[k+2]
    dn <- (ne - ns) / sr; dp <- (pe - ne) / sr
    if (dn < dur_neg[1] || dn > dur_neg[2]) return(NULL)
    if (dp < dur_pos[1] || dp > dur_pos[2]) return(NULL)
    np <- min(sig[ns:ne]); pp <- max(sig[ne:pe]); ptp <- pp - np
    if (ptp < amp[1] || ptp > amp[2]) return(NULL)
    c(ns, ne, pe, np, pp, ptp)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(matrix(numeric(0), ncol = 6))
  do.call(rbind, rows)
}

# Construct a synthetic signal with one clean SO
sr   <- 256
t    <- seq(0, 4, by = 1/sr)
sig  <- -sin(2 * pi * 1 * t) * 200   # 1 Hz, ±200 uV

cpp_cands <- detect_so_candidates_cpp(
  sig, sr, 0.1, 1.5, 0.1, 1.0, 50, 500
)
r_cands <- .so_r(sig, sr, c(0.1, 1.5), c(0.1, 1.0), c(50, 500))

check("same number of candidates",    nrow(cpp_cands), nrow(r_cands))
if (nrow(cpp_cands) > 0 && nrow(r_cands) > 0) {
  check("neg_start matches", cpp_cands[,1], r_cands[,1])
  check("pos_end   matches", cpp_cands[,3], r_cands[,3])
  check("ptp       matches", cpp_cands[,6], r_cands[,6], tol = 1e-8)
}

# Flat signal: no candidates
flat_cands <- detect_so_candidates_cpp(rep(0, 1000), sr, 0.1, 1.5, 0.1, 1.0, 50, 500)
check("flat signal -> 0 candidates", nrow(flat_cands), 0L)

# Very short signal: no candidates
short_cands <- detect_so_candidates_cpp(c(1, -1, 1), sr, 0.1, 1.5, 0.1, 1.0, 0, 1e6)
check("< 4 zero crossings -> 0 candidates", nrow(short_cands), 0L)

cat("\nDone.\n")
