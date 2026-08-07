# detect_qrs.R
# R port of pan_tompkin.m (Hooman Sedghamiz, Linköping University, 2018)
# See inst/licenses/pan_tompkin-LICENSE.txt

#' Detect QRS Complexes in an ECG Signal (Pan-Tompkins Algorithm)
#'
#' Detects R-peaks in a raw ECG signal using the Pan-Tompkins algorithm.
#' Implements bandpass filtering, derivative, squaring, moving-average
#' integration, and adaptive dual-threshold peak detection with T-wave
#' rejection and search-back for missed beats.
#'
#' @param ecg  Numeric vector. Raw ECG signal.
#' @param fs   Numeric. Sampling frequency in Hz. If `ecg` comes from
#'   [mrpheus::read_philips_physlog()], use `rec$HDR$sfreq`.
#'
#' @return A list of class `mrpheus_qrs` with components:
#' \describe{
#'   \item{`qrs_i`}{Integer vector of R-peak sample indices (1-based),
#'     referenced to the bandpass-filtered signal.}
#'   \item{`qrs_amp`}{Numeric vector of R-wave amplitudes at each detected
#'     peak (normalised bandpass signal units).}
#'   \item{`delay`}{Numeric. Sample delay introduced by the moving-average
#'     stage (half the 150 ms window length).}
#'   \item{`fs`}{Numeric. Sampling frequency passed to the function.}
#' }
#'
#' @details
#' Processing pipeline (Pan & Tompkins, 1985):
#' 1. **Bandpass filter** 5–15 Hz (Butterworth N=3). At fs = 200 Hz,
#'    implemented as separate 12 Hz low-pass then 5 Hz high-pass stages.
#' 2. **Derivative filter** — 5-point kernel interpolated to match fs.
#' 3. **Squaring** — non-linear emphasis of dominant peaks.
#' 4. **150 ms moving-average integration.**
#' 5. **Adaptive thresholding** — dual-threshold scheme with T-wave rejection
#'    and search-back for missed beats (initialised on 2-second training window).
#'
#' R-peak indices are returned relative to the bandpass-filtered signal.
#' Account for `delay` if you need indices aligned to the raw ECG.
#'
#' @references
#' Pan J, Tompkins WJ. A real-time QRS detection algorithm.
#' *IEEE Trans Biomed Eng.* 1985;32(3):230–236.
#' \doi{10.1109/TBME.1985.325532}
#'
#' @seealso [mrpheus::read_philips_physlog()], [mrpheus::compute_hr_signal()],
#'   [mrpheus::compute_hrv_sleep()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec  <- read_philips_physlog("sub-01_physlog.log")
#' qrs  <- detect_qrs(rec$C[, "v1raw"], fs = rec$HDR$sfreq)
#' qrs
#'
#' # Downstream HR signal
#' hr <- compute_hr_signal(qrs)
#' }
detect_qrs <- function(ecg, fs) {

  if (!is.numeric(ecg) || !is.vector(ecg))
    cli::cli_abort("`ecg` must be a numeric vector.")
  if (!is.numeric(fs) || length(fs) != 1L || fs <= 0)
    cli::cli_abort("`fs` must be a single positive number.")

  ecg <- as.double(ecg)

  delay         <- 0
  m_selected_RR <- 0
  mean_RR       <- 0

  # ===========================================================================
  # Stage 1 — Bandpass filtering
  # fs == 200: separate LP (12 Hz) then HP (5 Hz), both Butterworth N=3.
  # Otherwise: single 5–15 Hz bandpass, Butterworth N=3.
  # ===========================================================================
  if (fs == 200) {
    ecg <- ecg - mean(ecg)

    lp_filt <- gsignal::butter(3L, 12 * 2 / fs, type = "low")
    ecg_l   <- as.double(gsignal::filtfilt(lp_filt, ecg))
    ecg_l   <- ecg_l / max(abs(ecg_l))

    hp_filt <- gsignal::butter(3L, 5 * 2 / fs, type = "high")
    ecg_h   <- as.double(gsignal::filtfilt(hp_filt, ecg_l))
    ecg_h   <- ecg_h / max(abs(ecg_h))

  } else {
    bp_filt <- gsignal::butter(3L, c(5, 15) * 2 / fs, type = "pass")
    ecg_h   <- as.double(gsignal::filtfilt(bp_filt, ecg))
    ecg_h   <- ecg_h / max(abs(ecg_h))
  }

  # ===========================================================================
  # Stage 2 — Derivative filter
  # H(z) = (1/8T)(−z^−2 − 2z^−1 + 2z + z^2)
  # For fs != 200: interpolate the 5-point kernel to match the sampling rate.
  # ===========================================================================
  if (fs != 200) {
    int_c   <- (5 - 1) / (fs * (1 / 40))
    b_deriv <- stats::approx(
      x    = 1:5,
      y    = c(1, 2, 0, -2, -1) * (1 / 8) * fs,
      xout = seq(1, 5, by = int_c)
    )$y
  } else {
    b_deriv <- c(1, 2, 0, -2, -1) * (1 / 8) * fs
  }

  ecg_d <- as.double(gsignal::filtfilt(b_deriv, 1, ecg_h))
  ecg_d <- ecg_d / max(ecg_d)

  # ===========================================================================
  # Stage 3 — Squaring
  # ===========================================================================
  ecg_s <- ecg_d^2

  # ===========================================================================
  # Stage 4 — Moving-average integration (150 ms window)
  # Equivalent to MATLAB: conv(ecg_s, ones(1, win) / win)
  #
  # Implemented via cumulative sum rather than stats::convolve(): the two
  # are mathematically identical for a boxcar kernel (verified against the
  # FFT-based version on test signals), but stats::convolve() computes
  # this via FFT, whose performance depends on how well length(ecg_s) + win
  # factors for padding (see stats::nextn()) — for some real recording
  # lengths this degrades from milliseconds to effectively hanging. The
  # cumulative-sum form is unconditionally O(n).
  # ===========================================================================
  win    <- round(0.150 * fs)
  delay  <- delay + win / 2

  n_s     <- length(ecg_s)
  padded  <- c(rep(0, win - 1), ecg_s, rep(0, win - 1))
  cs      <- cumsum(c(0, padded))
  out_len <- n_s + win - 1
  ecg_m   <- (cs[(win + 1):(win + out_len)] - cs[1:out_len]) / win
  ecg_m   <- as.double(ecg_m)

  # ===========================================================================
  # Stage 5 — Adaptive thresholding
  # ===========================================================================

  min_dist <- round(0.2 * fs)
  peaks    <- .find_peaks_min_dist(ecg_m, min_dist)
  pks      <- peaks$pks
  locs     <- peaks$locs
  LLp      <- length(pks)

  # Initialise from 2-second training window
  train_n    <- min(2L * as.integer(fs), length(ecg_m))
  THR_SIG    <- max(ecg_m[1:train_n]) / 3
  THR_NOISE  <- mean(ecg_m[1:train_n]) / 2
  SIG_LEV    <- THR_SIG
  NOISE_LEV  <- THR_NOISE

  train_nh   <- min(2L * as.integer(fs), length(ecg_h))
  THR_SIG1   <- max(ecg_h[1:train_nh]) / 3
  THR_NOISE1 <- mean(ecg_h[1:train_nh]) / 2
  SIG_LEV1   <- THR_SIG1
  NOISE_LEV1 <- THR_NOISE1

  qrs_c       <- double(LLp)
  qrs_i       <- integer(LLp)
  qrs_i_raw   <- integer(LLp)
  qrs_amp_raw <- double(LLp)
  nois_c      <- double(LLp)
  nois_i      <- integer(LLp)

  ecg_h_len   <- length(ecg_h)
  ecg_m_len   <- length(ecg_m)
  Beat_C      <- 0L
  Beat_C1     <- 0L
  Noise_Count <- 0L

  # Defaults retained on loop fall-through (matches MATLAB behaviour)
  y_i      <- 0
  x_i      <- 1L
  ser_back <- 0L
  skip     <- 0L

  for (i in seq_len(LLp)) {

    # -------------------------------------------------------------------------
    # Locate peak in bandpass signal (150 ms window)
    # -------------------------------------------------------------------------
    win_start <- locs[i] - round(0.150 * fs)

    if (win_start >= 1L && locs[i] <= ecg_h_len) {
      sub  <- ecg_h[win_start:locs[i]]
      x_i  <- which.max(sub)
      y_i  <- sub[x_i]
    } else if (i == 1L) {
      sub      <- ecg_h[1L:min(locs[i], ecg_h_len)]
      x_i      <- which.max(sub)
      y_i      <- sub[x_i]
      ser_back <- 1L
    } else if (locs[i] >= ecg_h_len) {
      ws  <- max(1L, win_start)
      sub <- ecg_h[ws:ecg_h_len]
      x_i <- which.max(sub)
      y_i <- sub[x_i]
    }

    # -------------------------------------------------------------------------
    # RR tracking (after 9 confirmed beats)
    # -------------------------------------------------------------------------
    if (Beat_C >= 9L) {
      diffRR  <- diff(qrs_i[(Beat_C - 8L):Beat_C])
      mean_RR <- mean(diffRR)
      comp    <- qrs_i[Beat_C] - qrs_i[Beat_C - 1L]

      if (comp <= 0.92 * mean_RR || comp >= 1.16 * mean_RR) {
        THR_SIG  <- 0.5 * THR_SIG
        THR_SIG1 <- 0.5 * THR_SIG1
      } else {
        m_selected_RR <- mean_RR
      }
    }

    test_m <- if (m_selected_RR > 0) {
      m_selected_RR
    } else if (mean_RR > 0) {
      mean_RR
    } else {
      0
    }

    # -------------------------------------------------------------------------
    # Search-back for missed beat (gap > 1.66 × mean RR)
    # -------------------------------------------------------------------------
    if (test_m > 0 && Beat_C > 0L) {
      if ((locs[i] - qrs_i[Beat_C]) >= round(1.66 * test_m)) {

        ss <- qrs_i[Beat_C] + round(0.200 * fs)
        se <- locs[i] - round(0.200 * fs)

        if (ss < se && se <= ecg_m_len) {
          sub_m     <- ecg_m[ss:se]
          x_t       <- which.max(sub_m)
          pks_temp  <- sub_m[x_t]
          locs_temp <- ss + x_t - 1L

          if (pks_temp > THR_NOISE) {
            Beat_C        <- Beat_C + 1L
            qrs_c[Beat_C] <- pks_temp
            qrs_i[Beat_C] <- locs_temp

            ws2 <- locs_temp - round(0.150 * fs)
            if (ws2 >= 1L && locs_temp <= ecg_h_len) {
              sub_h <- ecg_h[ws2:locs_temp]
            } else {
              sub_h <- ecg_h[max(1L, ws2):min(locs_temp, ecg_h_len)]
            }
            x_i_t <- which.max(sub_h)
            y_i_t <- sub_h[x_i_t]

            if (y_i_t > THR_NOISE1) {
              Beat_C1              <- Beat_C1 + 1L
              qrs_i_raw[Beat_C1]   <- locs_temp - round(0.150 * fs) + (x_i_t - 1L)
              qrs_amp_raw[Beat_C1] <- y_i_t
              SIG_LEV1 <- 0.25 * y_i_t + 0.75 * SIG_LEV1
            }

            SIG_LEV <- 0.25 * pks_temp + 0.75 * SIG_LEV
          }
        }
      }
    }

    # -------------------------------------------------------------------------
    # Classify peak as QRS or noise
    # -------------------------------------------------------------------------
    if (pks[i] >= THR_SIG) {

      # T-wave rejection: within 360 ms of last QRS and shallower slope
      if (Beat_C >= 3L &&
          (locs[i] - qrs_i[Beat_C]) <= round(0.360 * fs)) {

        s1 <- max(1L, locs[i]        - round(0.075 * fs))
        s2 <- max(1L, qrs_i[Beat_C] - round(0.075 * fs))

        slope1 <- mean(diff(ecg_m[s1:locs[i]]))
        slope2 <- mean(diff(ecg_m[s2:qrs_i[Beat_C]]))

        if (abs(slope1) <= abs(0.5 * slope2)) {
          Noise_Count         <- Noise_Count + 1L
          nois_c[Noise_Count] <- pks[i]
          nois_i[Noise_Count] <- locs[i]
          skip                <- 1L
          NOISE_LEV1 <- 0.125 * y_i    + 0.875 * NOISE_LEV1
          NOISE_LEV  <- 0.125 * pks[i] + 0.875 * NOISE_LEV
        } else {
          skip <- 0L
        }
      }

      if (skip == 0L) {
        Beat_C        <- Beat_C + 1L
        qrs_c[Beat_C] <- pks[i]
        qrs_i[Beat_C] <- locs[i]

        if (y_i >= THR_SIG1) {
          Beat_C1 <- Beat_C1 + 1L
          qrs_i_raw[Beat_C1] <- if (ser_back == 1L) {
            x_i
          } else {
            locs[i] - round(0.150 * fs) + (x_i - 1L)
          }
          qrs_amp_raw[Beat_C1] <- y_i
          SIG_LEV1 <- 0.125 * y_i    + 0.875 * SIG_LEV1
        }
        SIG_LEV <- 0.125 * pks[i] + 0.875 * SIG_LEV
      }

    } else if (pks[i] >= THR_NOISE) {
      NOISE_LEV1 <- 0.125 * y_i    + 0.875 * NOISE_LEV1
      NOISE_LEV  <- 0.125 * pks[i] + 0.875 * NOISE_LEV
    } else {
      Noise_Count         <- Noise_Count + 1L
      nois_c[Noise_Count] <- pks[i]
      nois_i[Noise_Count] <- locs[i]
      NOISE_LEV1 <- 0.125 * y_i    + 0.875 * NOISE_LEV1
      NOISE_LEV  <- 0.125 * pks[i] + 0.875 * NOISE_LEV
    }

    # Update adaptive thresholds
    if (NOISE_LEV != 0 || SIG_LEV != 0) {
      THR_SIG   <- NOISE_LEV + 0.25 * abs(SIG_LEV - NOISE_LEV)
      THR_NOISE <- 0.5 * THR_SIG
    }
    if (NOISE_LEV1 != 0 || SIG_LEV1 != 0) {
      THR_SIG1   <- NOISE_LEV1 + 0.25 * abs(SIG_LEV1 - NOISE_LEV1)
      THR_NOISE1 <- 0.5 * THR_SIG1
    }

    skip     <- 0L
    ser_back <- 0L
  }

  # ---- Trim output buffers -------------------------------------------------
  if (Beat_C1 > 0L) {
    qrs_i_raw   <- qrs_i_raw[1:Beat_C1]
    qrs_amp_raw <- qrs_amp_raw[1:Beat_C1]
  } else {
    qrs_i_raw   <- integer(0)
    qrs_amp_raw <- double(0)
  }

  structure(
    list(
      qrs_i   = as.integer(qrs_i_raw),
      qrs_amp = qrs_amp_raw,
      delay   = delay,
      fs      = fs
    ),
    class = "mrpheus_qrs"
  )
}

#' @export
print.mrpheus_qrs <- function(x, ...) {
  n_beats    <- length(x$qrs_i)
  mean_hr    <- if (n_beats >= 2L) {
    round(60 / (mean(diff(x$qrs_i)) / x$fs))
  } else NA_real_

  cli::cli_h1("mrpheus QRS detection")
  cli::cli_inform(c(
    "i" = "R-peaks detected: {n_beats}",
    "i" = "Mean HR:          {mean_hr} bpm",
    "i" = "Filter delay:     {round(x$delay)} samples ({round(x$delay / x$fs, 3)} s)",
    "i" = "Sampling rate:    {x$fs} Hz"
  ))
  invisible(x)
}

# ---------------------------------------------------------------------------
# Internal: local maxima with minimum-distance constraint.
# Matches MATLAB findpeaks(..., 'MinPeakDistance', d): selects highest peaks
# such that no two are within d samples of each other; returns in location order.
# ---------------------------------------------------------------------------
.find_peaks_min_dist <- function(x, min_dist = 1L) {
  n <- length(x)
  if (n < 3L) return(list(pks = double(0), locs = integer(0)))

  is_peak  <- x[2:(n - 1L)] >= x[1:(n - 2L)] & x[2:(n - 1L)] >= x[3:n]
  locs_all <- which(is_peak) + 1L

  if (length(locs_all) == 0L)
    return(list(pks = double(0), locs = integer(0)))

  pks_all <- x[locs_all]
  ord     <- order(pks_all, decreasing = TRUE)
  n_cand  <- length(locs_all)

  # Pre-allocate the accepted-locations buffer ONCE, padded with a sentinel
  # larger than any real location. This lets findInterval() binary-search
  # the full fixed-size buffer every time with no per-candidate subsetting
  # or copying — subsetting a growing vector inside the loop (e.g.
  # `accepted[seq_len(n_accepted)]`) looks fine but still costs O(k) per
  # call, which reintroduces an O(n*k) blow-up over n candidates. This
  # version is genuinely O(n log n).
  sentinel   <- n + min_dist + 1L
  accepted   <- rep(sentinel, n_cand)
  n_accepted <- 0L
  keep       <- logical(n_cand)

  for (k in ord) {
    loc <- locs_all[k]

    pos <- findInterval(loc, accepted)
    too_close <-
      (pos >= 1L && (loc - accepted[pos]) < min_dist) ||
      (pos <  n_accepted && (accepted[pos + 1L] - loc) < min_dist)

    if (!too_close) {
      if (pos < n_accepted) {
        accepted[(pos + 2L):(n_accepted + 1L)] <- accepted[(pos + 1L):n_accepted]
      }
      accepted[pos + 1L] <- loc
      n_accepted <- n_accepted + 1L
      keep[k] <- TRUE
    }
  }

  final <- sort(locs_all[keep])
  list(pks = x[final], locs = final)
}
