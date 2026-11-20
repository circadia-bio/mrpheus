#' Detect slow oscillations
#'
#' Detects slow oscillations (slow waves) in an EEG channel with the algorithm
#' of `yasa.sw_detect()` in YASA (Vallat & Walker, 2021; YASA 0.8.1), and
#' reproduces its output: the same events, at the same sample positions, with
#' the same amplitudes (see Details).
#'
#' The channel is band-passed with an FIR filter (Hamming window,
#' zero phase, as `mne.filter.filter_data()`), applied to the continuous signal,
#' not epoch by epoch. Local minima and maxima of the filtered signal are found
#' (negative peaks with a trough between `amp_neg_uv`, positive peaks between
#' `amp_pos_uv`) and kept if they lie in the epochs given by `stages`. Each
#' negative peak is paired with the closest following positive peak, the
#' peak-to-peak amplitude must lie within `amp_ptp_threshold_uv`, and the
#' durations of the negative and positive phases (taken from the zero crossings
#' around each peak) must lie within `duration_neg_s` and `duration_pos_s`.
#' Every event that shares its start, and then every event that shares its end,
#' with another event is dropped, as in YASA, so two negative peaks paired with
#' the same positive peak are both discarded.
#'
#' @param psg An `mrpheus_psg` object from [mrpheus::prepare_psg()].
#' @param channel Character. EEG channel label. If `NULL` (default), the first
#'   non-bad EEG channel is used.
#' @param stages Integer vector or `NULL`. Epoch indices in which peaks are
#'   accepted (for example the N2 and N3 epochs). If `NULL`, all epochs are
#'   used. As in YASA, only the peaks have to lie in these epochs: an event can
#'   start or end outside them.
#' @param freq_so Numeric vector of length 2. Band-pass of the slow oscillation
#'   band (Hz), as YASA's `freq_sw`. Default `c(0.3, 1.5)`.
#' @param amp_ptp_threshold_uv Numeric vector of length 2. Minimum and maximum
#'   peak-to-peak amplitude (uV), exclusive. Default `c(75, 350)`.
#' @param duration_pos_s Numeric vector of length 2. Minimum and maximum
#'   duration of the positive phase (s), exclusive. Default `c(0.1, 1.0)`.
#' @param duration_neg_s Numeric vector of length 2. Minimum and maximum
#'   duration of the negative phase (s), exclusive. Default `c(0.3, 1.5)`.
#' @param amp_neg_uv Numeric vector of length 2, or `NULL`. Minimum and maximum
#'   absolute amplitude of the negative peak (uV), inclusive. Default
#'   `c(40, 200)`. `NULL` applies no limit.
#' @param amp_pos_uv Numeric vector of length 2, or `NULL`. Minimum and maximum
#'   amplitude of the positive peak (uV), inclusive. Default `c(10, 150)`.
#'   `NULL` applies no limit.
#' @param transition_hz Numeric. Width of both transition bands of the FIR
#'   filter (Hz). The -6 dB points lie half a transition band outside `freq_so`
#'   (0.2 and 1.6 Hz for the defaults). Default `0.2`, as in YASA.
#'
#' @details
#' Amplitudes must be in uV. As YASA, the channel is skipped (with a warning)
#' when the standard deviation of its central 90 % of samples is outside
#' 0.1-1000 uV, which usually means a wrong unit.
#'
#' The result was compared with `yasa.sw_detect()` on the Sleep-EDF night
#' SC4001E0 (8.8 h at 100 Hz, and a 6.5 h version at 256 Hz) with a stage mask
#' that cuts through events: the numbers of events are equal (2,236 and 1,000),
#' the event times are identical, and amplitudes agree to about 1e-12 uV. The
#' test suite checks the same agreement against stored YASA output on a
#' synthetic EEG at 100 and 256 Hz, and on a crafted signal that reaches the
#' edge cases of YASA's rules
#' (`tests/testthat/test-yasa-slow-oscillation-parity.R`).
#'
#' Not implemented: YASA's slow oscillation - spindle coupling and its
#' `remove_outliers` option. Where YASA's code has quirks they are reproduced:
#' a peak before the first zero crossing takes the last zero crossing of the
#' recording as its previous crossing (the event is then rejected by the
#' duration criteria).
#'
#' @return A tibble with one row per detected slow oscillation:
#' \describe{
#'   \item{epoch}{Integer. Epoch of the negative peak.}
#'   \item{start_s}{Numeric. Onset (s) relative to the start of `epoch`. Can be
#'     negative, or the offset can exceed the epoch length, when the event
#'     crosses an epoch boundary.}
#'   \item{end_s}{Numeric. Offset (s) relative to the start of `epoch`.}
#'   \item{duration_s}{Numeric. Duration of the slow oscillation (s).}
#'   \item{neg_peak_uv}{Numeric. Value of the negative peak in the filtered
#'     signal (uV, negative).}
#'   \item{pos_peak_uv}{Numeric. Value of the positive peak (uV).}
#'   \item{ptp_uv}{Numeric. Peak-to-peak amplitude (uV).}
#'   \item{channel}{Character.}
#'   \item{neg_peak_s}{Numeric. Time of the negative peak (s from the start of
#'     the first epoch).}
#'   \item{mid_crossing_s}{Numeric. Time of the zero crossing between the
#'     negative and the positive peak (s from the start of the first epoch).}
#'   \item{pos_peak_s}{Numeric. Time of the positive peak (s from the start of
#'     the first epoch).}
#'   \item{slope_uv_s}{Numeric. Slope between the negative peak and the mid
#'     crossing (uV/s).}
#' }
#'
#' @references
#' Vallat, R., & Walker, M. P. (2021). An open-source, high-performance tool
#' for automated sleep staging. *eLife*, 10, e70092.
#' \doi{10.7554/eLife.70092}
#'
#' @export
compute_slow_oscillations <- function(psg,
                                      channel              = NULL,
                                      stages               = NULL,
                                      freq_so              = c(0.3, 1.5),
                                      amp_ptp_threshold_uv = c(75, 350),
                                      duration_pos_s       = c(0.1, 1.0),
                                      duration_neg_s       = c(0.3, 1.5),
                                      amp_neg_uv           = c(40, 200),
                                      amp_pos_uv           = c(10, 150),
                                      transition_hz        = 0.2) {
  stopifnot(inherits(psg, "mrpheus_psg"))

  if (is.null(channel)) {
    channel <- psg$channel_map$label[
      psg$channel_map$type == "EEG" & !psg$channel_map$bad
    ][1]
  }

  sr       <- psg$channel_map$sample_rate[psg$channel_map$label == channel]
  n_epochs <- length(psg$epochs)
  spe      <- as.integer(round(psg$epoch_s * sr))
  epoch_s  <- spe / sr

  nyq <- sr / 2
  lo  <- freq_so[1] - transition_hz
  hi  <- freq_so[2] + transition_hz
  if (!(transition_hz > 0) || lo < 0 || hi > nyq || freq_so[1] >= freq_so[2]) {
    cli::cli_abort(c(
      "Invalid slow oscillation filter band.",
      "i" = "The FIR filter needs 0 <= {.code freq_so[1] - transition_hz} ({lo} Hz) and {.code freq_so[2] + transition_hz} ({hi} Hz) <= {nyq} Hz, with {.code transition_hz} > 0 and {.code freq_so[1] < freq_so[2]}."
    ))
  }

  # the continuous channel signal, and the sample mask of the accepted epochs
  sig <- unlist(lapply(psg$epochs, function(e) e[[channel]]), use.names = FALSE)
  if (length(sig) != n_epochs * spe)
    cli::cli_abort("Channel {.val {channel}} is not present, with {spe} samples, in every epoch.")
  keep_epoch <- if (is.null(stages)) rep(TRUE, n_epochs) else seq_len(n_epochs) %in% stages
  mask       <- rep(keep_epoch, each = spe)

  empty <- tibble::tibble(
    epoch = integer(), start_s = numeric(), end_s = numeric(), duration_s = numeric(),
    neg_peak_uv = numeric(), pos_peak_uv = numeric(), ptp_uv = numeric(),
    channel = character(), neg_peak_s = numeric(), mid_crossing_s = numeric(),
    pos_peak_s = numeric(), slope_uv_s = numeric()
  )

  # YASA skips a channel whose trimmed standard deviation is not 0.1-1000 uV
  trim_sd <- stats::sd(.trim_both(sig, 0.05))
  if (!(trim_sd > 0.1 && trim_sd < 1e3)) {
    cli::cli_warn(c(
      "Wrong data amplitude for channel {.val {channel}} (trimmed SD = {signif(trim_sd, 3)}).",
      "i" = "The unit of the data must be uV. The channel is skipped."
    ))
    return(empty)
  }

  ev <- .yasa_sw_detect(
    sig, sr, mask, freq_so, duration_neg_s, duration_pos_s,
    amp_neg_uv %||% c(-Inf, Inf), amp_pos_uv %||% c(-Inf, Inf),
    amp_ptp_threshold_uv, transition_hz
  )
  if (is.null(ev)) {
    cli::cli_alert_success("Detected 0 slow oscillations in channel {.val {channel}}.")
    return(empty)
  }

  epoch <- as.integer(ev$NegIdx %/% spe) + 1L
  out <- tibble::tibble(
    epoch          = epoch,
    start_s        = ev$Start - (epoch - 1L) * epoch_s,
    end_s          = ev$End   - (epoch - 1L) * epoch_s,
    duration_s     = ev$Duration,
    neg_peak_uv    = ev$ValNegPeak,
    pos_peak_uv    = ev$ValPosPeak,
    ptp_uv         = ev$PTP,
    channel        = channel,
    neg_peak_s     = ev$NegPeak,
    mid_crossing_s = ev$MidCrossing,
    pos_peak_s     = ev$PosPeak,
    slope_uv_s     = ev$Slope
  )
  cli::cli_alert_success(
    "Detected {nrow(out)} slow oscillations in channel {.val {channel}}."
  )
  out
}

# ── Internals ─────────────────────────────────────────────────────────────────

# scipy.stats.trimboth(x, proportiontocut): sorted values without the lowest
# and the highest `cut` share.
.trim_both <- function(x, cut) {
  n   <- length(x)
  low <- as.integer(cut * n)
  if (low == 0L) return(x)
  sort(x)[seq.int(low + 1L, n - low)]
}

# scipy.signal.find_peaks(x, height = c(hmin, hmax)): 0-based indices of the
# local maxima (a plateau counts once, at its middle) whose height lies within
# the limits (inclusive).
.sw_find_peaks <- function(x, height = NULL) {
  n <- length(x)
  if (n < 3L) return(integer(0))
  d    <- diff(x)                                    # d[i] = x[i + 1] - x[i]
  up   <- d[-length(d)] > 0                          # x[i] > x[i - 1], i = 2..n-1
  dn   <- d[-1L] < 0                                 # x[i + 1] < x[i]
  flat <- d[-1L] == 0
  i    <- which(up & dn) + 1L                        # 1-based indices of strict peaks
  fl   <- which(up & flat) + 1L                      # starts of plateaus
  if (length(fl)) {
    mids <- integer(0)
    for (s in fl) {
      j <- s
      while (j < n && x[j + 1L] == x[s]) j <- j + 1L # last sample of the plateau
      if (j < n && x[j + 1L] < x[s])                 # it drops afterwards: a peak
        mids <- c(mids, (s - 1L + j - 1L) %/% 2L + 1L)
    }
    i <- sort(c(i, mids))
  }
  if (!is.null(height) && length(i))
    i <- i[x[i] >= height[1L] & x[i] <= height[2L]]
  i - 1L
}

# np.round(x, 4) (multiply, round half to even, divide)
.np_round4 <- function(x) round(x * 1e4) / 1e4

# Port of the single-channel detection in yasa.sw_detect() (YASA 0.8.1):
# band-pass, then the events of the filtered signal.
.yasa_sw_detect <- function(x, sf, mask, freq_sw, dur_neg, dur_pos,
                            amp_neg, amp_pos, amp_ptp, transition) {
  h  <- .mne_fir_bandpass_design(sf, freq_sw[1L], freq_sw[2L], transition, transition)
  xf <- .mne_fir_apply(x, h)
  .sw_events(xf, sf, mask, dur_neg, dur_pos, amp_neg, amp_pos, amp_ptp)
}

# The detection proper, on the band-passed signal `xf`. `mask` is a logical
# vector with one entry per sample. Returns a data frame with YASA's columns
# (times in seconds from the first sample; `NegIdx` is the 0-based sample index
# of the negative peak), or NULL when nothing is found.
.sw_events <- function(xf, sf, mask, dur_neg, dur_pos, amp_neg, amp_pos, amp_ptp) {
  idx_neg <- .sw_find_peaks(-xf, amp_neg)
  idx_pos <- .sw_find_peaks(xf,  amp_pos)
  idx_neg <- idx_neg[mask[idx_neg + 1L]]
  idx_pos <- idx_pos[mask[idx_pos + 1L]]
  if (!length(idx_neg) || !length(idx_pos)) return(NULL)

  # the last peak must be a positive one
  if (idx_pos[length(idx_pos)] < idx_neg[length(idx_neg)])
    idx_pos <- c(idx_pos, idx_neg[length(idx_neg)] + 1L)

  # closest following positive peak of each negative peak (np.searchsorted, left)
  pk      <- findInterval(idx_neg, idx_pos, left.open = TRUE) + 1L
  idx_pos <- idx_pos[pk]

  sw_pt  <- abs(xf[idx_neg + 1L])
  sw_ptp <- sw_pt + xf[idx_pos + 1L]
  good   <- sw_ptp > amp_ptp[1L] & sw_ptp < amp_ptp[2L]
  if (!any(good)) return(NULL)
  sw_ptp  <- sw_ptp[good];  sw_pt  <- sw_pt[good]
  idx_neg <- idx_neg[good]; idx_pos <- idx_pos[good]

  # zero crossings of the filtered signal (x > 0 changes between i and i + 1)
  pos <- xf > 0
  zc  <- which(pos[-length(pos)] != pos[-1L]) - 1L
  last_peak <- max(idx_pos[length(idx_pos)], idx_neg[length(idx_neg)])
  if (!length(zc) || zc[length(zc)] < last_peak) zc <- c(zc, last_peak)

  at_zc <- function(k) zc[ifelse(k == 0L, length(zc), k)]   # zc[k - 1]; Python wraps at 0
  neg_sorted <- findInterval(idx_neg, zc, left.open = TRUE)
  prev_neg   <- at_zc(neg_sorted)   - idx_neg
  foll_neg   <- zc[neg_sorted + 1L] - idx_neg
  pos_sorted <- findInterval(idx_pos, zc, left.open = TRUE)
  prev_pos   <- at_zc(pos_sorted)   - idx_pos
  foll_pos   <- zc[pos_sorted + 1L] - idx_pos

  neg_dur <- (abs(prev_neg) + foll_neg) / sf
  pos_dur <- (abs(prev_pos) + foll_pos) / sf

  sw_start <- (idx_neg + prev_neg) / sf
  sw_end   <- (idx_pos + foll_pos) / sf
  sw_dur   <- .np_round4(sw_end - sw_start)
  sw_dur2  <- .np_round4(pos_dur + neg_dur)
  sw_mid   <- (idx_neg + foll_neg) / sf
  sw_neg_t <- idx_neg / sf
  sw_slope <- sw_pt / (sw_mid - sw_neg_t)

  ok <- prev_neg != 0 & foll_neg != 0 & prev_pos != 0 & foll_pos != 0 &
    sw_dur == sw_dur2 &
    sw_dur <= dur_neg[2L] + dur_pos[2L] & sw_dur >= dur_neg[1L] + dur_pos[1L] &
    neg_dur > dur_neg[1L] & neg_dur < dur_neg[2L] &
    pos_dur > dur_pos[1L] & pos_dur < dur_pos[2L] &
    sw_mid > sw_start & sw_mid < sw_end &
    is.finite(sw_slope) & sw_slope > 0
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(NULL)

  # YASA drops every event that shares its start, and then every event that
  # shares its end, with another event (drop_duplicates(keep = False)): two
  # negative peaks paired with the same positive peak are both discarded.
  keep <- which(ok)
  dup  <- function(v) duplicated(v) | duplicated(v, fromLast = TRUE)
  keep <- keep[!dup(sw_start[keep])]
  keep <- keep[!dup(sw_end[keep])]
  if (!length(keep)) return(NULL)

  data.frame(
    Start = sw_start[keep], NegPeak = sw_neg_t[keep], MidCrossing = sw_mid[keep],
    PosPeak = idx_pos[keep] / sf, End = sw_end[keep], Duration = sw_dur[keep],
    ValNegPeak = xf[idx_neg[keep] + 1L], ValPosPeak = xf[idx_pos[keep] + 1L],
    PTP = sw_ptp[keep], Slope = sw_slope[keep], Frequency = 1 / sw_dur[keep],
    NegIdx = idx_neg[keep]
  )
}
