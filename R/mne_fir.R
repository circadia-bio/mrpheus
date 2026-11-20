# R/mne_fir.R
#
# MNE-equivalent FIR band-pass, for any sampling rate:
#   mne.filter.filter_data(x, sfreq, l_freq, h_freq, method = "fir",
#     l_trans_bandwidth, h_trans_bandwidth, filter_length = "auto",
#     phase = "zero", fir_window = "hamming", fir_design = "firwin",
#     pad = "reflect_limited")
#
# The design (.mne_fir_bandpass_design) equals mne.filter.create_filter to about
# 1e-16 (checked at 100, 128, 200, 256 and 500 Hz); the filtered signal
# (.mne_fir_apply) equals filter_data to about 5e-13 uV on real EEG. Used by
# compute_slow_oscillations(). The staging code filters with stored 100 Hz
# coefficients instead (inst/filters/).

# np.sinc
.sinc <- function(x) {
  out <- rep(1, length(x))
  nz  <- x != 0
  out[nz] <- sin(pi * x[nz]) / (pi * x[nz])
  out
}

# scipy.signal.firwin(numtaps, cutoff, window = "hamming", pass_zero = TRUE):
# a low-pass with unit gain at DC. `cutoff` is a fraction of the Nyquist frequency.
.firwin_lowpass <- function(numtaps, cutoff) {
  m   <- seq_len(numtaps) - 1 - 0.5 * (numtaps - 1)
  h   <- cutoff * .sinc(cutoff * m)
  win <- 0.54 - 0.46 * cos(2 * pi * (seq_len(numtaps) - 1) / (numtaps - 1))   # symmetric
  h   <- h * win
  h / sum(h)
}

# mne.filter.create_filter(None, sfreq, l_freq, h_freq, l_trans_bandwidth,
# h_trans_bandwidth, filter_length = "auto", method = "fir", phase = "zero",
# fir_window = "hamming", fir_design = "firwin") for a band-pass.
.mne_fir_bandpass_design <- function(sfreq, l_freq, h_freq, l_trans, h_trans) {
  len_factor <- 3.3                                  # Hamming window
  nyq   <- sfreq / 2
  f_p1  <- l_freq;  f_p2 <- h_freq
  f_s1  <- l_freq - l_trans
  f_s2  <- h_freq + h_trans
  if (f_s1 < 0 || f_s2 > nyq || l_freq >= h_freq)
    stop("invalid band-pass design")

  # filter_length = "auto": len_factor / min(transition) seconds, made odd
  N <- as.integer(round(len_factor / min(l_trans, h_trans) * sfreq))
  N <- N + (N - 1L) %% 2L

  freq <- c(f_s1, f_p1, f_p2, f_s2) / nyq
  gain <- c(0, 1, 1, 0)
  if (f_s2 != nyq) { freq <- c(freq, 1); gain <- c(gain, 0) }
  if (f_s1 != 0)   { freq <- c(0, freq); gain <- c(0, gain) }

  h <- numeric(N)
  prev_freq <- freq[length(freq)]
  prev_gain <- gain[length(gain)]
  if (prev_gain == 1) h[(N %/% 2L) + 1L] <- 1
  rev_f <- rev(freq)[-1L]
  rev_g <- rev(gain)[-1L]
  for (k in seq_along(rev_f)) {
    this_freq <- rev_f[k]; this_gain <- rev_g[k]
    if (this_gain != prev_gain) {
      transition <- (prev_freq - this_freq) / 2
      this_N <- as.integer(round(len_factor / transition))
      this_N <- this_N + 1L - this_N %% 2L
      if (this_N > N) stop("requested filter length too short for the transition band")
      this_h <- .firwin_lowpass(this_N, (prev_freq + this_freq) / 2)
      offset <- (N - this_N) %/% 2L
      rng    <- seq.int(offset + 1L, N - offset)
      if (this_gain == 0) h[rng] <- h[rng] - this_h else h[rng] <- h[rng] + this_h
    }
    prev_gain <- this_gain
    prev_freq <- this_freq
  }
  h
}

# Full linear convolution of x with h by overlap-add FFT.
.fft_conv_full <- function(x, h) {
  M <- length(h); N <- length(x)
  n_fft <- 2^ceiling(log2(max(8 * M, 1024)))
  L <- n_fft - M + 1L
  H <- stats::fft(c(h, numeric(n_fft - M)))
  out <- numeric(N + M - 1L)
  for (s in seq.int(1L, N, by = L)) {
    e   <- min(s + L - 1L, N)
    seg <- x[s:e]
    y   <- Re(stats::fft(stats::fft(c(seg, numeric(n_fft - length(seg)))) * H,
                         inverse = TRUE)) / n_fft
    ix  <- s:(s + length(seg) + M - 2L)
    out[ix] <- out[ix] + y[seq_along(ix)]
  }
  out
}

# Zero-phase application of a symmetric (odd-length) FIR, as MNE's
# _overlap_add_filter with pad = "reflect_limited".
.mne_fir_apply <- function(x, h) {
  M      <- length(h)
  N      <- length(x)
  n_edge <- max(min(M, N) - 1L, 0L)
  xp     <- .mne_smart_pad(x, c(n_edge, n_edge))
  delay  <- (M - 1L) %/% 2L
  y      <- .fft_conv_full(xp, h)
  y      <- y[seq.int(delay + 1L, delay + length(xp))]
  y[seq.int(n_edge + 1L, n_edge + N)]
}
