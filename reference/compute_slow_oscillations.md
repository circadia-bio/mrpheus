# Detect slow oscillations

Detects slow oscillations (SOs) in EEG channels using a zero-crossing
approach, following the algorithm described in Molle et al. (2002). The
defaults (band, duration, trough, peak and peak-to-peak limits) are
those of `yasa.sw_detect()` in YASA (Vallat & Walker, 2021), so that
densities are comparable in magnitude. The detection itself is not
identical to YASA's: this function pairs zero crossings into a negative
and a following positive half-wave, whereas YASA pairs each negative
peak with the nearest following positive peak, so event counts and
timings differ. Epochs are filtered one at a time.

## Usage

``` r
compute_slow_oscillations(
  psg,
  channel = NULL,
  stages = NULL,
  freq_so = c(0.3, 1.5),
  amp_ptp_threshold_uv = c(75, 350),
  duration_pos_s = c(0.1, 1),
  duration_neg_s = c(0.3, 1.5),
  amp_neg_uv = c(40, 200),
  amp_pos_uv = c(10, 150),
  transition_hz = 0.2
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- channel:

  Character. EEG channel label. If `NULL` (default), the first non-bad
  EEG channel is used.

- stages:

  Integer vector or `NULL`. Epoch indices restricted to N2/N3. If
  `NULL`, detection runs across all epochs.

- freq_so:

  Numeric vector of length 2. SO frequency band (Hz), as in YASA's
  `freq_sw`. Default `c(0.3, 1.5)`. See `transition_hz`.

- amp_ptp_threshold_uv:

  Numeric vector of length 2. Min and max acceptable peak-to-peak
  amplitude (uV). Default `c(75, 350)`.

- duration_pos_s:

  Numeric vector of length 2. Min and max positive half- wave duration
  (s). Default `c(0.1, 1.0)`.

- duration_neg_s:

  Numeric vector of length 2. Min and max negative half- wave duration
  (s). Default `c(0.3, 1.5)`.

- amp_neg_uv:

  Numeric vector of length 2, or `NULL`. Min and max absolute amplitude
  of the negative peak (uV). Default `c(40, 200)`. `NULL` applies no
  limit.

- amp_pos_uv:

  Numeric vector of length 2, or `NULL`. Min and max amplitude of the
  positive peak (uV). Default `c(10, 150)`. `NULL` applies no limit.

- transition_hz:

  Numeric. Width of the filter transition band (Hz). YASA band-passes
  with an FIR filter whose -6 dB points lie half a transition band
  outside `freq_so` (0.2 and 1.6 Hz for the defaults). Here the
  band-pass is a zero-phase Butterworth high-pass / low-pass cascade
  whose -6 dB points are placed at `freq_so[1] - transition_hz / 2` and
  `freq_so[2] + transition_hz / 2`. Default `0.2`.

## Value

A tibble with one row per detected slow oscillation:

- epoch:

  Integer.

- start_s:

  Numeric. Onset (s) relative to epoch start.

- end_s:

  Numeric. Offset (s) relative to epoch start.

- duration_s:

  Numeric.

- neg_peak_uv:

  Numeric. Negative peak amplitude (uV).

- pos_peak_uv:

  Numeric. Positive peak amplitude (uV).

- ptp_uv:

  Numeric. Peak-to-peak amplitude (uV).

- channel:

  Character.

## Details

Amplitudes are assumed to be in uV. The high-pass is a 4th-order
Butterworth in transfer-function form; its coefficients become
ill-conditioned at very low cutoffs and high sample rates (rcond about
2e-14 at 1000 Hz for a 0.2 Hz cutoff, against a failure limit near
2e-16), so very high sample rates should be downsampled first.

## References

Molle, M., Marshall, L., Gais, S., & Born, J. (2002). Grouping of
spindle activity during slow oscillations in human non-rapid eye
movement sleep. *Journal of Neuroscience*, 22(24), 10941-10947.

Vallat, R., & Walker, M. P. (2021). An open-source, high-performance
tool for automated sleep staging. *eLife*, 10, e70092.
[doi:10.7554/eLife.70092](https://doi.org/10.7554/eLife.70092)
