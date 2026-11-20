# Detect slow oscillations

Detects slow oscillations (slow waves) in an EEG channel with the
algorithm of `yasa.sw_detect()` in YASA (Vallat & Walker, 2021; YASA
0.8.1), and reproduces its output: the same events, at the same sample
positions, with the same amplitudes (see Details).

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

  Integer vector or `NULL`. Epoch indices in which peaks are accepted
  (for example the N2 and N3 epochs). If `NULL`, all epochs are used. As
  in YASA, only the peaks have to lie in these epochs: an event can
  start or end outside them.

- freq_so:

  Numeric vector of length 2. Band-pass of the slow oscillation band
  (Hz), as YASA's `freq_sw`. Default `c(0.3, 1.5)`.

- amp_ptp_threshold_uv:

  Numeric vector of length 2. Minimum and maximum peak-to-peak amplitude
  (uV), exclusive. Default `c(75, 350)`.

- duration_pos_s:

  Numeric vector of length 2. Minimum and maximum duration of the
  positive phase (s), exclusive. Default `c(0.1, 1.0)`.

- duration_neg_s:

  Numeric vector of length 2. Minimum and maximum duration of the
  negative phase (s), exclusive. Default `c(0.3, 1.5)`.

- amp_neg_uv:

  Numeric vector of length 2, or `NULL`. Minimum and maximum absolute
  amplitude of the negative peak (uV), inclusive. Default `c(40, 200)`.
  `NULL` applies no limit.

- amp_pos_uv:

  Numeric vector of length 2, or `NULL`. Minimum and maximum amplitude
  of the positive peak (uV), inclusive. Default `c(10, 150)`. `NULL`
  applies no limit.

- transition_hz:

  Numeric. Width of both transition bands of the FIR filter (Hz). The -6
  dB points lie half a transition band outside `freq_so` (0.2 and 1.6 Hz
  for the defaults). Default `0.2`, as in YASA.

## Value

A tibble with one row per detected slow oscillation:

- epoch:

  Integer. Epoch of the negative peak.

- start_s:

  Numeric. Onset (s) relative to the start of `epoch`. Can be negative,
  or the offset can exceed the epoch length, when the event crosses an
  epoch boundary.

- end_s:

  Numeric. Offset (s) relative to the start of `epoch`.

- duration_s:

  Numeric. Duration of the slow oscillation (s).

- neg_peak_uv:

  Numeric. Value of the negative peak in the filtered signal (uV,
  negative).

- pos_peak_uv:

  Numeric. Value of the positive peak (uV).

- ptp_uv:

  Numeric. Peak-to-peak amplitude (uV).

- channel:

  Character.

- neg_peak_s:

  Numeric. Time of the negative peak (s from the start of the first
  epoch).

- mid_crossing_s:

  Numeric. Time of the zero crossing between the negative and the
  positive peak (s from the start of the first epoch).

- pos_peak_s:

  Numeric. Time of the positive peak (s from the start of the first
  epoch).

- slope_uv_s:

  Numeric. Slope between the negative peak and the mid crossing (uV/s).

## Details

The channel is band-passed with an FIR filter (Hamming window, zero
phase, as `mne.filter.filter_data()`), applied to the continuous signal,
not epoch by epoch. Local minima and maxima of the filtered signal are
found (negative peaks with a trough between `amp_neg_uv`, positive peaks
between `amp_pos_uv`) and kept if they lie in the epochs given by
`stages`. Each negative peak is paired with the closest following
positive peak, the peak-to-peak amplitude must lie within
`amp_ptp_threshold_uv`, and the durations of the negative and positive
phases (taken from the zero crossings around each peak) must lie within
`duration_neg_s` and `duration_pos_s`. Every event that shares its
start, and then every event that shares its end, with another event is
dropped, as in YASA, so two negative peaks paired with the same positive
peak are both discarded.

Amplitudes must be in uV. As YASA, the channel is skipped (with a
warning) when the standard deviation of its central 90 % of samples is
outside 0.1-1000 uV, which usually means a wrong unit.

The result was compared with `yasa.sw_detect()` on the Sleep-EDF night
SC4001E0 (8.8 h at 100 Hz, and a 6.5 h version at 256 Hz) with a stage
mask that cuts through events: the numbers of events are equal (2,236
and 1,000), the event times are identical, and amplitudes agree to about
1e-12 uV. The test suite checks the same agreement against stored YASA
output on a synthetic EEG at 100 and 256 Hz, and on a crafted signal
that reaches the edge cases of YASA's rules
(`tests/testthat/test-yasa-slow-oscillation-parity.R`).

Not implemented: YASA's slow oscillation - spindle coupling and its
`remove_outliers` option. Where YASA's code has quirks they are
reproduced: a peak before the first zero crossing takes the last zero
crossing of the recording as its previous crossing (the event is then
rejected by the duration criteria).

## References

Vallat, R., & Walker, M. P. (2021). An open-source, high-performance
tool for automated sleep staging. *eLife*, 10, e70092.
[doi:10.7554/eLife.70092](https://doi.org/10.7554/eLife.70092)
