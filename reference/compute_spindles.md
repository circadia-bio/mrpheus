# Detect sleep spindles

Detects sleep spindles in an EEG channel with the algorithm of
`yasa.spindles_detect()` in YASA (Vallat & Walker, 2021; YASA 0.8.1),
which follows Lacourse et al. (2019), and reproduces its output: the
same events, at the same sample positions, with the same properties (see
Details).

## Usage

``` r
compute_spindles(
  psg,
  channel = NULL,
  stages = NULL,
  freq_spindle = c(12, 15),
  freq_broad = c(1, 30),
  duration_s = c(0.5, 2),
  min_distance_ms = 500,
  thresh_rel_pow = 0.2,
  thresh_corr = 0.65,
  thresh_rms = 1.5
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

  Integer vector or `NULL`. Epoch indices in which spindles may occur
  (for example the N1, N2 and N3 epochs). If `NULL`, all epochs are
  used. As in YASA, the RMS threshold is computed within these epochs
  and samples outside them cannot start a spindle, but an event can
  extend a little beyond them (the smoothing spreads over 0.1 s).

- freq_spindle:

  Numeric vector of length 2. Spindle band (Hz), as YASA's `freq_sp`.
  Default `c(12, 15)`.

- freq_broad:

  Numeric vector of length 2. Broadband filter (Hz). Its transition
  bands are chosen as MNE does for `filter_data()`: 25 % of the cutoff,
  at least 2 Hz, and no more than the cutoff for the high-pass. Default
  `c(1, 30)`.

- duration_s:

  Numeric vector of length 2. Minimum and maximum duration of a spindle
  (s), exclusive. Default `c(0.5, 2)`.

- min_distance_ms:

  Numeric or `NULL`. Spindles separated by less than this (ms) are
  merged. Default `500`. `NULL` or `0` switches merging off.

- thresh_rel_pow:

  Numeric or `NULL`. Minimum relative sigma power. Default `0.20`.
  `NULL` switches the criterion off.

- thresh_corr:

  Numeric or `NULL`. Minimum moving correlation. Default `0.65`. `NULL`
  switches the criterion off.

- thresh_rms:

  Numeric or `NULL`. Number of trimmed standard deviations above the
  mean for the moving RMS. Default `1.5`. `NULL` switches the criterion
  off. At least one criterion must be given.

## Value

A tibble with one row per detected spindle:

- epoch:

  Integer. Epoch of the first sample of the spindle.

- start_s:

  Numeric. Onset (s) relative to the start of `epoch`.

- end_s:

  Numeric. Offset (s) relative to the start of `epoch`; can exceed the
  epoch length when the spindle crosses an epoch boundary.

- duration_s:

  Numeric. Duration (s).

- peak_freq_hz:

  Numeric. Median instantaneous frequency (Hz) of the sigma-filtered
  signal within the spindle (YASA's `Frequency`).

- rms_uv:

  Numeric. Root mean square (uV) of the detrended broadband signal
  within the spindle (YASA's `RMS`).

- channel:

  Character.

- peak_s:

  Numeric. Time of the most prominent peak (s from the start of the
  first epoch; YASA's `Peak`).

- amplitude_uv:

  Numeric. Peak-to-peak amplitude of the detrended broadband signal (uV;
  YASA's `Amplitude`).

- amp_filtered_uv:

  Numeric. Peak-to-peak amplitude of the sigma-filtered signal (uV;
  YASA's `AmpFiltered`).

- abs_power:

  Numeric. Median log10 instantaneous power of the sigma signal (YASA's
  `AbsPower`).

- rel_power:

  Numeric. Median relative sigma power (YASA's `RelPower`).

- oscillations:

  Integer. Number of oscillations: peaks at least 60 ms apart.

- symmetry:

  Numeric. Position of the most prominent peak within the spindle, from
  0 to 1 (YASA's `Symmetry`).

## Details

The channel is band-passed twice with FIR filters (Hamming window, zero
phase, as `mne.filter.filter_data()`) applied to the continuous signal:
a broadband filter (`freq_broad`, default 1-30 Hz) and a sigma filter
(`freq_spindle`, default 12-15 Hz, transition bands of 1.5 Hz). Three
criteria are evaluated at every sample: the relative power of the sigma
band in the broadband signal (a short-time Fourier transform with 2 s
windows every 0.2 s, linearly interpolated), the moving correlation
between the sigma and the broadband signal (0.3 s windows) and the
moving RMS of the sigma signal (0.3 s windows; the threshold is the mean
plus `thresh_rms` trimmed standard deviations within the accepted
epochs, at most 10 uV). The number of criteria met at each sample is
smoothed over 0.1 s and a sample belongs to a spindle where at least all
but one of the criteria are met. Events closer than `min_distance_ms`
are merged, and those whose duration is within `duration_s` are kept.

Amplitudes must be in uV. As YASA, the channel is skipped (with a
warning) when the standard deviation of its central 90 % of samples is
outside 0.1-1000 uV, which usually means a wrong unit.

The result was compared with `yasa.spindles_detect()` on the Sleep-EDF
night SC4001E0 (8.8 h at 100 Hz, and a 6.5 h version at 256 Hz) with a
stage mask that cuts through events: the numbers of events are equal (73
and 78), the start and end times are identical, and every property
agrees to about 1e-11. The test suite checks the same agreement against
stored YASA output on a synthetic EEG with spindles at 100 and 256 Hz
(`tests/testthat/test-yasa-spindle-parity.R`).

Not implemented: YASA's `multi_only` and `remove_outliers` options, and
its spindle - slow oscillation coupling.

## References

Lacourse, K., Delfrate, J., Beaudry, J., Peppard, P., & Warby, S. C.
(2019). A sleep spindle detection algorithm that emulates human expert
spindle scoring. *Journal of Neuroscience Methods*, 316, 3-11.

Vallat, R., & Walker, M. P. (2021). An open-source, high-performance
tool for automated sleep staging. *eLife*, 10, e70092.
[doi:10.7554/eLife.70092](https://doi.org/10.7554/eLife.70092)
