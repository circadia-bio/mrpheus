# Detect slow oscillations

Detects slow oscillations (SOs) in EEG channels using a zero-crossing
approach in the delta band (0.5–2 Hz), following the algorithm described
in Mölle et al. (2002) and implemented in YASA (Vallat & Walker, 2021).

## Usage

``` r
compute_slow_oscillations(
  psg,
  channel = NULL,
  stages = NULL,
  freq_so = c(0.5, 2),
  amp_ptp_threshold_uv = c(75, 500),
  duration_pos_s = c(0.1, 1),
  duration_neg_s = c(0.1, 1.5)
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

  Numeric vector of length 2. SO frequency band (Hz). Default
  `c(0.5, 2)`.

- amp_ptp_threshold_uv:

  Numeric vector of length 2. Min and max acceptable peak-to-peak
  amplitude (µV). Default `c(75, 500)`.

- duration_pos_s:

  Numeric vector of length 2. Min and max positive half- wave duration
  (s). Default `c(0.1, 1.0)`.

- duration_neg_s:

  Numeric vector of length 2. Min and max negative half- wave duration
  (s). Default `c(0.1, 1.5)`.

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

  Numeric. Negative peak amplitude (µV).

- pos_peak_uv:

  Numeric. Positive peak amplitude (µV).

- ptp_uv:

  Numeric. Peak-to-peak amplitude (µV).

- channel:

  Character.

## References

Mölle, M., Marshall, L., Gais, S., & Born, J. (2002). Grouping of
spindle activity during slow oscillations in human non-rapid eye
movement sleep. *Journal of Neuroscience*, 22(24), 10941–10947.

Vallat, R., & Walker, M. P. (2021). An open-source, high-performance
tool for automated sleep staging. *eLife*, 10, e70092.
[doi:10.7554/eLife.70092](https://doi.org/10.7554/eLife.70092)
