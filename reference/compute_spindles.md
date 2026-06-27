# Detect sleep spindles

Detects sleep spindles in EEG channels using a band-pass / RMS envelope
approach, closely following the algorithm in Lacourse et al. (2019) and
implemented in YASA (Vallat & Walker, 2021). Spindles are identified as
transient bursts of 11–16 Hz activity during NREM sleep.

## Usage

``` r
compute_spindles(
  psg,
  channel = NULL,
  stages = NULL,
  freq_spindle = c(11, 16),
  rms_window_s = 0.3,
  min_duration_s = 0.5,
  max_duration_s = 3,
  threshold_sd = 1.5
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

  Integer vector or `NULL`. Epoch indices (1-based) restricted to NREM
  sleep. If `NULL`, spindles are detected across all epochs.

- freq_spindle:

  Numeric vector of length 2. Spindle frequency band (Hz). Default
  `c(11, 16)`.

- rms_window_s:

  Numeric. Moving RMS window length in seconds. Default `0.3`.

- min_duration_s:

  Numeric. Minimum spindle duration in seconds. Default `0.5`.

- max_duration_s:

  Numeric. Maximum spindle duration in seconds. Default `3.0`.

- threshold_sd:

  Numeric. RMS threshold as a multiple of the channel SD (within sigma
  band). Default `1.5`.

## Value

A tibble with one row per detected spindle:

- epoch:

  Integer. Epoch in which the spindle was detected.

- start_s:

  Numeric. Spindle onset relative to epoch start (seconds).

- end_s:

  Numeric. Spindle offset relative to epoch start (seconds).

- duration_s:

  Numeric. Spindle duration in seconds.

- peak_freq_hz:

  Numeric. Peak frequency within the spindle.

- rms_uv:

  Numeric. Mean RMS amplitude within the spindle (µV).

- channel:

  Character. Channel label.

## References

Lacourse, K., Yetton, B., Mednick, S., & Bhatt, D. L. (2019). *Massive
online sleep staging: A polysomnography data repository*. Journal of
Sleep Research.

Vallat, R., & Walker, M. P. (2021). An open-source, high-performance
tool for automated sleep staging. *eLife*, 10, e70092.
[doi:10.7554/eLife.70092](https://doi.org/10.7554/eLife.70092)
