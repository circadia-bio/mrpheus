# Compute HRV metrics across sleep stages

Extracts R-peak positions from the ECG channel, computes standard time-
and frequency-domain heart rate variability (HRV) metrics, and
stratifies results by sleep stage.

## Usage

``` r
compute_hrv_sleep(
  psg,
  staging = NULL,
  ecg_channel = NULL,
  min_rr_ms = 300,
  max_rr_ms = 2000
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- staging:

  Tibble from
  [`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
  or `NULL`. If `NULL`, HRV is computed across all epochs without stage
  stratification.

- ecg_channel:

  Character or `NULL`. ECG channel label. If `NULL` (default), the first
  non-bad ECG channel is used.

- min_rr_ms:

  Numeric. Minimum physiologically plausible RR interval (ms). Default
  `300` (200 bpm ceiling).

- max_rr_ms:

  Numeric. Maximum physiologically plausible RR interval (ms). Default
  `2000` (30 bpm floor).

## Value

A tibble with one row per sleep stage (or one row if `staging` is
`NULL`):

- stage:

  Character. AASM stage or `"ALL"`.

- n_epochs:

  Integer. Number of epochs in this stage.

- mean_rr_ms:

  Numeric. Mean RR interval (ms).

- sdnn_ms:

  Numeric. SDNN — SD of all NN intervals.

- rmssd_ms:

  Numeric. RMSSD — root mean square of successive differences.

- lf_power:

  Numeric. LF band power (0.04–0.15 Hz).

- hf_power:

  Numeric. HF band power (0.15–0.4 Hz).

- lf_hf_ratio:

  Numeric. LF/HF ratio.
