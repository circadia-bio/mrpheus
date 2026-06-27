# Detect artefact epochs in a PSG recording

Flags epochs containing likely artefacts based on amplitude thresholds,
excessive high-frequency (muscle) power, and movement contamination.
Operates on the EEG channels by default. Returns a logical vector (one
value per epoch) and optionally attaches per-epoch artefact metrics for
inspection.

## Usage

``` r
detect_artifacts(
  psg,
  channels = NULL,
  amp_threshold_uv = 500,
  hf_band = c(40, 100),
  hf_percentile = 0.99,
  verbose = TRUE
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- channels:

  Character vector. Channel labels to evaluate. If `NULL` (default), all
  non-bad EEG channels are used.

- amp_threshold_uv:

  Numeric. Peak-to-peak amplitude threshold in µV. Epochs exceeding this
  in any selected channel are flagged. Default `500`.

- hf_band:

  Numeric vector of length 2. Frequency band (Hz) considered
  high-frequency / muscle contamination. Default `c(40, 100)`.

- hf_percentile:

  Numeric. Epochs whose HF power exceeds this percentile (across all
  epochs) are flagged. Default `0.99`.

- verbose:

  Logical. Print summary. Default `TRUE`.

## Value

A tibble with one row per epoch and columns:

- epoch:

  Integer. Epoch index (1-based).

- artefact:

  Logical. `TRUE` if the epoch is flagged as artefact.

- reason:

  Character. Comma-separated reasons for flagging, or `NA`.

- peak_to_peak_uv:

  Numeric. Maximum peak-to-peak amplitude across selected channels.

- hf_power_db:

  Numeric. Mean HF band power (dB) across selected channels.
