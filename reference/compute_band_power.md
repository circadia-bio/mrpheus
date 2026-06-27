# Compute EEG band power per epoch

Estimates power spectral density (PSD) using Welch's method and
integrates power within standard EEG frequency bands (delta, theta,
alpha, sigma, beta, gamma) for each epoch and each specified channel.
Mirrors the band-power feature extraction used by YASA's staging
pipeline.

## Usage

``` r
compute_band_power(
  psg,
  channels = NULL,
  bands = list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 13), sigma = c(13, 16),
    beta = c(16, 30), gamma = c(30, 40)),
  relative = FALSE,
  window_s = 4,
  overlap = 0.5
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- channels:

  Character vector. EEG channel labels. If `NULL` (default), all non-bad
  EEG channels are used.

- bands:

  Named list of length-2 numeric vectors defining frequency bands.
  Default:

      list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 13),
           sigma = c(13, 16), beta = c(16, 30), gamma = c(30, 40))

- relative:

  Logical. If `TRUE`, return relative band power (band / total power in
  0.5–40 Hz). Default `FALSE`.

- window_s:

  Numeric. Welch window length in seconds. Default `4`.

- overlap:

  Numeric in \[0, 1). Fractional overlap between Welch windows. Default
  `0.5`.

## Value

A tibble with columns `epoch`, `channel`, one column per band, and
`total_power`. Units are µV²/Hz (or dimensionless if `relative = TRUE`).

## Examples

``` r
if (FALSE) { # \dontrun{
bp <- compute_band_power(psg)
bp <- compute_band_power(psg, channels = "EEG Fpz-Cz", relative = TRUE)
} # }
```
