# Compute EEG band power per epoch

Estimates power spectral density (PSD) using Welch's method and
integrates power within standard EEG frequency bands for each epoch and
channel.

## Usage

``` r
compute_band_power(
  psg,
  channels = NULL,
  bands = list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 12), sigma = c(12, 16),
    beta = c(16, 30), gamma = c(30, 45)),
  relative = FALSE,
  win_sec = 4,
  noverlap = 200L,
  nfft = 1024L
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- channels:

  Character vector. Channel labels to include. If `NULL` (default), all
  non-bad EEG channels are used.

- bands:

  Named list of length-2 numeric vectors defining frequency bands in Hz.
  Default matches YASA's `bandpower()`:

      list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 12),
           sigma = c(12, 16), beta = c(16, 30), gamma = c(30, 45))

- relative:

  Logical. If `TRUE`, each band power is divided by the sum of all band
  powers (dimensionless). Matches YASA `relative = True`. Default
  `FALSE`.

- win_sec:

  Numeric. Welch window length in seconds. Default `4`.

- noverlap:

  Integer. Number of overlapping samples between Welch windows. Default
  `200`. Must be less than `win_sec * sample_rate`.

- nfft:

  Integer. FFT length. Default `1024`. If smaller than the window length
  in samples it is automatically raised to the next power of 2.

## Value

A tibble with columns `epoch`, `channel`, one column per named band, and
`total_power` (sum of all band powers, in V^2/Hz units before any
relative scaling).

## Details

The Welch implementation matches `scipy.signal.welch` with YASA's
default parameters (`win_sec = 4`, `noverlap = 200`, `nfft = 1024`, Hann
window, constant detrend, mean averaging). Relative power matches YASA's
`relative = True` behaviour: each band is divided by the sum of all band
powers (not the total PSD integral).

## Examples

``` r
if (FALSE) { # \dontrun{
bp <- compute_band_power(psg)
bp <- compute_band_power(psg, channels = "EEG Fpz-Cz", relative = TRUE)

# Custom bands
bp <- compute_band_power(
  psg,
  bands = list(slow = c(0.5, 1), delta = c(1, 4), theta = c(4, 8))
)
} # }
```
