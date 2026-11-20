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
    beta = c(16, 30), gamma = c(30, 40)),
  relative = TRUE,
  win_sec = 4,
  noverlap = NULL,
  nfft = NULL,
  window = c("hamming", "hann"),
  average = c("median", "mean")
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
           sigma = c(12, 16), beta = c(16, 30), gamma = c(30, 40))

- relative:

  Logical. If `TRUE`, each band power is divided by the PSD integral
  between the lowest and highest band edge (dimensionless). Matches YASA
  `relative = True`. Default `TRUE`, as in YASA; `FALSE` gives absolute
  power.

- win_sec:

  Numeric. Welch window length in seconds. Default `4`.

- noverlap:

  Integer or `NULL`. Number of overlapping samples between Welch
  windows. `NULL` (default) uses half the window, as
  `scipy.signal.welch` does. Must be less than `win_sec * sample_rate`.

- nfft:

  Integer or `NULL`. FFT length. `NULL` (default) uses the window length
  in samples. If smaller than the window length it is automatically
  raised to the next power of 2.

- window:

  Character. Welch window, `"hamming"` (default, YASA's) or `"hann"`.
  Both are the periodic form that `scipy.signal.welch` uses.

- average:

  Character. How the segment periodograms are combined: `"median"`
  (default, YASA's; bias-corrected as in scipy) or `"mean"`.

## Value

A tibble with columns `epoch`, `channel`, one column per named band, and
`total_power` (the PSD integral between the lowest and highest band
edge, in V^2 units before any relative scaling; YASA's `TotalAbsPow`).

## Details

The defaults follow `yasa.bandpower()` (YASA 0.8.1): the PSD is
`scipy.signal.welch()` with a 4 s window (`win_sec = 4`), a periodic
Hamming window, median averaging over the overlapping segments (scipy's
bias-corrected median), an overlap of half a window, an FFT length equal
to the window length, and constant detrending. Band power is the Simpson
integral of the PSD, and the default bands are YASA's, including gamma
at 30-40 Hz. Relative power matches YASA's `relative = True`: each band
is divided by the integral of the PSD between the lowest and the highest
band edge (0.5-40 Hz by default). With bands that tile that range and
edges on FFT bins this equals the sum of the band powers; otherwise it
is larger.

To use another estimate, change `window` (`"hann"`), `average`
(`"mean"`), `noverlap` or `nfft`.

## Examples

``` r
if (FALSE) { # \dontrun{
bp <- compute_band_power(psg)
bp <- compute_band_power(psg, channels = "EEG Fpz-Cz", relative = FALSE)

# Custom bands
bp <- compute_band_power(
  psg,
  bands = list(slow = c(0.5, 1), delta = c(1, 4), theta = c(4, 8))
)

# The mean of Hann-windowed segments instead of YASA's defaults
bp <- compute_band_power(psg, window = "hann", average = "mean")
} # }
```
