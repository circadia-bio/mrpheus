# Compute a time-frequency spectrogram

Generates a short-time Fourier transform (STFT) based spectrogram for a
single PSG channel across all epochs. Returns power (µV²/Hz) as a matrix
with time (epochs) on rows and frequency on columns.

## Usage

``` r
compute_spectrogram(
  psg,
  channel = NULL,
  freq_range = c(0, 40),
  window_s = 2,
  overlap = 0.5,
  db = TRUE
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- channel:

  Character. A single channel label.

- freq_range:

  Numeric vector of length 2. Frequency range to return (Hz). Default
  `c(0, 40)`.

- window_s:

  Numeric. STFT window length in seconds. Default `2`.

- overlap:

  Numeric in \[0, 1). Fractional window overlap. Default `0.5`.

- db:

  Logical. Return power in decibels (`10 * log10(power)`)? Default
  `TRUE`.

## Value

A list of class `mrpheus_spectrogram` with:

- power:

  Numeric matrix. Rows = epochs, columns = frequency bins.

- freqs:

  Numeric vector. Frequency bin centres (Hz).

- epochs:

  Integer vector. Epoch indices.

- channel:

  Character. The channel label.

- db:

  Logical. Whether power is in dB.
