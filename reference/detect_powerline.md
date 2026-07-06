# Auto-detect powerline frequency

Estimates whether dominant powerline interference is at 50 Hz
(Europe/UK) or 60 Hz (North America) by comparing PSD power at both
frequencies using Welch's method on the first 10 seconds of the signal.

## Usage

``` r
detect_powerline(signal, sr)
```

## Arguments

- signal:

  Numeric vector. Should be an EEG or broadband channel with sufficient
  bandwidth.

- sr:

  Numeric. Sampling rate in Hz. Must be \> 100.

## Value

Integer. `50L` or `60L`.

## Details

If the signal's Nyquist limit is below 60 Hz, returns `50L` immediately.

## Examples

``` r
if (FALSE) { # \dontrun{
detect_powerline(eeg_signal, sr = 256)
} # }
```
