# Apply a notch (bandstop) filter to a signal

Designs a 2nd-order Butterworth bandstop filter centred at `freq` Hz
with bandwidth `bandwidth_hz` and applies it zero-phase using
[`gsignal::filtfilt()`](https://rdrr.io/pkg/gsignal/man/filtfilt.html).
Because the filter has only ~9 coefficients, it is safe to use on long
recordings without memory issues.

## Usage

``` r
notch_filter(signal, sr, freq, bandwidth_hz = 2, harmonics = TRUE)
```

## Arguments

- signal:

  Numeric vector.

- sr:

  Numeric. Sampling rate in Hz.

- freq:

  Numeric. Centre frequency in Hz (e.g. `50` or `60`).

- bandwidth_hz:

  Numeric. Full bandwidth of the notch in Hz. Default `2`.

- harmonics:

  Logical. Also notch harmonics (2×, 3×, ...) below Nyquist − 5 Hz.
  Default `TRUE`.

## Value

Filtered numeric vector, same length as `signal`.

## Details

By default also applies the filter at all harmonics of `freq` below the
Nyquist limit, which is the standard approach for powerline noise
removal.

## Examples

``` r
if (FALSE) { # \dontrun{
clean <- notch_filter(eeg_signal, sr = 256, freq = 50)
clean <- notch_filter(eeg_signal, sr = 256, freq = 60, harmonics = FALSE)
} # }
```
