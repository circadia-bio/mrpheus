# Apply a bandpass filter to a signal

Designs a 4th-order Butterworth bandpass filter and applies it
zero-phase using
[`gsignal::filtfilt()`](https://rdrr.io/pkg/gsignal/man/filtfilt.html).
Handles edge cases gracefully:

## Usage

``` r
bandpass_filter(signal, sr, low_hz, high_hz)
```

## Arguments

- signal:

  Numeric vector.

- sr:

  Numeric. Sampling rate in Hz.

- low_hz:

  Numeric. Low cutoff frequency in Hz. Use `0` or negative to apply a
  low-pass filter only.

- high_hz:

  Numeric. High cutoff frequency in Hz.

## Value

Filtered numeric vector, same length as `signal`.

## Details

- `low_hz <= 0` → low-pass at `high_hz`

- `high_hz >= Nyquist` → high-pass at `low_hz`

- Otherwise → bandpass between `low_hz` and `high_hz`

The `high_hz` cutoff is clamped to 99 % of Nyquist to avoid instability.
Because this uses a short IIR filter, it is safe on long recordings.

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard EEG bandpass
filtered <- bandpass_filter(eeg_signal, sr = 256, low_hz = 0.3, high_hz = 35)

# High-pass only (baseline drift removal)
hp <- bandpass_filter(eeg_signal, sr = 256, low_hz = 0.1, high_hz = Inf)
} # }
```
