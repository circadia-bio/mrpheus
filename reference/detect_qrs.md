# Detect QRS Complexes in an ECG Signal (Pan-Tompkins Algorithm)

Detects R-peaks in a raw ECG signal using the Pan-Tompkins algorithm.
Implements bandpass filtering, derivative, squaring, moving-average
integration, and adaptive dual-threshold peak detection with T-wave
rejection and search-back for missed beats.

## Usage

``` r
detect_qrs(ecg, fs)
```

## Arguments

- ecg:

  Numeric vector. Raw ECG signal.

- fs:

  Numeric. Sampling frequency in Hz. If `ecg` comes from
  [`read_philips_physlog()`](https://mrpheus.circadia-lab.uk/reference/read_philips_physlog.md),
  use `rec$HDR$sfreq`.

## Value

A list of class `mrpheus_qrs` with components:

- `qrs_i`:

  Integer vector of R-peak sample indices (1-based), referenced to the
  bandpass-filtered signal.

- `qrs_amp`:

  Numeric vector of R-wave amplitudes at each detected peak (normalised
  bandpass signal units).

- `delay`:

  Numeric. Sample delay introduced by the moving-average stage (half the
  150 ms window length).

- `fs`:

  Numeric. Sampling frequency passed to the function.

## Details

Processing pipeline (Pan & Tompkins, 1985):

1.  **Bandpass filter** 5–15 Hz (Butterworth N=3). At fs = 200 Hz,
    implemented as separate 12 Hz low-pass then 5 Hz high-pass stages.

2.  **Derivative filter** — 5-point kernel interpolated to match fs.

3.  **Squaring** — non-linear emphasis of dominant peaks.

4.  **150 ms moving-average integration.**

5.  **Adaptive thresholding** — dual-threshold scheme with T-wave
    rejection and search-back for missed beats (initialised on 2-second
    training window).

R-peak indices are returned relative to the bandpass-filtered signal.
Account for `delay` if you need indices aligned to the raw ECG.

## References

Pan J, Tompkins WJ. A real-time QRS detection algorithm. *IEEE Trans
Biomed Eng.* 1985;32(3):230–236.
[doi:10.1109/TBME.1985.325532](https://doi.org/10.1109/TBME.1985.325532)

## See also

[`read_philips_physlog()`](https://mrpheus.circadia-lab.uk/reference/read_philips_physlog.md),
[`compute_hr_signal()`](https://mrpheus.circadia-lab.uk/reference/compute_hr_signal.md),
[`compute_hrv_sleep()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_sleep.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rec  <- read_philips_physlog("sub-01_physlog.log")
qrs  <- detect_qrs(rec$C[, "v1raw"], fs = rec$HDR$sfreq)
qrs

# Downstream HR signal
hr <- compute_hr_signal(qrs)
} # }
```
