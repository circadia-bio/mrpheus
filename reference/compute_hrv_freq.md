# Frequency-domain HRV (HF power), locked to each infant's own respiratory rate, as a proxy for respiratory sinus arrhythmia (RSA).

Adults have a standard HF band (0.15-0.4 Hz) because adult respiratory
rate is fairly stable across that range. Neonates breathe far faster
(roughly 40-70 breaths/min at rest, i.e. ~0.7-1.2 Hz, and higher still
if tachypnoeic), well above the adult band. Using the adult HF band on
infant data would largely miss the respiratory-linked HRV signal
entirely. This instead detects each infant's own dominant respiratory
frequency directly from the physlog's `resp` channel, and integrates
RR-interval spectral power in a narrow band centered on that frequency.

## Usage

``` r
compute_hrv_freq(
  qrs,
  resp_signal,
  fs,
  resample_hz = 4,
  band_halfwidth_hz = 0.1,
  resp_search_range_hz = c(0.3, 2.5)
)
```

## Arguments

- qrs:

  An `mrpheus_qrs` object (output of
  [`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md)),
  giving R-peak sample indices and sampling rate.

- resp_signal:

  Numeric vector, the raw respiration channel, same sampling rate and
  duration as the signal QRS detection was run on.

- fs:

  Sampling rate (Hz) of `resp_signal` (should match the ECG's original
  sampling rate; physlog channels share one clock).

- resample_hz:

  Frequency (Hz) to resample the RR-interval series to before spectral
  estimation. Default 4 Hz, standard for HRV frequency analysis and
  comfortably above twice the plausible infant respiratory range
  (Nyquist-safe up to 2 Hz / 120 breaths/min).

- band_halfwidth_hz:

  Half-width (Hz) of the HF integration band, centered on the detected
  respiratory frequency. Default 0.1 Hz.

- resp_search_range_hz:

  Plausible neonatal respiratory frequency range (Hz) to search for the
  dominant peak in `resp_signal`. Default c(0.3, 2.5) (18-150
  breaths/min), deliberately wide to avoid missing tachypnoeic or
  bradypnoeic infants.

## Value

A list: `hf_power` (integrated PSD power in the respiration-locked
band), `resp_freq_hz` (detected dominant respiratory frequency),
`band_low`/`band_high` (integration band bounds), `n_beats_used`.

## Details

Verified against synthetic data: recovers an injected respiratory
frequency exactly (0 Hz error against a known 0.9 Hz test signal), and
HF power scales with injected RSA amplitude as expected (power scaling
roughly with amplitude^2).
