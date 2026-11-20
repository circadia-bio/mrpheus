# Time-resolved (sliding-window) respiration-locked HF power.

Same underlying method as
[`compute_hrv_freq()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_freq.md)
(RR intervals resampled to an even grid, spectral power integrated in a
band centered on the infant's own respiratory frequency), but instead of
averaging across the whole recording into one number, computes HF power
in sliding windows, giving a time series that can be compared against
another time-varying signal (e.g. the Kuramoto synchrony trace from
`dynR::leida_pipeline()`) within the same recording.

## Usage

``` r
compute_hrv_freq_windowed(
  qrs,
  resp_signal,
  fs,
  window_sec = 60,
  step_sec = 15,
  resample_hz = 4,
  band_halfwidth_hz = 0.1,
  resp_search_range_hz = c(0.3, 2.5)
)
```

## Arguments

- qrs:

  An `mrpheus_qrs` object (output of
  [`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md)).

- resp_signal:

  Numeric vector, the raw respiration channel.

- fs:

  Sampling rate (Hz) of `resp_signal`.

- window_sec:

  Window length in seconds. Default 60s: long enough for ~0.017 Hz
  frequency resolution (1/window_sec), comfortably finer than the 0.2
  Hz-wide integration band, short enough to give reasonable temporal
  resolution over a ~15 minute scan.

- step_sec:

  Step between window centers in seconds. Default 15s (75% overlap
  between consecutive windows).

- resample_hz:

  RR-interval resampling rate. Default 4 Hz.

- band_halfwidth_hz:

  Half-width of the integration band. Default 0.1 Hz.

- resp_search_range_hz:

  Plausible respiratory frequency range for the whole-recording
  frequency detection. Default c(0.3, 2.5).

## Value

A list: `time` (window center times, seconds, from recording start),
`hf_power` (HF power at each window), `resp_freq_hz`, `band_low`,
`band_high`.

## Details

The respiratory frequency (and therefore the integration band) is
detected once from the *whole* recording's `resp` channel and held fixed
across all windows, rather than re-detected per window: a single
spectral estimate from the whole recording is far more stable than one
from a 60-second segment, and respiratory rate is generally stable
within a single scan session for a sleeping/sedated infant. If that
assumption doesn't hold for a given recording (e.g. a clear rate change
partway through), this will misrepresent windows far from the global
rate, worth checking the QC plot of `resp_freq_hz` if a session's
within-scan behaviour is in question.

Verified against synthetic data: a ramping RSA amplitude trajectory
(rising vagal tone over a 15-minute recording) is recovered with r=0.998
between windowed HF power and the true underlying amplitude-squared
trajectory.
