# Changelog

## mrpheus 0.1.2

#### AASM sleep staging — validated parity with YASA

The
[`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
pipeline now achieves validated numerical parity with the YASA Python
reference implementation (Vallat & Walker, 2021), producing 77.8%
epoch-level agreement against scored ground truth on the Sleep-EDF
Cassette dataset — matching YASA’s own 77.7% on the same recording.

**Feature extraction fixes (149-feature LightGBM matrix):**

- Filter: replaced Butterworth with exact 825-tap MNE Hamming FIR
  coefficients bundled in `inst/filters/mne_bandpass_100hz.csv`; applied
  as a single zero-phase FFT pass (not `filtfilt`, which squares the
  response).
- EMG resampling: ported MNE’s `_resample_fft` algorithm exactly —
  npad=‘auto’ pads to the next power of 2 for fast FFTs, reflect-limited
  padding, correct conjugate-Nyquist placement in the zero-padded
  frequency domain.
- Welch PSD: periodic (DFT-even) Hamming window matching
  `scipy.signal.get_window('hamming', N, fftbins=True)`; segment
  detrending (`detrend='constant'`); bias-corrected median averaging.
- Band integration: composite Simpson’s rule matching
  `scipy.integrate.simpson` (odd N: standard 1/3; even N: 1/3 on first
  N-3 points + 3/8 on last 4); inclusive frequency bounds.
- Petrosian FD: fixed to count local extrema (sign changes in
  consecutive first differences) rather than zero crossings of x —
  matching antropy’s convention.
- Rolling normalisation: `_c7min_norm` now uses a C++ triangular rolling
  mean via `roll_triang_mean_cpp`; `_p2min_norm` uses right-aligned
  [`zoo::rollapply`](https://rdrr.io/pkg/zoo/man/rollapply.html); both
  followed by `robust_scale`.

**Model prediction fixes:**

- Feature ordering:
  [`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
  reads the model’s `feature_names=` line from the bundled `.txt` file
  and reorders columns to match the EEG \| time \| EOG \| EMG order the
  model was trained with (LightGBM R passes matrices by position, not by
  name).
- LightGBM v4 compatibility:
  [`predict()`](https://rdrr.io/r/stats/predict.html) now returns a
  matrix directly; removed the deprecated `reshape=TRUE` argument.
- Label mapping: class indices follow sklearn `LabelEncoder`
  alphabetical order — N1=0, N2=1, N3=2, REM=3, W=4.

**Performance:**

- Three inner loops ported to Rcpp C++: `perm_entropy_cpp`,
  `higuchi_fd_cpp`, `roll_triang_mean_cpp`. Full-day staging (2650
  epochs × 3 channels) completes in ~40 s.

**Documentation:**

- New vignette: *Automatic AASM Sleep Staging* — full pipeline API with
  PhysioNet download instructions.
- New pkgdown article: *Sleep Staging: Live Walkthrough (SC4001E0)* —
  fully executed on a real Sleep-EDF cassette recording with posterior
  probability figures.
- 56 new unit tests covering all staging feature helpers, Rcpp
  functions, and edge cases. `test-staging-features.R` now has 149
  passing tests.

## mrpheus 0.1.1

#### Scope

- Package reframed from PSG-only to **raw physiological signal analysis
  for biological rhythms research**, encompassing PSG, MRI-concurrent
  physiology, and EEG.
- DESCRIPTION, README, pkgdown site, and vignettes updated accordingly.

#### New functions

- [`read_philips_physlog()`](https://mrpheus.circadia-lab.uk/reference/read_philips_physlog.md)
  — reads Philips PMU `.log` files (wBTU wireless VCG, wired ECG, or
  custom sampling rate). Returns a `mrpheus_physlog` S3 object with
  signal matrix, event markers, and header metadata.
- [`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md)
  — Pan-Tompkins QRS detector (Pan & Tompkins, 1985). Bandpass filter →
  derivative → squaring → 150 ms moving average → adaptive
  dual-threshold with T-wave rejection and search-back. Returns a
  `mrpheus_qrs` S3 object.
- [`compute_hr_signal()`](https://mrpheus.circadia-lab.uk/reference/compute_hr_signal.md)
  — converts R-peak indices to a sample-by-sample instantaneous heart
  rate trace (bpm).

#### Other changes

- Staging model (`inst/models/yasa_staging.txt`) is now bundled with the
  package — no Python setup required after installation.
- `data-raw/fetch_yasa_model.py` updated for YASA 0.7.0.
- New vignette: *MRI physiological signal processing*.
- New test suite: `test-read_philips_physlog.R`, `test-detect_qrs.R`,
  `test-compute_hr_signal.R` with MATLAB bit-perfect fixture comparison.
- `boldR` added to ecosystem table in docs and vignettes.

## mrpheus 0.1.0

- Initial release. Package scaffolded with full PSG analysis pipeline.
- [`read_edf()`](https://mrpheus.circadia-lab.uk/reference/read_edf.md)
  /
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md)
  — EDF/EDF+ ingestion, channel inventory, epoch segmentation, and
  bad-channel flagging.
- [`detect_artifacts()`](https://mrpheus.circadia-lab.uk/reference/detect_artifacts.md)
  — amplitude and high-frequency artefact detection.
- [`compute_band_power()`](https://mrpheus.circadia-lab.uk/reference/compute_band_power.md)
  — Welch PSD with δ/θ/α/σ/β/γ bands per epoch.
- [`compute_spectrogram()`](https://mrpheus.circadia-lab.uk/reference/compute_spectrogram.md)
  — STFT-based time-frequency spectrogram.
- [`compute_spindles()`](https://mrpheus.circadia-lab.uk/reference/compute_spindles.md)
  — RMS envelope spindle detection (Lacourse et al., 2019).
- [`compute_slow_oscillations()`](https://mrpheus.circadia-lab.uk/reference/compute_slow_oscillations.md)
  — zero-crossing SO detection (Mölle et al., 2002).
- [`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
  — automatic AASM staging via pre-trained LightGBM model ported from
  YASA (Vallat & Walker, 2021). Feature extraction parity with the
  Python pipeline is in progress.
- [`export_hypnogram()`](https://mrpheus.circadia-lab.uk/reference/export_hypnogram.md)
  — returns a `mrpheus_hypnogram` object ready for
  `hypnor::new_hypnogram()` once `hypnor` is available.
- [`detect_apneas()`](https://mrpheus.circadia-lab.uk/reference/detect_apneas.md),
  [`compute_ahi()`](https://mrpheus.circadia-lab.uk/reference/compute_ahi.md),
  [`compute_odi()`](https://mrpheus.circadia-lab.uk/reference/compute_odi.md)
  — respiratory stubs, full implementation pending.
- [`compute_hrv_sleep()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_sleep.md)
  — HRV stub, full implementation pending.
- `palette_orpheus` — 8-colour palette extracted from the Roman mosaic
  *Orpheus Charming the Animals* (3rd century AD, Palermo Archaeological
  Museum).
- Hex sticker: olive background, ivory lyre, ochre motion lines,
  vermillion border.
- pkgdown site deploying via Netlify at `mrpheus.circadia-lab.uk`.
