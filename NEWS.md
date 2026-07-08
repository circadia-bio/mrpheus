# mrpheus 0.1.5

### Rcpp hot-path optimisation (Tiers 2-4)

All remaining per-recording hot paths ported to C++:

- `rowmedian_cpp`: row-wise median of the Welch periodogram matrix, replacing
  `apply(pgrams, 1L, stats::median)` — eliminates ~750 000 R `median()` calls
  per recording (251 freq bins x ~3 000 epochs)
- `roll_right_mean_cpp`: right-aligned partial rolling mean (k=4), replacing
  `zoo::rollapply` for the `_p2min_norm` normalisation step (49 calls/recording)
- `robust_scale_cpp`: `(x - median) / (q95 - q5)` with type-7 quantile and
  NA-awareness, replacing `median()` + `quantile()` for 49 normalisation calls
- `roll_rms_cpp`: centered rolling RMS replacing `zoo::rollapply` in
  `compute_spindles()` — runs on the full concatenated epoch signal per channel
- `detect_so_candidates_cpp`: zero-crossing scan + duration/amplitude filtering
  for slow oscillation detection, replacing the R `which()` + `lapply` loop
  in `compute_slow_oscillations()`

New file `src/event_detection.cpp` houses the event detection hot paths.
Removed `zoo` from Imports — no longer used anywhere in the package.

All new functions validated with `devtools::check()`: 0 errors | 0 warnings | 0 notes.

# mrpheus 0.1.4

### Rcpp hot-path optimisation (Tier 1)

Four more per-epoch feature functions ported to C++ for the staging pipeline:

- `nzc_cpp`: zero-crossing count — single O(N) pass replacing `sum(diff(sign(x)) != 0)`
- `petrosian_fd_cpp`: local extrema count + log formula — replaces `diff` + logical sum in R
- `hjorth_cpp`: mobility and complexity — single O(N) pass accumulating all three variances
  (x, diff(x), diff(diff(x))) without materialising intermediate vectors
- `stat_features_cpp`: std, IQR (type-7 quantile), skewness, excess kurtosis — replaces
  `sd()`, `IQR()` (which sorts), and two `mean()` calls per epoch

All four are called ~3 000 times per recording (once per epoch × channel). Fixed a latent
bug where the pure-R `vapply` definition of `.roll_triang_mean` in the normalisation
section was shadowing `roll_triang_mean_cpp` on every `load_all()`, meaning the Rcpp
version was never actually called. Now a single Rcpp-backed definition exists.

All 22 parity checks pass (`dev/test_tier1_cpp.R`). `0 errors | 0 warnings | 0 notes`.

# mrpheus 0.1.3

### PSG preprocessing pipeline

Full preprocessing pipeline from EDF to analysis-ready epochs, matching the
companion Python/MNE notebook workflow.

* `preprocess_psg()` — continuous-signal preprocessing pipeline: channel
  renaming, DC removal, powerline notch filter (auto-detected 50/60 Hz +
  harmonics), and channel-type-specific bandpass filtering (EEG/EOG/EMG/ECG).
  Operates on the full continuous signal before re-epoching to avoid filter
  discontinuities at epoch boundaries.
* `remove_dc()`, `detect_powerline()`, `notch_filter()`, `bandpass_filter()` —
  exported signal-level functions; each operates on a plain numeric vector
  independently of the PSG pipeline.

### EOG artefact correction

* `correct_eog_regression()` — removes ocular artefacts via multiple linear
  regression of each EEG channel on the EOG reference channels. Equivalent to
  MNE's SSP projection approach without requiring blink event detection.
* `correct_eog_ica()` — ICA-based EOG removal via `fastICA`. Components whose
  time courses correlate with EOG channels above a threshold (default 0.35,
  matching the notebook) are identified and subtracted.

### Spectral analysis

* `compute_band_power()` — rewritten Welch PSD matching `scipy.signal.welch`
  with YASA defaults (`win_sec = 4`, `noverlap = 200`, `nfft = 1024`, Hann
  window, constant detrend, mean averaging). Band integration now uses
  composite Simpson's rule (`.scipy_simpson`) matching `scipy.integrate.simpson`.
  Parity validated against `yasa.bandpower()`: absolute power MARE < 0.08 %,
  relative power MARE < 0.8 % across all six standard bands.
* `compute_temporal_bandpower()` — sliding-window band power analysis across
  the full recording night. Returns a long-format tibble with one row per
  (window, band) including `time_hours`, `dominant_stage`, `power`, and
  `relative_power`. Mirrors `calculate_temporal_bandpower()` from the
  companion notebook.

### Dependencies

* Added `fastICA (>= 1.2-0)` to `Imports` for `correct_eog_ica()`.

---

# mrpheus 0.1.2

### AASM sleep staging — validated parity with YASA

The `stage_epochs()` pipeline now achieves validated numerical parity with
the YASA Python reference implementation (Vallat & Walker, 2021), producing
77.8% epoch-level agreement against scored ground truth on the Sleep-EDF
Cassette dataset — matching YASA's own 77.7% on the same recording.

**Feature extraction fixes (149-feature LightGBM matrix):**

* Filter: replaced Butterworth with exact 825-tap MNE Hamming FIR coefficients
  bundled in `inst/filters/mne_bandpass_100hz.csv`; applied as a single
  zero-phase FFT pass (not `filtfilt`, which squares the response).
* EMG resampling: ported MNE's `_resample_fft` algorithm exactly — npad='auto'
  pads to the next power of 2 for fast FFTs, reflect-limited padding, correct
  conjugate-Nyquist placement in the zero-padded frequency domain.
* Welch PSD: periodic (DFT-even) Hamming window matching
  `scipy.signal.get_window('hamming', N, fftbins=True)`; segment detrending
  (`detrend='constant'`); bias-corrected median averaging.
* Band integration: composite Simpson's rule matching `scipy.integrate.simpson`
  (odd N: standard 1/3; even N: 1/3 on first N-3 points + 3/8 on last 4);
  inclusive frequency bounds.
* Petrosian FD: fixed to count local extrema (sign changes in consecutive
  first differences) rather than zero crossings of x — matching antropy's
  convention.
* Rolling normalisation: `_c7min_norm` now uses a C++ triangular rolling
  mean via `roll_triang_mean_cpp`; `_p2min_norm` uses right-aligned
  `zoo::rollapply`; both followed by `robust_scale`.

**Model prediction fixes:**

* Feature ordering: `stage_epochs()` reads the model's `feature_names=` line
  from the bundled `.txt` file and reorders columns to match the EEG | time |
  EOG | EMG order the model was trained with (LightGBM R passes matrices by
  position, not by name).
* LightGBM v4 compatibility: `predict()` now returns a matrix directly;
  removed the deprecated `reshape=TRUE` argument.
* Label mapping: class indices follow sklearn `LabelEncoder` alphabetical
  order — N1=0, N2=1, N3=2, REM=3, W=4.

**Performance:**

* Three inner loops ported to Rcpp C++: `perm_entropy_cpp`,
  `higuchi_fd_cpp`, `roll_triang_mean_cpp`. Full-day staging
  (2650 epochs × 3 channels) completes in ~40 s.

**Documentation:**

* New vignette: *Automatic AASM Sleep Staging* — full pipeline API with
  PhysioNet download instructions.
* New pkgdown article: *Sleep Staging: Live Walkthrough (SC4001E0)* —
  fully executed on a real Sleep-EDF cassette recording with posterior
  probability figures.
* 56 new unit tests covering all staging feature helpers, Rcpp functions,
  and edge cases. `test-staging-features.R` now has 149 passing tests.

# mrpheus 0.1.1

### Scope

* Package reframed from PSG-only to **raw physiological signal analysis for
  biological rhythms research**, encompassing PSG, MRI-concurrent physiology,
  and EEG.
* DESCRIPTION, README, pkgdown site, and vignettes updated accordingly.

### New functions

* `read_philips_physlog()` — reads Philips PMU `.log` files (wBTU wireless
  VCG, wired ECG, or custom sampling rate). Returns a `mrpheus_physlog`
  S3 object with signal matrix, event markers, and header metadata.
* `detect_qrs()` — Pan-Tompkins QRS detector (Pan & Tompkins, 1985).
  Bandpass filter → derivative → squaring → 150 ms moving average →
  adaptive dual-threshold with T-wave rejection and search-back.
  Returns a `mrpheus_qrs` S3 object.
* `compute_hr_signal()` — converts R-peak indices to a sample-by-sample
  instantaneous heart rate trace (bpm).

### Other changes

* Staging model (`inst/models/yasa_staging.txt`) is now bundled with the
  package — no Python setup required after installation.
* `data-raw/fetch_yasa_model.py` updated for YASA 0.7.0.
* New vignette: *MRI physiological signal processing*.
* New test suite: `test-read_philips_physlog.R`, `test-detect_qrs.R`,
  `test-compute_hr_signal.R` with MATLAB bit-perfect fixture comparison.
* `boldR` added to ecosystem table in docs and vignettes.

# mrpheus 0.1.0

* Initial release. Package scaffolded with full PSG analysis pipeline.
* `read_edf()` / `prepare_psg()` — EDF/EDF+ ingestion, channel inventory,
  epoch segmentation, and bad-channel flagging.
* `detect_artifacts()` — amplitude and high-frequency artefact detection.
* `compute_band_power()` — Welch PSD with δ/θ/α/σ/β/γ bands per epoch.
* `compute_spectrogram()` — STFT-based time-frequency spectrogram.
* `compute_spindles()` — RMS envelope spindle detection (Lacourse et al., 2019).
* `compute_slow_oscillations()` — zero-crossing SO detection (Mölle et al., 2002).
* `stage_epochs()` — automatic AASM staging via pre-trained LightGBM model
  ported from YASA (Vallat & Walker, 2021). Feature extraction parity with
  the Python pipeline is in progress.
* `export_hypnogram()` — returns a `mrpheus_hypnogram` object ready for
  `hypnor::new_hypnogram()` once `hypnor` is available.
* `detect_apneas()`, `compute_ahi()`, `compute_odi()` — respiratory stubs,
  full implementation pending.
* `compute_hrv_sleep()` — HRV stub, full implementation pending.
* `palette_orpheus` — 8-colour palette extracted from the Roman mosaic
  *Orpheus Charming the Animals* (3rd century AD, Palermo Archaeological Museum).
* Hex sticker: olive background, ivory lyre, ochre motion lines, vermillion border.
* pkgdown site deploying via Netlify at `mrpheus.circadia-lab.uk`.
