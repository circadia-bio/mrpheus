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
