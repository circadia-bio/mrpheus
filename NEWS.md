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
