# Changelog

## mrpheus (development version)

### Initial scaffolding

- Package scaffolded with full PSG analysis pipeline.
- [`read_edf()`](https://mrpheus.circadia-lab.uk/reference/read_edf.md)
  /
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md)
  — EDF/EDF+ ingestion, channel inventory, epoch segmentation,
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
  — automatic AASM staging via pre-trained LightGBM model (ported from
  YASA; Vallat & Walker, 2021). Feature extraction parity with the
  Python pipeline is in progress.
- [`export_hypnogram()`](https://mrpheus.circadia-lab.uk/reference/export_hypnogram.md)
  — handoff to `hypnor`.
- [`detect_apneas()`](https://mrpheus.circadia-lab.uk/reference/detect_apneas.md),
  [`compute_ahi()`](https://mrpheus.circadia-lab.uk/reference/compute_ahi.md),
  [`compute_odi()`](https://mrpheus.circadia-lab.uk/reference/compute_odi.md)
  — respiratory stubs.
- [`compute_hrv_sleep()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_sleep.md)
  — HRV stub.
- `palette_orpheus` — 8-colour palette extracted from the Roman mosaic
  *Orpheus Charming the Animals* (3rd century AD, Palermo Archaeological
  Museum).
- Hex sticker: olive background, ivory lyre, ochre motion lines,
  vermillion border.
- pkgdown site deploying via Netlify at `mrpheus.circadia-lab.uk`.
