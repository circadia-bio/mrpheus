# Package index

## Ingestion

Read and prepare PSG recordings.

- [`read_edf()`](https://mrpheus.circadia-lab.uk/reference/read_edf.md)
  : Read an EDF or EDF+ recording
- [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md)
  : Prepare a PSG recording for analysis

## Artefacts

Detect and flag artefact epochs.

- [`detect_artifacts()`](https://mrpheus.circadia-lab.uk/reference/detect_artifacts.md)
  : Detect artefact epochs in a PSG recording

## Spectral analysis

Power spectral density and time-frequency representations.

- [`compute_band_power()`](https://mrpheus.circadia-lab.uk/reference/compute_band_power.md)
  : Compute EEG band power per epoch
- [`compute_spectrogram()`](https://mrpheus.circadia-lab.uk/reference/compute_spectrogram.md)
  : Compute a time-frequency spectrogram

## Event detection

Automatic detection of sleep microstructure events.

- [`compute_spindles()`](https://mrpheus.circadia-lab.uk/reference/compute_spindles.md)
  : Detect sleep spindles
- [`compute_slow_oscillations()`](https://mrpheus.circadia-lab.uk/reference/compute_slow_oscillations.md)
  : Detect slow oscillations

## Sleep staging

Automatic AASM epoch staging and hypnogram export.

- [`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
  : Automatic AASM sleep staging
- [`export_hypnogram()`](https://mrpheus.circadia-lab.uk/reference/export_hypnogram.md)
  : Export a staged hypnogram to hypnor

## Respiratory

Apnea/hypopnea detection and respiratory indices.

- [`detect_apneas()`](https://mrpheus.circadia-lab.uk/reference/detect_apneas.md)
  : Detect respiratory events (apneas and hypopneas)
- [`compute_ahi()`](https://mrpheus.circadia-lab.uk/reference/compute_ahi.md)
  : Compute Apnea-Hypopnea Index (AHI)
- [`compute_odi()`](https://mrpheus.circadia-lab.uk/reference/compute_odi.md)
  : Compute Oxygen Desaturation Index (ODI)

## Cardiac

HRV metrics across sleep stages.

- [`compute_hrv_sleep()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_sleep.md)
  : Compute HRV metrics across sleep stages

## Data

Bundled datasets and palettes.

- [`palette_orpheus`](https://mrpheus.circadia-lab.uk/reference/palette_orpheus.md)
  : Orpheus mosaic palette

## Package

Package-level documentation.

- [`mrpheus`](https://mrpheus.circadia-lab.uk/reference/mrpheus-package.md)
  [`mrpheus-package`](https://mrpheus.circadia-lab.uk/reference/mrpheus-package.md)
  : mrpheus: Polysomnography Signal Analysis for Sleep Research
