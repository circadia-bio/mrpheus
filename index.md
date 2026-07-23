# 🪉 mrpheus

**Raw physiological signal analysis for biological rhythms research.**

[![r-universe](https://circadia-bio.r-universe.dev/badges/mrpheus)](https://circadia-bio.r-universe.dev/mrpheus)
[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3)](https://www.r-project.org/)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://mrpheus.circadia-lab.uk/LICENSE)
[![Circadia
Lab](https://img.shields.io/badge/circadia--bio-GitHub-014370)](https://github.com/circadia-bio)

------------------------------------------------------------------------

## 📖 What is mrpheus?

`mrpheus` is the raw signal layer of the [Circadia
Lab](https://github.com/circadia-bio) R ecosystem. It ingests and
processes multi-modal physiological recordings — polysomnography
(EDF/EDF+), MRI-concurrent physiological logs (Philips PMU), and EEG —
and extracts features across three biological rhythm domains:
**cardiac** (QRS detection, HRV), **respiratory** (apnoea detection,
respiratory metrics), and **neural** (sleep spindles, slow oscillations,
automatic AASM sleep staging).

Staged hypnograms and derived metrics pass downstream to `hypnor` and
`syncR` for cross-modal linking.

The name is spelled as *Morpheus* but carries a silent **m**, pronounced
as *Orpheus* — carrying both myths at once. Morpheus, god of dreams,
gives the package its subject; Orpheus, who descended into the
underworld to navigate the unconscious, gives it its spirit.

------------------------------------------------------------------------

## ✨ Features

- 📂 **EDF/EDF+ ingestion** —
  [`read_edf()`](https://mrpheus.circadia-lab.uk/reference/read_edf.md),
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md):
  channel inventory, epoch segmentation, bad-channel flagging
- 🧹 **Preprocessing pipeline** —
  [`preprocess_psg()`](https://mrpheus.circadia-lab.uk/reference/preprocess_psg.md):
  channel renaming, DC removal, powerline notch (50/60 Hz + harmonics),
  per-channel-type bandpass filtering; operates on the full continuous
  signal before re-epoching
- 🔊 **Signal filters** —
  [`remove_dc()`](https://mrpheus.circadia-lab.uk/reference/remove_dc.md),
  [`detect_powerline()`](https://mrpheus.circadia-lab.uk/reference/detect_powerline.md),
  [`notch_filter()`](https://mrpheus.circadia-lab.uk/reference/notch_filter.md),
  [`bandpass_filter()`](https://mrpheus.circadia-lab.uk/reference/bandpass_filter.md):
  exported functions that work on any plain numeric vector
- 👁️ **EOG artefact correction** —
  [`correct_eog_regression()`](https://mrpheus.circadia-lab.uk/reference/correct_eog_regression.md)
  (linear regression, no blink detection required),
  [`correct_eog_ica()`](https://mrpheus.circadia-lab.uk/reference/correct_eog_ica.md)
  (FastICA, threshold-based component rejection)
- 🚫 **Artefact detection** —
  [`detect_artifacts()`](https://mrpheus.circadia-lab.uk/reference/detect_artifacts.md):
  amplitude and muscle contamination flagging
- 📊 **Spectral analysis** —
  [`compute_band_power()`](https://mrpheus.circadia-lab.uk/reference/compute_band_power.md)
  (δ/θ/α/σ/β/γ per epoch, YASA-validated Welch PSD),
  [`compute_temporal_bandpower()`](https://mrpheus.circadia-lab.uk/reference/compute_temporal_bandpower.md)
  (sliding-window band power across the full night),
  [`compute_spectrogram()`](https://mrpheus.circadia-lab.uk/reference/compute_spectrogram.md)
- 🔁 **Sleep event detection** —
  [`compute_spindles()`](https://mrpheus.circadia-lab.uk/reference/compute_spindles.md),
  [`compute_slow_oscillations()`](https://mrpheus.circadia-lab.uk/reference/compute_slow_oscillations.md)
- 🛏️ **Automatic AASM staging** —
  [`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md):
  pre-trained LightGBM model ported from
  [YASA](https://github.com/raphaelvallat/yasa) (Vallat & Walker, 2021);
  77.8 % epoch accuracy on Sleep-EDF Cassette
- 🫁 **Respiratory metrics** —
  [`detect_apneas()`](https://mrpheus.circadia-lab.uk/reference/detect_apneas.md),
  [`compute_ahi()`](https://mrpheus.circadia-lab.uk/reference/compute_ahi.md),
  [`compute_odi()`](https://mrpheus.circadia-lab.uk/reference/compute_odi.md)
- 💓 **Cardiac rhythm** —
  [`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md)
  (Pan-Tompkins),
  [`compute_hr_signal()`](https://mrpheus.circadia-lab.uk/reference/compute_hr_signal.md),
  [`compute_hrv_sleep()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_sleep.md)
- 🏥 **Philips PMU ingestion** —
  [`read_philips_physlog()`](https://mrpheus.circadia-lab.uk/reference/read_philips_physlog.md):
  wBTU / wired / custom presets, event markers, scan-window alignment
- 🔗 **Ecosystem handoffs** —
  [`export_hypnogram()`](https://mrpheus.circadia-lab.uk/reference/export_hypnogram.md)
  → `hypnor`; derived metrics → `syncR`

------------------------------------------------------------------------

## 🗂️ Project Structure

    mrpheus/
    ├── R/
    │   ├── mrpheus-package.R           # package doc + palette_orpheus
    │   ├── read_edf.R                  # read_edf(), print.mrpheus_edf
    │   ├── prepare_psg.R               # prepare_psg(), print.mrpheus_psg
    │   ├── signal_filters.R            # remove_dc(), detect_powerline(),
    │   │                               #   notch_filter(), bandpass_filter()
    │   ├── preprocess_psg.R            # preprocess_psg()
    │   ├── correct_eog.R               # correct_eog_regression(), correct_eog_ica()
    │   ├── detect_artifacts.R          # detect_artifacts()
    │   ├── compute_band_power.R        # compute_band_power()
    │   ├── compute_temporal_bandpower.R# compute_temporal_bandpower()
    │   ├── compute_spectrogram.R       # compute_spectrogram()
    │   ├── compute_spindles.R          # compute_spindles()
    │   ├── compute_slow_oscillations.R
    │   ├── compute_hr_signal.R         # compute_hr_signal()
    │   ├── compute_hrv_sleep.R         # compute_hrv_sleep()
    │   ├── stage_epochs.R              # stage_epochs()
    │   ├── staging_features.R          # internal YASA feature extraction
    │   ├── export_hypnogram.R          # export_hypnogram()
    │   ├── detect_qrs.R                # detect_qrs(), print.mrpheus_qrs
    │   ├── detect_apneas.R             # detect_apneas(), compute_ahi(), compute_odi()
    │   └── read_philips_physlog.R      # read_philips_physlog()
    ├── src/
    │   ├── staging_features.cpp        # Rcpp hot paths: staging pipeline
    │   └── event_detection.cpp         # Rcpp hot paths: spindles, slow oscillations
    ├── inst/
    │   ├── extdata/                    # example_physlog.log
    │   ├── filters/                    # MNE FIR + resample poly coefficients
    │   └── models/                     # serialised LightGBM staging model
    ├── data-raw/                       # model extraction, parity validation scripts
    ├── tests/testthat/
    ├── vignettes/
    │   ├── getting-started.Rmd
    │   ├── mri-physiology.Rmd
    │   ├── sleep-staging.Rmd
    │   ├── eeg-preprocessing.Rmd
    │   └── spectral-analysis.Rmd
    ├── _pkgdown.yml
    └── DESCRIPTION

------------------------------------------------------------------------

## 🚀 Getting Started

### Installation

Install from [r-universe](https://circadia-bio.r-universe.dev)
(recommended — pre-built binaries, no compiler required):

``` r

install.packages(
  "mrpheus",
  repos = c("https://circadia-bio.r-universe.dev", "https://cloud.r-project.org")
)
```

Or install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("circadia-bio/mrpheus")
```

### Prerequisites

The staging model requires LightGBM. If the r-universe binary does not
install it automatically:

``` r

install.packages("lightgbm")   # CRAN binary
```

### Full PSG pipeline

``` r

library(mrpheus)

# 1. Ingest
rec <- read_edf("psg_001.edf")
psg <- prepare_psg(rec)

# 2. Preprocess
psg <- preprocess_psg(
  psg,
  channel_rename = c("C3-A2" = "C3", "C4-A1" = "C4")  # if needed
)

# 3. EOG correction
psg <- correct_eog_regression(psg)   # or correct_eog_ica()

# 4. Stage
stages <- stage_epochs(psg)

# 5. Band power
bp  <- compute_band_power(psg, relative = TRUE)
tbp <- compute_temporal_bandpower(psg, hypno = stages$stage)

# 6. Events
sp  <- compute_spindles(psg)
so  <- compute_slow_oscillations(psg)

# 7. Export to hypnor
hypnogram <- export_hypnogram(stages, start_time = rec$header$startTime)
```

------------------------------------------------------------------------

## 📦 Dependencies

| Package | Version | Purpose |
|----|----|----|
| `edfReader` | \>= 1.2.1 | EDF/EDF+ file ingestion |
| `gsignal` | \>= 0.3.5 | IIR filter design, Welch PSD |
| `pracma` | \>= 2.4.2 | Numerical utilities |
| `lightgbm` | \>= 4.0.0 | YASA staging model inference |
| `fastICA` | \>= 1.2-0 | ICA decomposition for EOG correction |
| `dplyr` | \>= 1.1.0 | Tabular manipulation |
| `tibble` | \>= 3.2.0 | Tidy output structures |
| `cli` | \>= 3.6.0 | Diagnostic messages |
| `rlang` | \>= 1.1.0 | Error handling |
| `Rcpp` | \>= 1.0.10 | C++ hot-path implementations (feature extraction, rolling statistics, event detection) |

------------------------------------------------------------------------

## 📄 Citation

If you use the automatic sleep staging feature, please also cite the
YASA paper:

``` bibtex
@article{vallat2021,
  author  = {Vallat, Raphael and Walker, Matthew P},
  title   = {An open-source, high-performance tool for automated sleep staging},
  journal = {eLife},
  volume  = {10},
  pages   = {e70092},
  year    = {2021},
  doi     = {10.7554/eLife.70092}
}
```

------------------------------------------------------------------------

## 👥 Authors

| Role | Name | Affiliation |
|----|----|----|
| Author, maintainer | Lucas França | Northumbria University, Circadia Lab |
| Author | Mario Leocadio-Miguel | Circadia Lab |

------------------------------------------------------------------------

## 🤝 Related Tools

- ⌚️ [**zeitR**](https://zeitr.circadia-lab.uk) — wrist actigraphy
  analysis and circadian metrics
- 😵‍💫 [**hypnor**](https://github.com/circadia-bio/hypnor) — hypnogram
  handling, sleep architecture metrics, and visualisation
- 🔄 [**syncR**](https://github.com/circadia-bio/syncR) — unified
  participant-indexed database (actigraphy + PSG + diary)
- 🧲 [**boldR**](https://github.com/circadia-bio/boldR) — fMRIPrep BOLD
  derivatives → parcellated analysis
- 🧮 [**tallieR**](https://github.com/circadia-bio/tallieR) —
  sociodemographics and questionnaires
- 🛌 [**slumbR**](https://github.com/circadia-bio/slumbR) — sleep diary
  processing
- 🎨 [**circadia**](https://github.com/circadia-bio/circadia) — shared
  visual identity (palettes, themes)
- 🔬 [**circadia-bio**](https://github.com/circadia-bio) — the Circadia
  Lab GitHub organisation

------------------------------------------------------------------------

## 📄 Licence

Released under the [MIT
License](https://mrpheus.circadia-lab.uk/LICENSE).

Copyright © Lucas França, Mario Leocadio-Miguel, 2025

------------------------------------------------------------------------

> **Staging model attribution:** The LightGBM model bundled in
> `inst/models/yasa_staging.txt` was originally trained as part of YASA
> (Vallat & Walker, *eLife*, 2021;
> [doi:10.7554/eLife.70092](https://doi.org/10.7554/eLife.70092)) and is
> redistributed here under YASA’s BSD 3-Clause License with attribution.
