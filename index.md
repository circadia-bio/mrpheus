# 🪉 mrpheus

**Raw polysomnography signal analysis for sleep and circadian
research.**

[![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3)](https://www.r-project.org/)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://mrpheus.circadia-lab.uk/LICENSE)
[![Circadia
Lab](https://img.shields.io/badge/circadia--bio-GitHub-014370)](https://github.com/circadia-bio)

------------------------------------------------------------------------

## 📖 What is mrpheus?

`mrpheus` is the deepest layer of the [Circadia
Lab](https://github.com/circadia-bio) R ecosystem — it ingests raw
multi-channel PSG recordings (EDF/EDF+) and produces the two outputs
that feed everything downstream: staged hypnograms (passed to `hypnor`)
and PSG-derived metrics (passed to `syncR`).

The name is spelled as *Morpheus* but carries a silent **m**, so it is
pronounced as *Orpheus* — a portmanteau of both myths simultaneously.
Morpheus, the god of dreams, gives the package its subject matter.
Orpheus, who descended into the underworld to navigate the unconscious,
gives it its spirit. Sleep is both.

> *“He calmed the savage beasts with song, as mrpheus calms your raw
> signals into stages.”*

### The Orpheus Mosaic Palette

The package ships a data object `palette_orpheus` — 8 colours extracted
from the [Roman mosaic *Orpheus Charming the
Animals*](https://en.wikipedia.org/wiki/Orpheus) (3rd century AD,
Palermo Archaeological Museum), the image that inspired the package
name:

| Name       | Hex       | Source                  |
|------------|-----------|-------------------------|
| sand       | `#CDB992` | warm background tessera |
| vermillion | `#B83E2C` | Orpheus’s robe          |
| olive      | `#6A7840` | tree and vegetation     |
| umber      | `#7C5432` | animal fur and earth    |
| bistre     | `#3C2212` | shadows and outlines    |
| ochre      | `#B07C3A` | warm amber accent       |
| slate      | `#6C8284` | birds; dusty teal-grey  |
| ivory      | `#EAD6AA` | highlights              |

``` r

scales::show_col(mrpheus::palette_orpheus)
```

------------------------------------------------------------------------

## ✨ Features

- 📂 **EDF/EDF+ ingestion** —
  [`read_edf()`](https://mrpheus.circadia-lab.uk/reference/read_edf.md)
  /
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md):
  channel inventory, epoch segmentation, bad-channel flagging
- 🚫 **Artefact detection** —
  [`detect_artifacts()`](https://mrpheus.circadia-lab.uk/reference/detect_artifacts.md):
  amplitude and muscle contamination flagging
- 📊 **Spectral analysis** —
  [`compute_band_power()`](https://mrpheus.circadia-lab.uk/reference/compute_band_power.md)
  (δ/θ/α/σ/β/γ per epoch),
  [`compute_spectrogram()`](https://mrpheus.circadia-lab.uk/reference/compute_spectrogram.md)
- 🔁 **Sleep event detection** —
  [`compute_spindles()`](https://mrpheus.circadia-lab.uk/reference/compute_spindles.md),
  [`compute_slow_oscillations()`](https://mrpheus.circadia-lab.uk/reference/compute_slow_oscillations.md)
- 🛏️ **Automatic AASM staging** —
  [`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md):
  uses a pre-trained LightGBM model ported from
  [YASA](https://github.com/raphaelvallat/yasa) (Vallat & Walker, 2021);
  see [Staging model](#staging-model) below
- 🫁 **Respiratory metrics** —
  [`detect_apneas()`](https://mrpheus.circadia-lab.uk/reference/detect_apneas.md),
  [`compute_ahi()`](https://mrpheus.circadia-lab.uk/reference/compute_ahi.md),
  [`compute_odi()`](https://mrpheus.circadia-lab.uk/reference/compute_odi.md)
- 💓 **Cardiac HRV** —
  [`compute_hrv_sleep()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_sleep.md):
  HRV stratified by sleep stage
- 🔗 **Ecosystem handoffs** —
  [`export_hypnogram()`](https://mrpheus.circadia-lab.uk/reference/export_hypnogram.md)
  → `hypnor`; PSG metrics → `syncR`

------------------------------------------------------------------------

## 🗂️ Project Structure

    mrpheus/
    ├── R/
    │   ├── mrpheus-package.R       # package doc + palette_orpheus data doc
    │   ├── read_edf.R              # read_edf(), print.mrpheus_edf
    │   ├── prepare_psg.R           # prepare_psg(), print.mrpheus_psg
    │   ├── detect_artifacts.R      # detect_artifacts()
    │   ├── compute_band_power.R    # compute_band_power()
    │   ├── compute_spectrogram.R   # compute_spectrogram()
    │   ├── compute_spindles.R      # compute_spindles()
    │   ├── compute_slow_oscillations.R
    │   ├── stage_epochs.R          # stage_epochs(), .extract_staging_features()
    │   ├── export_hypnogram.R      # export_hypnogram()
    │   ├── detect_apneas.R         # detect_apneas(), compute_ahi(), compute_odi()
    │   └── compute_hrv_sleep.R     # compute_hrv_sleep()
    ├── inst/models/
    │   └── yasa_staging.txt        # serialised LightGBM model (see data-raw/)
    ├── data-raw/
    │   ├── fetch_yasa_model.py     # extracts model from YASA Python package
    │   └── palette_orpheus.R       # generates palette_orpheus data object
    ├── tests/testthat/
    ├── _pkgdown.yml
    └── DESCRIPTION

------------------------------------------------------------------------

## 🚀 Getting Started

### Prerequisites

``` r

# R >= 4.1
install.packages(c("edfReader", "gsignal", "pracma", "lightgbm",
                   "dplyr", "tibble", "cli", "rlang"))
```

### Installation

``` r

remotes::install_github("circadia-bio/mrpheus")
```

### Staging model

The LightGBM staging model is not bundled in the repo by default (it
requires Python + YASA to extract). Run once after cloning:

``` bash
pip install yasa lightgbm
python data-raw/fetch_yasa_model.py
# then: usethis::use_data() or copy to inst/models/yasa_staging.txt
```

The model weights are cross-language — a model serialised by Python’s
`lgb.Booster.save_model()` loads identically in R’s
[`lightgbm::lgb.load()`](https://rdrr.io/pkg/lightgbm/man/lgb.load.html).

### Basic pipeline

``` r

library(mrpheus)

# 1. Ingest
rec <- read_edf("psg_001.edf")
psg <- prepare_psg(rec)

# 2. Artefacts
art <- detect_artifacts(psg)

# 3. Band power
bp  <- compute_band_power(psg, relative = TRUE)

# 4. Events
sp  <- compute_spindles(psg)
so  <- compute_slow_oscillations(psg)

# 5. Stage
staging <- stage_epochs(psg, artefacts = art)

# 6. Export to hypnor
hypnogram <- export_hypnogram(staging, start_time = rec$header$startTime)
```

------------------------------------------------------------------------

## 📦 Dependencies

| Package     | Version | Purpose                       |
|-------------|---------|-------------------------------|
| `edfReader` | ≥ 1.2.1 | EDF/EDF+ file ingestion       |
| `gsignal`   | ≥ 0.3.5 | Digital filtering, Welch PSD  |
| `pracma`    | ≥ 2.4.2 | Numerical integration (trapz) |
| `lightgbm`  | ≥ 4.0.0 | YASA staging model inference  |
| `dplyr`     | ≥ 1.1.0 | Tabular manipulation          |
| `tibble`    | ≥ 3.2.0 | Tidy output structures        |
| `cli`       | ≥ 3.6.0 | Diagnostic messages           |
| `rlang`     | ≥ 1.1.0 | Error handling                |

------------------------------------------------------------------------

## 👥 Authors

| Role | Name | Affiliation |
|----|----|----|
| Author, maintainer | Lucas França | Northumbria University, Circadia Lab |
| Author | Mario Leocadio-Miguel | Circadia Lab |

------------------------------------------------------------------------

## 🤝 Related Tools

- 🌊 [**zeitR**](https://github.com/circadia-bio/zeitR) — wrist
  actigraphy analysis and circadian metrics
- 😴 [**hypnor**](https://github.com/circadia-bio/hypnor) — hypnogram
  handling, plotting, and architecture metrics
- 🔗 [**syncR**](https://github.com/circadia-bio/syncR) — unified
  participant-indexed database (actigraphy + sleep diary + PSG)
- 📋 [**tallieR**](https://github.com/circadia-bio/tallieR) —
  sociodemographics and questionnaires
- 📓 [**slumbR**](https://github.com/circadia-bio/slumbR) — sleep diary
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
