# Getting started with mrpheus

## What is mrpheus?

**mrpheus** is an R package for raw polysomnography (PSG) signal
analysis. It ingests EDF/EDF+ recordings, detects sleep events, computes
spectral features, and stages each 30-second epoch using a pre-trained
LightGBM model ported from [YASA](https://github.com/raphaelvallat/yasa)
(Vallat & Walker, 2021).

The name is spelled as *Morpheus* but carries a silent **m**, pronounced
as *Orpheus* — carrying both myths at once. Morpheus, god of dreams,
gives the package its subject; Orpheus, who descended into the
underworld to navigate the unconscious, gives it its spirit.

mrpheus is part of the [Circadia Lab](https://circadia-lab.uk) R
ecosystem:

| Package | Purpose |
|----|----|
| **mrpheus** | Raw PSG signal analysis and automatic sleep staging |
| [hypnor](https://github.com/circadia-bio/hypnor) | Hypnogram handling and sleep architecture metrics |
| [zeitR](https://zeitr.circadia-lab.uk) | Wrist actigraphy analysis and circadian metrics |
| [syncR](https://github.com/circadia-bio/syncR) | Unified participant-indexed database |
| [slumbR](https://github.com/circadia-bio/slumbR) | Sleep diary processing |

------------------------------------------------------------------------

## Installation

``` r

# install.packages("pak")
pak::pak("circadia-bio/mrpheus")
```

### Staging model

The automatic sleep staging model is not bundled in the package by
default — it must be extracted from YASA once after installation:

``` bash
pip install yasa lightgbm
python data-raw/fetch_yasa_model.py
```

This saves `inst/models/yasa_staging.txt` — a cross-language LightGBM
model that loads identically in R. See
[`?stage_epochs`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
for details.

------------------------------------------------------------------------

## Reading a recording

The entry point is
[`read_edf()`](https://mrpheus.circadia-lab.uk/reference/read_edf.md),
which wraps `edfReader` with consistent output formatting:

``` r

library(mrpheus)

rec <- read_edf("recordings/psg_001.edf")
rec
```

    ## mrpheus EDF recording
    ## ℹ Path: recordings/psg_001.edf
    ## ℹ Duration: 7.98 hours
    ## ℹ Channels: 14
    ## # A tibble: 14 × 3
    ##   label              sample_rate transducer
    ##   <chr>                    <dbl> <chr>
    ## 1 EEG Fpz-Cz                 100 AgAgCl electrode
    ## 2 EEG Pz-Oz                  100 AgAgCl electrode
    ## 3 EOG horizontal              100 AgAgCl electrode
    ## ...

Pass `only_header = TRUE` to inspect channels without loading signal
data:

``` r

hdr <- read_edf("recordings/psg_001.edf", only_header = TRUE)
```

------------------------------------------------------------------------

## Preparing a recording

[`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md)
segments the recording into 30-second epochs, classifies channels by
type (EEG, EOG, EMG, ECG, RESP), and flags flat or bad channels:

``` r

psg <- prepare_psg(rec)
psg
```

    ## mrpheus PSG
    ## ℹ Epochs: 958 × 30 s
    ## ℹ Recording: 7.98 h
    ## # A tibble: 14 × 4
    ##   label              type  sample_rate   bad
    ##   <chr>              <chr>       <dbl> <lgl>
    ## 1 EEG Fpz-Cz         EEG           100 FALSE
    ## 2 EEG Pz-Oz          EEG           100 FALSE
    ## 3 EOG horizontal      EOG           100 FALSE
    ## ...

Channel classification uses regex patterns which can be overridden:

``` r

psg <- prepare_psg(rec,
  eeg_pattern  = "EEG|C3|C4|F3|F4|O1|O2|Fpz|Pz",
  resp_pattern = "Thor|Abdo|Flow|SpO2|airflow"
)
```

------------------------------------------------------------------------

## Artefact detection

[`detect_artifacts()`](https://mrpheus.circadia-lab.uk/reference/detect_artifacts.md)
flags epochs with excessive amplitude or high-frequency (muscle)
contamination:

``` r

art <- detect_artifacts(psg)
```

    ## ℹ Artefact epochs: 12 / 958 (1.3%)

``` r

# Inspect flagged epochs
art[art$artefact, ]
```

    ## # A tibble: 12 × 5
    ##   epoch artefact reason    peak_to_peak_uv hf_power_db
    ##   <int> <lgl>    <chr>               <dbl>       <dbl>
    ## 1    34 TRUE     amplitude            524.        -18.2
    ## 2    87 TRUE     high_freq            312.          6.4
    ## ...

------------------------------------------------------------------------

## Spectral analysis

### Band power

[`compute_band_power()`](https://mrpheus.circadia-lab.uk/reference/compute_band_power.md)
estimates PSD via Welch’s method and integrates within standard EEG
bands per epoch per channel:

``` r

bp <- compute_band_power(psg, relative = TRUE)
bp
```

    ## # A tibble: 958 × 9
    ##   epoch channel    delta theta alpha sigma  beta gamma total_power
    ##   <int> <chr>      <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>       <dbl>
    ## 1     1 EEG Fpz-Cz 0.412 0.183 0.142 0.089 0.112 0.062          NA
    ## ...

### Spectrogram

[`compute_spectrogram()`](https://mrpheus.circadia-lab.uk/reference/compute_spectrogram.md)
returns a time × frequency power matrix:

``` r

sg <- compute_spectrogram(psg, channel = "EEG Fpz-Cz", db = TRUE)
sg
```

    ## mrpheus spectrogram
    ## ℹ Channel: EEG Fpz-Cz
    ## ℹ Epochs: 958
    ## ℹ Frequency range: 0-40 Hz (161 bins)
    ## ℹ Units: dB

------------------------------------------------------------------------

## Sleep event detection

### Spindles

[`compute_spindles()`](https://mrpheus.circadia-lab.uk/reference/compute_spindles.md)
detects sleep spindles via RMS envelope thresholding in the sigma band
(11–16 Hz):

``` r

sp <- compute_spindles(psg, channel = "EEG Fpz-Cz")
```

    ## ✔ Detected 312 spindles in channel 'EEG Fpz-Cz'.

``` r

sp
```

    ## # A tibble: 312 × 7
    ##   epoch start_s end_s duration_s peak_freq_hz rms_uv channel
    ##   <int>   <dbl> <dbl>      <dbl>        <dbl>  <dbl> <chr>
    ## 1   142    4.23  5.11       0.88         13.2    18.4 EEG Fpz-Cz
    ## ...

### Slow oscillations

[`compute_slow_oscillations()`](https://mrpheus.circadia-lab.uk/reference/compute_slow_oscillations.md)
uses zero-crossing detection in the delta band (0.5–2 Hz), following
Mölle et al. (2002):

``` r

so <- compute_slow_oscillations(psg, channel = "EEG Fpz-Cz")
```

    ## ✔ Detected 188 slow oscillations in channel 'EEG Fpz-Cz'.

------------------------------------------------------------------------

## Automatic sleep staging

[`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
stages each epoch using the pre-trained LightGBM model. Pass the
artefact table to exclude flagged epochs:

``` r

staging <- stage_epochs(
  psg,
  eeg_channel = "EEG Fpz-Cz",
  eog_channel = "EOG horizontal",
  emg_channel = "EMG chin",
  artefacts   = art
)
staging
```

    ## # A tibble: 958 × 7
    ##   epoch stage prob_W prob_N1 prob_N2 prob_N3 prob_REM
    ##   <int> <chr>  <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
    ## 1     1 W      0.921   0.042   0.022   0.008    0.007
    ## 2     2 W      0.874   0.081   0.031   0.007    0.007
    ## 3     3 N1     0.312   0.448   0.185   0.032    0.023
    ## ...

Stages follow AASM nomenclature: `W`, `N1`, `N2`, `N3`, `REM`. Artefact
epochs are coded `NA`.

> **Note:**
> [`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
> requires the staging model at `inst/models/yasa_staging.txt`. Feature
> parity with YASA’s Python pipeline is still being validated — see the
> source of `.extract_staging_features()` for the current implementation
> status.

------------------------------------------------------------------------

## Exporting to hypnor

[`export_hypnogram()`](https://mrpheus.circadia-lab.uk/reference/export_hypnogram.md)
converts the staging output into a `hypnor_hypnogram` object for
downstream architecture analysis:

``` r

hypnogram <- export_hypnogram(
  staging,
  epoch_s        = 30,
  start_time     = rec$header$startTime,
  participant_id = "P001"
)

# hypnor functions now work directly on the hypnogram
hypnor::plot_hypnogram(hypnogram)
hypnor::compute_architecture(hypnogram)
```

------------------------------------------------------------------------

## Full pipeline

Putting it all together:

``` r

library(mrpheus)

# 1. Ingest
rec      <- read_edf("recordings/psg_001.edf")
psg      <- prepare_psg(rec)

# 2. Quality
art      <- detect_artifacts(psg)

# 3. Spectral
bp       <- compute_band_power(psg, relative = TRUE)
sg       <- compute_spectrogram(psg, channel = "EEG Fpz-Cz")

# 4. Events
sp       <- compute_spindles(psg)
so       <- compute_slow_oscillations(psg)

# 5. Stage
staging  <- stage_epochs(psg, artefacts = art)

# 6. Hand off
hypnogram <- export_hypnogram(staging, start_time = rec$header$startTime)
```

------------------------------------------------------------------------

## References

Vallat, R., & Walker, M. P. (2021). An open-source, high-performance
tool for automated sleep staging. *eLife*, 10, e70092.
<https://doi.org/10.7554/eLife.70092>

Mölle, M., Marshall, L., Gais, S., & Born, J. (2002). Grouping of
spindle activity during slow oscillations in human non-rapid eye
movement sleep. *Journal of Neuroscience*, 22(24), 10941–10947.
