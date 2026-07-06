# PSG Preprocessing Pipeline

**mrpheus** provides a full PSG preprocessing pipeline through
[`preprocess_psg()`](https://mrpheus.circadia-lab.uk/reference/preprocess_psg.md).
This vignette walks through the standard workflow for a Sleep-EDF
recording, then shows how to use the individual filter functions
independently when you need finer control.

> **Note** — All code blocks use `eval = FALSE` because they require a
> local EDF file. For a fully-executed walkthrough see the companion
> article: [Preprocessing: Live
> Walkthrough](https://mrpheus.circadia-lab.uk/articles/preprocessing-demo.html)

------------------------------------------------------------------------

## Data

The examples use **SC4001E0-PSG.edf** from the Sleep-EDF Cassette study
(Kemp et al., 2000), a 22-hour ambulatory PSG of a healthy young adult.
Freely available from
[PhysioNet](https://physionet.org/content/sleep-edfx/1.0.0/) under CC BY
1.0 (free account required).

``` python
# Download via MNE Python (handles authentication automatically)
import mne
files = mne.datasets.sleep_physionet.age.fetch_data(subjects=[0], recording=[1])
print(files[0])  # path to SC4001E0-PSG.edf
```

------------------------------------------------------------------------

## Standard pipeline

The typical workflow is four steps:

``` r

library(mrpheus)

# 1. Read EDF
rec <- read_edf("SC4001E0-PSG.edf")

# 2. Epoch and classify channels
psg <- prepare_psg(rec)
psg$channel_map
```

The Sleep-EDF file uses linked-ear channel names (`EEG Fpz-Cz`,
`EEG Pz-Oz`). These are already matched by
[`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md)’s
default EEG pattern, but if your recording uses notation like `C3-A2`
you can pass a rename map at the preprocessing step:

``` r

# 3. Preprocess: rename, remove DC, notch, bandpass
psg_clean <- preprocess_psg(
  psg,
  channel_rename = c(
    "C3-A2" = "C3", "C4-A1" = "C4",
    "O1-A2" = "O1", "O2-A1" = "O2"
  )
)
```

[`preprocess_psg()`](https://mrpheus.circadia-lab.uk/reference/preprocess_psg.md)
does the following in order, on the **full continuous signal** before
re-epoching (this avoids filter discontinuities at epoch boundaries):

| Step | Default behaviour |
|----|----|
| Channel renaming | Optional; updates `channel_map` and all internal references |
| DC removal | Subtracts per-channel mean |
| Powerline detection | Compares PSD power at 50 and 60 Hz; picks the dominant one |
| Notch filter | 2nd-order Butterworth bandstop at powerline frequency + harmonics |
| Bandpass filter | EEG 0.3–35 Hz · EOG 0.3–15 Hz · EMG 10–99 Hz · ECG 0.5–40 Hz |

``` r

# 4. Automatic AASM sleep staging on the cleaned PSG
stages <- stage_epochs(psg_clean)
head(stages)
```

------------------------------------------------------------------------

## Customising the pipeline

All parameters have sensible defaults but can be overridden:

``` r

psg_clean <- preprocess_psg(
  psg,
  powerline_freq  = 50L,          # skip auto-detection
  notch_harmonics = TRUE,         # also notch 100, 150 Hz, …
  notch_bw_hz     = 3,            # wider notch bandwidth
  eeg_bandpass    = c(0.5, 40),   # broader EEG range
  eog_bandpass    = c(0.3, 15),
  emg_bandpass    = c(10, 99),
  ecg_bandpass    = c(0.5, 40),
  dc              = TRUE,
  verbose         = TRUE
)
```

To skip channel-type-specific filtering entirely, set a bandpass to the
full range for that type:

``` r

# Keep raw EMG (no bandpass)
psg_clean <- preprocess_psg(psg, emg_bandpass = c(0, Inf))
```

------------------------------------------------------------------------

## Using the signal-level functions independently

Each step in
[`preprocess_psg()`](https://mrpheus.circadia-lab.uk/reference/preprocess_psg.md)
delegates to an exported function that works on a plain numeric vector.
These are useful when:

- You want to filter a single channel without building a full PSG
  object.
- You are building your own pipeline outside the `mrpheus_psg`
  structure.
- You need to apply filtering to signals loaded from another source
  (e.g. CSV, MATLAB `.mat`, or a custom reader).

### DC removal

``` r

eeg <- rec$signals[["EEG Fpz-Cz"]]$signal
eeg <- remove_dc(eeg)
```

### Powerline detection

``` r

sr   <- rec$channels$sample_rate[rec$channels$label == "EEG Fpz-Cz"]
freq <- detect_powerline(eeg, sr)
#> ℹ Detected powerline frequency: 50 Hz
```

### Notch filtering

``` r

# Remove 50 Hz and harmonics (100, 150 Hz, …)
eeg_notched <- notch_filter(eeg, sr = sr, freq = 50)

# Remove 60 Hz only, no harmonics
eeg_notched <- notch_filter(eeg, sr = sr, freq = 60, harmonics = FALSE)

# Wider notch bandwidth (default is 2 Hz)
eeg_notched <- notch_filter(eeg, sr = sr, freq = 50, bandwidth_hz = 4)
```

### Bandpass filtering

``` r

# Standard EEG passband
eeg_bp <- bandpass_filter(eeg, sr = sr, low_hz = 0.3, high_hz = 35)

# High-pass only (baseline drift removal) — set high_hz beyond Nyquist
eeg_hp <- bandpass_filter(eeg, sr = sr, low_hz = 0.1, high_hz = Inf)

# Low-pass only — set low_hz to 0
eeg_lp <- bandpass_filter(eeg, sr = sr, low_hz = 0, high_hz = 35)
```

------------------------------------------------------------------------

## A note on filter design

The notch and bandpass filters both use **zero-phase IIR (Butterworth)**
designs applied via
[`gsignal::filtfilt()`](https://rdrr.io/pkg/gsignal/man/filtfilt.html):

- Notch: 2nd-order Butterworth bandstop (≈9 coefficients). Safe on long
  recordings — no memory issues.
- Bandpass: 4th-order Butterworth. The `high_hz` cutoff is automatically
  clamped to 99 % of the Nyquist frequency to prevent instability.

Zero-phase (`filtfilt`) filtering introduces no phase distortion but
doubles the effective filter order. If you prefer a causal (single-pass)
filter for any reason, use
[`gsignal::filter()`](https://rdrr.io/pkg/gsignal/man/filter.html)
directly with the same Butterworth coefficients.

The staging pipeline
([`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md))
applies its own separate MNE-matched FIR bandpass internally —
[`preprocess_psg()`](https://mrpheus.circadia-lab.uk/reference/preprocess_psg.md)
filtering and staging filtering are independent.

------------------------------------------------------------------------

## References

Kemp B, Zwinderman A, Tuk B, Kamphuisen H, Oberyé J (2000). Analysis of
a sleep-dependent neuronal feedback loop: the slow-wave microcontinuity
of the EEG. *IEEE Transactions on Biomedical Engineering*, 47(9),
1185–1194.
