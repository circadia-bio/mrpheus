# Preprocess a PSG recording

Applies the standard sleep EEG preprocessing pipeline to an
`mrpheus_psg` object: optional channel renaming, DC offset removal,
powerline notch filtering (plus harmonics), and channel-type-specific
bandpass filtering.

## Usage

``` r
preprocess_psg(
  psg,
  channel_rename = NULL,
  dc = TRUE,
  powerline_freq = NULL,
  notch_harmonics = TRUE,
  notch_bw_hz = 2,
  eeg_bandpass = c(0.3, 35),
  eog_bandpass = c(0.3, 15),
  emg_bandpass = c(10, 99),
  ecg_bandpass = c(0.5, 40),
  verbose = TRUE
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- channel_rename:

  Named character vector mapping old channel labels to new labels, e.g.
  `c("C3-A2" = "C3", "O2-A1" = "O2")`. Applied before filtering; updates
  `channel_map` labels accordingly. Default `NULL`.

- dc:

  Logical. Remove per-channel DC offset via
  [`remove_dc()`](https://mrpheus.circadia-lab.uk/reference/remove_dc.md).
  Default `TRUE`.

- powerline_freq:

  Integer or `NULL`. Powerline frequency in Hz (`50L` or `60L`). `NULL`
  (default) triggers auto-detection via
  [`detect_powerline()`](https://mrpheus.circadia-lab.uk/reference/detect_powerline.md)
  on the first clean EEG channel.

- notch_harmonics:

  Logical. Notch harmonics of `powerline_freq` up to Nyquist - 5 Hz.
  Default `TRUE`.

- notch_bw_hz:

  Numeric. Full bandwidth of each notch in Hz. Default `2`.

- eeg_bandpass:

  Numeric vector of length 2. Bandpass limits (Hz) applied to all EEG
  channels. Default `c(0.3, 35)`.

- eog_bandpass:

  Numeric vector of length 2. Bandpass limits (Hz) applied to all EOG
  channels. Default `c(0.3, 15)`.

- emg_bandpass:

  Numeric vector of length 2. Bandpass limits (Hz) applied to all EMG
  channels. Default `c(10, 99)`.

- ecg_bandpass:

  Numeric vector of length 2. Bandpass limits (Hz) applied to all ECG
  channels. Default `c(0.5, 40)`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A new `mrpheus_psg` object with filtered signals and re-segmented
epochs. The `edf$signals` entries are replaced in-place with the
filtered continuous signals so that subsequent epoch access is
consistent.

## Details

Filtering is applied to the full **continuous** signal stored in
`psg$edf$signals` before re-epoching, which avoids filter
discontinuities at epoch boundaries. The returned object is a new
`mrpheus_psg` with the same epoch length and channel map as the input.

The individual filter steps
([`remove_dc()`](https://mrpheus.circadia-lab.uk/reference/remove_dc.md),
[`detect_powerline()`](https://mrpheus.circadia-lab.uk/reference/detect_powerline.md),
[`notch_filter()`](https://mrpheus.circadia-lab.uk/reference/notch_filter.md),
[`bandpass_filter()`](https://mrpheus.circadia-lab.uk/reference/bandpass_filter.md))
are exported separately and can be called on any numeric vector without
constructing a PSG object.

## Examples

``` r
if (FALSE) { # \dontrun{
rec  <- read_edf("data/psg_001.edf")
psg  <- prepare_psg(rec)

# Auto-detect powerline, apply standard filters
psg_clean <- preprocess_psg(psg)

# Rename linked-ear channels to 10-20 names, fix powerline explicitly
psg_clean <- preprocess_psg(
  psg,
  channel_rename = c("C3-A2" = "C3", "C4-A1" = "C4",
                     "O1-A2" = "O1", "O2-A1" = "O2"),
  powerline_freq = 50L
)
} # }
```
