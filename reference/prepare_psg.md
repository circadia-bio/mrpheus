# Prepare a PSG recording for analysis

Takes an `mrpheus_edf` object and segments it into standard epochs,
performs a channel inventory (classifying signals by type), and flags
channels that appear flat or likely bad. This is the standard entry
point before any downstream analyses (spectral, event detection,
staging).

## Usage

``` r
prepare_psg(
  edf,
  epoch_s = 30,
  eeg_pattern = "EEG|C3|C4|F3|F4|O1|O2|Fpz|Pz",
  eog_pattern = "EOG|ROC|LOC",
  emg_pattern = "EMG|chin|Chin",
  ecg_pattern = "ECG|EKG",
  resp_pattern = "Thor|Abdo|Flow|SpO2|airflow",
  flat_threshold = 1e-06
)
```

## Arguments

- edf:

  An `mrpheus_edf` object from
  [`read_edf()`](https://mrpheus.circadia-lab.uk/reference/read_edf.md).

- epoch_s:

  Numeric. Epoch length in seconds. Default `30` (standard AASM epoch).

- eeg_pattern:

  Character. Regex pattern to identify EEG channels. Default
  `"EEG|C3|C4|F3|F4|O1|O2|Fpz|Pz"`.

- eog_pattern:

  Character. Regex pattern to identify EOG channels. Default
  `"EOG|ROC|LOC"`.

- emg_pattern:

  Character. Regex pattern to identify EMG channels. Default
  `"EMG|chin|Chin"`.

- ecg_pattern:

  Character. Regex pattern to identify ECG/EKG channels. Default
  `"ECG|EKG"`.

- resp_pattern:

  Character. Regex pattern to identify respiratory channels. Default
  `"Thor|Abdo|Flow|SpO2|airflow"`.

- flat_threshold:

  Numeric. Variance below this value flags a channel as flat/bad.
  Default `1e-6`.

## Value

A list of class `mrpheus_psg` with components:

- edf:

  The original `mrpheus_edf` object.

- epochs:

  List. Each element is one epoch (30 s by default), itself a named list
  of channel vectors.

- n_epochs:

  Integer. Total number of complete epochs.

- epoch_s:

  Numeric. Epoch duration in seconds.

- channel_map:

  Data frame. Channel label, detected type (EEG/EOG/EMG/ECG/RESP/OTHER),
  sample rate, and `bad` flag.

## Examples

``` r
if (FALSE) { # \dontrun{
rec  <- read_edf("data/psg_001.edf")
psg  <- prepare_psg(rec)
psg$channel_map
} # }
```
