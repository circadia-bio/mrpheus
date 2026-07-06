# Remove EOG artefacts using linear regression

For each EEG channel, fits a multiple linear regression on all EOG
reference channels and subtracts the fitted EOG component, leaving the
residual as the cleaned signal. This removes both transient blink
artefacts and slow ocular drift without requiring blink event detection.

## Usage

``` r
correct_eog_regression(
  psg,
  eog_channels = NULL,
  eeg_channels = NULL,
  verbose = TRUE
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- eog_channels:

  Character vector. EOG channel labels to use as regressors. If `NULL`
  (default), all channels of type `"EOG"` in `psg$channel_map` are used.

- eeg_channels:

  Character vector. EEG channel labels to clean. If `NULL` (default),
  all non-bad `"EEG"` channels are used.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A new `mrpheus_psg` with cleaned EEG signals and re-segmented epochs.
EOG and other channels are unchanged.

## Details

This is the pure-R equivalent of MNE's `compute_proj_eog()` / SSP
projection pipeline. Both approaches project out the subspace spanned by
the EOG channels; regression does so directly on the continuous signal.

## Examples

``` r
if (FALSE) { # \dontrun{
rec   <- read_edf("psg.edf")
psg   <- prepare_psg(rec) |> preprocess_psg()
clean <- correct_eog_regression(psg)
} # }
```
