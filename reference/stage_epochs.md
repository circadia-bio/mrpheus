# Automatic AASM sleep staging

Stages each 30-second epoch using a pre-trained LightGBM model
originally developed for YASA (Vallat & Walker, 2021) and shipped as a
cross-language serialised model in `inst/models/yasa_staging.txt`.
Features are computed in R to match the Python feature extraction
pipeline exactly; bit-exact parity is validated in the package test
suite.

## Usage

``` r
stage_epochs(
  psg,
  eeg_channel = NULL,
  eog_channel = NULL,
  emg_channel = NULL,
  artefacts = NULL,
  model_path = system.file("models", "yasa_staging.txt", package = "mrpheus")
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- eeg_channel:

  Character. Central EEG channel (e.g. `"EEG Fpz-Cz"`). If `NULL`
  (default), the first non-bad EEG channel is used.

- eog_channel:

  Character or `NULL`. EOG channel. If `NULL` (default), the first
  non-bad EOG channel is used (EOG features are omitted if none found).

- emg_channel:

  Character or `NULL`. Chin EMG channel. If `NULL` (default), the first
  non-bad EMG channel is used (EMG features are omitted if none found).

- artefacts:

  Tibble or `NULL`. Output of
  [`detect_artifacts()`](https://mrpheus.circadia-lab.uk/reference/detect_artifacts.md).
  Artefact epochs are assigned `NA` stage and excluded from the model.
  If `NULL`, all epochs are staged.

- model_path:

  Character. Path to the serialised LightGBM model. Defaults to the
  bundled model at
  `system.file("models/yasa_staging.txt", package = "mrpheus")`.

## Value

A tibble with one row per epoch:

- epoch:

  Integer.

- stage:

  Character. AASM stage: `W`, `N1`, `N2`, `N3`, `REM`, or `NA`
  (artefact).

- prob_W, prob_N1, prob_N2, prob_N3, prob_REM:

  Numeric. Posterior class probabilities from the LightGBM model.

## Details

Stages returned follow standard AASM nomenclature: `W` (wake), `N1`,
`N2`, `N3`, `REM`.

## References

Vallat, R., & Walker, M. P. (2021). An open-source, high-performance
tool for automated sleep staging. *eLife*, 10, e70092.
[doi:10.7554/eLife.70092](https://doi.org/10.7554/eLife.70092)

## See also

[`export_hypnogram()`](https://mrpheus.circadia-lab.uk/reference/export_hypnogram.md)
