# Automatic AASM sleep staging

Stages each 30-second epoch with the pre-trained LightGBM classifiers
that ship with YASA (Vallat & Walker, 2021), converted to LightGBM text
files in `inst/models/`. As in YASA, the classifier is chosen from the
channels that are supplied: EEG only, EEG + EOG, EEG + EMG or EEG +
EOG + EMG, each with or without the demographic features (`age`,
`male`). Features are computed in R to match the Python feature
extraction and are rounded to single precision before prediction, as
YASA does. See the `NEWS.md` entry for the agreement that was measured
against YASA and for its limits.

Stages returned follow standard AASM nomenclature: `W` (wake), `N1`,
`N2`, `N3`, `REM`.

## Usage

``` r
stage_epochs(
  psg,
  eeg_channel = NULL,
  eog_channel = NULL,
  emg_channel = NULL,
  artefacts = NULL,
  model_path = NULL,
  age = NULL,
  male = NULL
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

  Character, `NA` or `NULL`. EOG channel. If `NULL` (default), the first
  non-bad EOG channel is used. Without an EOG channel (`NA`, or none
  found) the EEG-only or EEG + EMG classifier is used.

- emg_channel:

  Character, `NA` or `NULL`. Chin EMG channel. If `NULL` (default), the
  first non-bad EMG channel is used. Without an EMG channel (`NA`, or
  none found) the classifier without EMG is used.

- artefacts:

  Tibble or `NULL`. Output of
  [`detect_artifacts()`](https://mrpheus.circadia-lab.uk/reference/detect_artifacts.md).
  Artefact epochs are assigned `NA` stage. The classifier still sees
  every epoch, so the rolling-window features of neighbouring epochs are
  unchanged. If `NULL`, all epochs are staged.

- model_path:

  Character or `NULL`. Path to a serialised LightGBM model (plain text
  or `.gz`), for example one of your own. `NULL` (default) picks the
  bundled YASA classifier that matches the channels and `age` / `male`.

- age:

  Numeric or `NULL`. Age in years, between 0 and 120. Together with
  `male` it selects the `+demo` classifiers. Truncated to an integer, as
  in YASA.

- male:

  Logical, `0` or `1`, or `NULL`. Sex (`1` / `TRUE` for male). Must be
  supplied together with `age`.

## Value

A tibble with one row per epoch:

- epoch:

  Integer.

- stage:

  Character. AASM stage: `W`, `N1`, `N2`, `N3`, `REM`, or `NA`
  (artefact).

- prob_W, prob_N1, prob_N2, prob_N3, prob_REM:

  Numeric. Posterior class probabilities from the LightGBM model.

The classifier that was used is attached as the attribute
`"classifier"`.

## References

Vallat, R., & Walker, M. P. (2021). An open-source, high-performance
tool for automated sleep staging. *eLife*, 10, e70092.
[doi:10.7554/eLife.70092](https://doi.org/10.7554/eLife.70092)

## See also

[`export_hypnogram()`](https://mrpheus.circadia-lab.uk/reference/export_hypnogram.md)
