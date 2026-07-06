# Remove EOG artefacts using ICA

Decomposes the EEG into independent components via FastICA, identifies
components whose time courses correlate with the EOG channels above
`threshold`, subtracts their contribution from the EEG, and returns a
cleaned `mrpheus_psg`. Mirrors MNE's ICA EOG rejection pipeline with
default settings matching the companion Python notebook
(`n_components = min(6, n_eeg)`, `threshold = 0.35`).

## Usage

``` r
correct_eog_ica(
  psg,
  eog_channels = NULL,
  eeg_channels = NULL,
  n_components = NULL,
  threshold = 0.35,
  fun = "logcosh",
  verbose = TRUE
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- eog_channels:

  Character vector. EOG channel labels used as the artifact reference.
  If `NULL` (default), all `"EOG"` channels are used.

- eeg_channels:

  Character vector. EEG channels to decompose and clean. If `NULL`
  (default), all non-bad `"EEG"` channels are used.

- n_components:

  Integer or `NULL`. Number of ICA components. `NULL` (default) uses
  `min(6L, n_eeg_channels)`, matching the notebook.

- threshold:

  Numeric. Absolute Pearson correlation threshold above which a
  component is flagged as EOG-related. Default `0.35`.

- fun:

  Character. Contrast function passed to
  [`fastICA::fastICA()`](https://rdrr.io/pkg/fastICA/man/fastICA.html).
  `"logcosh"` (default) or `"exp"`. The R method is used internally for
  robustness across channel counts; this is slower than the C backend
  but avoids matrix-conformality errors on small channel sets.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A new `mrpheus_psg` with cleaned EEG signals and re-segmented epochs.
EOG and other channels are unchanged.

## Details

**Note:** FastICA (this implementation) and MNE's default InfoMax ICA
find different decompositions; the identified components and final
signal will therefore differ from MNE output, but artifact removal
performance is comparable.

## Examples

``` r
if (FALSE) { # \dontrun{
rec   <- read_edf("psg.edf")
psg   <- prepare_psg(rec) |> preprocess_psg()
clean <- correct_eog_ica(psg)

# Stricter threshold
clean <- correct_eog_ica(psg, threshold = 0.5)
} # }
```
