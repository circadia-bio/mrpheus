# Compute temporal band power across a PSG recording

Applies a sliding window over the scored epochs, concatenates the signal
within each window, estimates the Welch PSD averaged across channels,
and integrates power within standard EEG frequency bands. Returns a
long-format tibble suitable for plotting temporal power dynamics across
the night.

## Usage

``` r
compute_temporal_bandpower(
  psg,
  hypno,
  channels = NULL,
  window_epochs = 10L,
  step_epochs = NULL,
  bands = list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 12), sigma = c(12, 16),
    beta = c(16, 30), gamma = c(30, 45)),
  win_sec = 4,
  noverlap = NULL,
  nfft = NULL
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- hypno:

  Integer vector of length `psg$n_epochs`. Sleep stage codes following
  the YASA convention: `0` = Wake, `1` = N1, `2` = N2, `3` = N3, `4` =
  REM, `-1` = Artefact. Used to assign a dominant stage to each window.

- channels:

  Character vector. Channel labels whose PSDs are averaged. If `NULL`
  (default), all non-bad EEG channels are used.

- window_epochs:

  Integer. Number of 30-s epochs per window. Default `10` (5 minutes).

- step_epochs:

  Integer. Window step size in epochs. Default `NULL`, which sets the
  step to `window_epochs %/% 2` (50 % overlap).

- bands:

  Named list of length-2 numeric vectors (Hz). Default matches
  [`compute_band_power()`](https://mrpheus.circadia-lab.uk/reference/compute_band_power.md):

      list(delta = c(0.5, 4), theta = c(4, 8), alpha = c(8, 12),
           sigma = c(12, 16), beta = c(16, 30), gamma = c(30, 45))

- win_sec:

  Numeric. Welch window length in seconds. Default `4`.

- noverlap:

  Integer or `NULL`. Welch overlap in samples. `NULL` (default) sets it
  to `win_sec / 2 * sample_rate`, i.e. 50 % of the Welch window —
  matching the notebook's `noverlap = int(2 * sfreq)`.

- nfft:

  Integer or `NULL`. FFT length. `NULL` (default) uses the next power of
  2 above `win_sec * sample_rate`.

## Value

A tibble with columns:

- time_hours:

  Numeric. Start time of the window in hours from recording onset.

- epoch_start:

  Integer. First epoch index in the window (1-based).

- epoch_end:

  Integer. Last epoch index in the window (1-based).

- dominant_stage:

  Integer. Most common stage code across the window (artefact epochs
  excluded). `NA` if all epochs are artefacts.

- band:

  Character. Band name.

- power:

  Numeric. Mean band power across channels (V^2/Hz).

- relative_power:

  Numeric. Band power as a fraction of total band-summed power at this
  time point.

## Details

Mirrors `calculate_temporal_bandpower()` from the companion Python
notebook: windows step by 50 % of their length (default), Welch overlap
defaults to half the window length in samples (proportional to sample
rate), and relative power is the share of each band in the total
band-summed power.

## Examples

``` r
if (FALSE) { # \dontrun{
rec   <- read_edf("SC4001E0-PSG.edf")
psg   <- prepare_psg(rec) |> preprocess_psg()
hypno <- stage_epochs(psg)$stage

tbp <- compute_temporal_bandpower(psg, hypno)

# Plot delta power across the night
library(ggplot2)
tbp |>
  dplyr::filter(band == "delta") |>
  ggplot(aes(time_hours, relative_power)) +
  geom_line() +
  labs(x = "Time (hours)", y = "Relative delta power")
} # }
```
