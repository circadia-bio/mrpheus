# Detect respiratory events (apneas and hypopneas)

Detects apneas and hypopneas in respiratory airflow and effort signals,
returning event-level metadata and summary indices. Detection follows
AASM 2012/2017 criteria: \>= 90% signal reduction for \>= 10 s (apnea),
or \>= 30% reduction with \>= 3% SpO2 desaturation or arousal for \>= 10
s (hypopnea).

## Usage

``` r
detect_apneas(
  psg,
  airflow_channel = NULL,
  spo2_channel = NULL,
  min_duration_s = 10,
  apnea_threshold = 0.9,
  hypopnea_threshold = 0.3,
  desaturation_threshold = 3
)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- airflow_channel:

  Character. Airflow channel label (e.g. nasal pressure or thermistor).
  If `NULL` (default), the first RESP channel is used.

- spo2_channel:

  Character or `NULL`. SpO2 channel for hypopnea scoring.

- min_duration_s:

  Numeric. Minimum event duration in seconds. Default `10`.

- apnea_threshold:

  Numeric. Fractional reduction required for apnea classification.
  Default `0.90`.

- hypopnea_threshold:

  Numeric. Fractional reduction required for hypopnea classification.
  Default `0.30`.

- desaturation_threshold:

  Numeric. SpO2 drop (percentage points) required to confirm hypopnea.
  Default `3`.

## Value

A list with:

- events:

  Tibble. One row per event: `epoch`, `start_s`, `end_s`, `duration_s`,
  `type` (`"apnea"` / `"hypopnea"`), `desaturation`.

- summary:

  Tibble. `n_apneas`, `n_hypopneas`, `n_events`.

## See also

[`compute_ahi()`](https://mrpheus.circadia-lab.uk/reference/compute_ahi.md),
[`compute_odi()`](https://mrpheus.circadia-lab.uk/reference/compute_odi.md)
