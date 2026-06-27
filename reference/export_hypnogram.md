# Export a staged hypnogram for use with hypnor

Prepares the staging tibble produced by
[`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
for downstream use with the `hypnor` package. Attaches recording
metadata as attributes and returns a tibble of class
`mrpheus_hypnogram`, which `hypnor::new_hypnogram()` accepts directly
once `hypnor` is installed.

## Usage

``` r
export_hypnogram(
  staging,
  epoch_s = 30,
  start_time = NULL,
  participant_id = NULL
)
```

## Arguments

- staging:

  A tibble from
  [`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
  with columns `epoch`, `stage`, and optional probability columns.

- epoch_s:

  Numeric. Epoch duration in seconds. Must match the `epoch_s` used in
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).
  Default `30`.

- start_time:

  POSIXct or `NULL`. Recording start time. Used to compute clock-time
  axes in `hypnor` plots. If `NULL`, epochs are indexed from 0.

- participant_id:

  Character or `NULL`. Optional identifier passed through to `hypnor`
  and `syncR`.

## Value

A tibble of class `mrpheus_hypnogram` with columns `epoch`, `stage`, and
any probability columns from the staging model. Metadata (`epoch_s`,
`start_time`, `participant_id`, `source`, `resolution`) are attached as
attributes and forwarded to `hypnor::new_hypnogram()` when `hypnor` is
available.

## See also

[`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
