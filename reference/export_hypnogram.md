# Export a staged hypnogram to hypnor

Converts the staging tibble produced by
[`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
into a `hypnor_hypnogram` object ready for architecture metric
computation and visualisation in the `hypnor` package. This is the
primary handoff between `mrpheus` and the rest of the Circadia Lab
ecosystem.

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

An object of class `hypnor_hypnogram` (defined in `hypnor`). If `hypnor`
is not installed, returns the staging tibble invisibly with a message.

## See also

[`stage_epochs()`](https://mrpheus.circadia-lab.uk/reference/stage_epochs.md)
