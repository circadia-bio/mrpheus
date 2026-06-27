# Compute Apnea-Hypopnea Index (AHI)

Calculates the AHI from respiratory event data and total sleep time.

## Usage

``` r
compute_ahi(events, tst_hours)
```

## Arguments

- events:

  Output of
  [`detect_apneas()`](https://mrpheus.circadia-lab.uk/reference/detect_apneas.md).

- tst_hours:

  Numeric. Total sleep time in hours.

## Value

Numeric scalar. AHI (events per hour).
