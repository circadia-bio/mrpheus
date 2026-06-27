# Compute Oxygen Desaturation Index (ODI)

Calculates the ODI (number of \>= 3% SpO2 desaturations per hour of
sleep) from the SpO2 channel.

## Usage

``` r
compute_odi(psg, spo2_channel, tst_hours, threshold = 3)
```

## Arguments

- psg:

  An `mrpheus_psg` object from
  [`prepare_psg()`](https://mrpheus.circadia-lab.uk/reference/prepare_psg.md).

- spo2_channel:

  Character. SpO2 channel label.

- tst_hours:

  Numeric. Total sleep time in hours.

- threshold:

  Numeric. Desaturation threshold (percentage points). Default `3`.

## Value

Numeric scalar. ODI (desaturations per hour).
