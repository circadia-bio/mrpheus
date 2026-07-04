# Compute Instantaneous Heart Rate from R-Peak Indices

Converts R-peak sample indices (from
[`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md))
into a sample-by-sample instantaneous heart rate trace by filling each
inter-beat interval with the constant HR derived from that interval's
length.

## Usage

``` r
compute_hr_signal(qrs)
```

## Arguments

- qrs:

  A `mrpheus_qrs` object from
  [`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md),
  or a named list with at least `qrs_i` (integer vector of R-peak
  indices) and `fs` (sampling frequency in Hz).

## Value

Numeric vector of length `max(qrs$qrs_i)`. Each sample contains the
instantaneous heart rate in beats per minute for the inter-beat interval
it falls within. Samples before the first detected R-peak are 0.

## See also

[`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md),
[`compute_hrv_sleep()`](https://mrpheus.circadia-lab.uk/reference/compute_hrv_sleep.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rec <- read_philips_physlog("sub-01_physlog.log")
qrs <- detect_qrs(rec$C[, "v1raw"], fs = rec$HDR$sfreq)
hr  <- compute_hr_signal(qrs)

plot(seq_along(hr) / rec$HDR$sfreq, hr, type = "l",
     xlab = "Time (s)", ylab = "HR (bpm)")
} # }
```
