# Remove DC offset from a signal

Subtracts the signal mean (per-channel DC offset removal). Should be
applied before filtering to prevent filter transients from a non-zero
baseline.

## Usage

``` r
remove_dc(signal)
```

## Arguments

- signal:

  Numeric vector.

## Value

Numeric vector with mean removed, same length as `signal`.

## Examples

``` r
x <- c(100.1, 100.5, 99.8, 100.3)
remove_dc(x)
#> [1] -0.075  0.325 -0.375  0.125
```
