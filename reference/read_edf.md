# Read an EDF or EDF+ recording

Reads a European Data Format (EDF or EDF+) file and returns a structured
list containing signal data, channel metadata, and recording header.
Wraps
[`edfReader::readEdfHeader()`](https://rdrr.io/pkg/edfReader/man/readEdfHeader.html)
and
[`edfReader::readEdfSignals()`](https://rdrr.io/pkg/edfReader/man/readEdfSignals.html)
with consistent output formatting for the mrpheus pipeline.

## Usage

``` r
read_edf(path, channels = NULL, only_header = FALSE)
```

## Arguments

- path:

  Character. Path to an `.edf` or `.edf+` file.

- channels:

  Character vector or `NULL`. Channel labels to import. If `NULL`
  (default), all channels are imported.

- only_header:

  Logical. If `TRUE`, return only the header without reading signal
  data. Useful for quick channel inspection. Default `FALSE`.

## Value

A list of class `mrpheus_edf` with components:

- header:

  Data frame. Recording metadata (patient info, start time, number of
  signals, etc.).

- signals:

  Named list of numeric vectors, one per channel.

- channels:

  Data frame. Channel-level metadata: label, sample rate, physical
  min/max, digital min/max, transducer type, prefiltering.

- duration_s:

  Numeric. Total recording duration in seconds.

- path:

  Character. Resolved path to the source file.

## Examples

``` r
if (FALSE) { # \dontrun{
rec <- read_edf("data/psg_001.edf")
rec <- read_edf("data/psg_001.edf", channels = c("EEG Fpz-Cz", "EOG horizontal"))
rec$channels
} # }
```
