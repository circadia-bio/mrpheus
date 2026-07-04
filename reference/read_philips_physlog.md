# Read a Philips MRI Physiological Log File

Reads a Philips PMU physiological log (`.log`) file recorded alongside
an MRI acquisition and returns a structured object containing the signal
matrix, event markers, and recording metadata. Supports all Philips PMU
hardware variants via the `system` preset.

## Usage

``` r
read_philips_physlog(
  path,
  system = c("wBTU", "wired", "custom"),
  sfreq = NULL,
  channels = NULL,
  skipprep = FALSE
)
```

## Arguments

- path:

  Character. Path to the Philips `.log` file.

- system:

  Philips PMU hardware preset. One of:

  `"wBTU"` (default)

  :   Wireless VCG system, 496 Hz. 10-channel layout: v1raw, v2raw, v1,
      v2, ppu, resp, gx, gy, gz, mark. Used with the Philips
      Achieva/Ingenia dStream wireless body telemetry unit.

  `"wired"`

  :   Older wired ECG system, 500 Hz. Fewer channels (no accelerometer);
      layout read from file header.

  `"custom"`

  :   Supply `sfreq` explicitly.

- sfreq:

  Numeric or `NULL`. Sampling frequency in Hz. Overrides the `system`
  preset when supplied. Required when `system = "custom"`.

- channels:

  Character vector or `NULL`. Channel names to return. `NULL` (default)
  returns all channels. `"none"` returns no signal data (markers only).

- skipprep:

  Logical. If `TRUE`, skip samples recorded before the blank `#` line
  that marks the end of the preparation phase. Default `FALSE`.

## Value

A list of class `mrpheus_physlog` with components:

- `C`:

  Integer matrix (samples × channels). Column order matches the file
  header.

- `M`:

  Two-column integer matrix: `marker` (bit-encoded value) and `index`
  (1-based sample number) for every non-zero marker.

- `I`:

  Named list of 1-based sample index vectors for each event type:
  VcgOnset, PpuOnset, TriggerResp, Measurement, ScannerStart,
  ScannerStop, TriggerExt, Calibration, RefTriggerVcg.

- `HDR`:

  List of header metadata: ID, DATETIME, STATS, DockableTable,
  COLUMN_NAMES, system, sfreq.

- `path`:

  Character. Resolved path to the source file.

## Details

Bit-encoded marker values:

    0x0001  ECG trigger       0x0002  PPU trigger      0x0004  respiration
    0x0008  slice onset       0x0010  scanner start    0x0020  scanner stop
    0x0040  external trigger  0x0080  calibration      0x8000  ref ECG trigger

Use the scanner-stop marker (more reliable than start) to time-lock
data:

    stopidx  <- tail(data$I$ScannerStop, 1) - enddelay
    startidx <- stopidx - round(n_vols * TR * data$HDR$sfreq)
    epoch    <- data$C[startidx:stopidx, ]

## See also

[`detect_qrs()`](https://mrpheus.circadia-lab.uk/reference/detect_qrs.md),
[`compute_hr_signal()`](https://mrpheus.circadia-lab.uk/reference/compute_hr_signal.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Philips Achieva/Ingenia with wireless VCG (default)
rec <- read_philips_physlog("sub-01_ses-01_physlog.log")
rec

# Older wired scanner
rec <- read_philips_physlog("sub-01_physlog.log", system = "wired")

# Custom sampling rate
rec <- read_philips_physlog("sub-01_physlog.log", system = "custom", sfreq = 400)

# ECG channels only
rec <- read_philips_physlog("sub-01_physlog.log", channels = c("v1raw", "v2raw"))
} # }
```
