#' Read a Philips MRI Physiological Log File
#'
#' Reads a Philips PMU physiological log (`.log`) file recorded alongside an
#' MRI acquisition and returns a structured object containing the signal
#' matrix, event markers, and recording metadata. Supports all Philips PMU
#' hardware variants via the `system` preset.
#'
#' @param path     Character. Path to the Philips `.log` file.
#' @param system   Philips PMU hardware preset. One of:
#'   \describe{
#'     \item{`"wBTU"` (default)}{Wireless VCG system, 496 Hz. 10-channel
#'       layout: v1raw, v2raw, v1, v2, ppu, resp, gx, gy, gz, mark.
#'       Used with the Philips Achieva/Ingenia dStream wireless body
#'       telemetry unit.}
#'     \item{`"wired"`}{Older wired ECG system, 500 Hz. Fewer channels
#'       (no accelerometer); layout read from file header.}
#'     \item{`"custom"`}{Supply `sfreq` explicitly.}
#'   }
#' @param sfreq    Numeric or `NULL`. Sampling frequency in Hz. Overrides the
#'   `system` preset when supplied. Required when `system = "custom"`.
#' @param channels Character vector or `NULL`. Channel names to return.
#'   `NULL` (default) returns all channels. `"none"` returns no signal data
#'   (markers only).
#' @param skipprep Logical. If `TRUE`, skip samples recorded before the blank
#'   `#` line that marks the end of the preparation phase. Default `FALSE`.
#'
#' @return A list of class `mrpheus_physlog` with components:
#' \describe{
#'   \item{`C`}{Integer matrix (samples × channels). Column order matches the
#'     file header.}
#'   \item{`M`}{Two-column integer matrix: `marker` (bit-encoded value) and
#'     `index` (1-based sample number) for every non-zero marker.}
#'   \item{`I`}{Named list of 1-based sample index vectors for each event
#'     type: VcgOnset, PpuOnset, TriggerResp, Measurement, ScannerStart,
#'     ScannerStop, TriggerExt, Calibration, RefTriggerVcg.}
#'   \item{`HDR`}{List of header metadata: ID, DATETIME, STATS,
#'     DockableTable, COLUMN_NAMES, system, sfreq.}
#'   \item{`path`}{Character. Resolved path to the source file.}
#' }
#'
#' @details
#' Bit-encoded marker values:
#' ```
#' 0x0001  ECG trigger       0x0002  PPU trigger      0x0004  respiration
#' 0x0008  slice onset       0x0010  scanner start    0x0020  scanner stop
#' 0x0040  external trigger  0x0080  calibration      0x8000  ref ECG trigger
#' ```
#'
#' Use the scanner-stop marker (more reliable than start) to time-lock data:
#' ```r
#' stopidx  <- tail(data$I$ScannerStop, 1) - enddelay
#' startidx <- stopidx - round(n_vols * TR * data$HDR$sfreq)
#' epoch    <- data$C[startidx:stopidx, ]
#' ```
#'
#' @seealso [mrpheus::detect_qrs()], [mrpheus::compute_hr_signal()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Philips Achieva/Ingenia with wireless VCG (default)
#' rec <- read_philips_physlog("sub-01_ses-01_physlog.log")
#' rec
#'
#' # Older wired scanner
#' rec <- read_philips_physlog("sub-01_physlog.log", system = "wired")
#'
#' # Custom sampling rate
#' rec <- read_philips_physlog("sub-01_physlog.log", system = "custom", sfreq = 400)
#'
#' # ECG channels only
#' rec <- read_philips_physlog("sub-01_physlog.log", channels = c("v1raw", "v2raw"))
#' }
read_philips_physlog <- function(
    path,
    system   = c("wBTU", "wired", "custom"),
    sfreq    = NULL,
    channels = NULL,
    skipprep = FALSE
) {

  system <- match.arg(system)

  if (!file.exists(path))
    cli::cli_abort("File not found: {.path {path}}")

  # ---- Resolve sampling frequency ----------------------------------------
  sfreq_presets <- c(wBTU = 496, wired = 500)

  if (!is.null(sfreq)) {
    if (!is.numeric(sfreq) || length(sfreq) != 1L || sfreq <= 0)
      cli::cli_abort("`sfreq` must be a single positive number.")
    sfreq_resolved <- sfreq
  } else if (system == "custom") {
    cli::cli_abort('`sfreq` must be supplied when `system = "custom"`.')
  } else {
    sfreq_resolved <- sfreq_presets[[system]]
  }

  lines <- readLines(path, warn = FALSE)

  # ---- Header regex patterns -----------------------------------------------
  p1 <- "##\\s*(.*),\\s*Release\\s*(\\w+)\\s*\\(SWID (\\d+)\\)"
  p2 <- "##.*?(\\d{2})-(\\d{2})-(\\d{4})\\s+(\\d{2}):(\\d{2}):(\\d{2})"
  p3 <- "##\\s+([-0-9 ]+)$"
  p4 <- "##\\s*Dockable table\\s*=\\s*(\\w+)"

  hdr        <- list()
  hdr_lines  <- character(0)
  col_names  <- NULL
  data_start <- NULL

  for (i in seq_along(lines)) {
    line <- lines[[i]]
    if (startsWith(line, "##")) {
      hdr_lines <- c(hdr_lines, line)
    } else if (startsWith(line, "#")) {
      col_names  <- regmatches(line, gregexpr("\\w+", line))[[1L]]
      data_start <- i + 1L
      break
    }
  }

  if (is.null(data_start))
    cli::cli_abort("Could not find column-names line in: {.path {path}}")

  # ---- Parse header fields -------------------------------------------------
  if (length(hdr_lines) >= 1L) {
    m <- regmatches(hdr_lines[1L], regexec(p1, hdr_lines[1L]))[[1L]]
    if (length(m) > 1L)
      hdr$ID <- list(Site = m[2L], Release = m[3L], SWID = as.numeric(m[4L]))
  }
  if (length(hdr_lines) >= 2L) {
    m <- regmatches(hdr_lines[2L], regexec(p2, hdr_lines[2L]))[[1L]]
    if (length(m) > 1L)
      hdr$DATETIME <- list(
        year  = as.integer(m[4L]), month = as.integer(m[3L]),
        day   = as.integer(m[2L]), hour  = as.integer(m[5L]),
        min   = as.integer(m[6L]), sec   = as.integer(m[7L])
      )
  }
  if (length(hdr_lines) >= 3L) {
    m <- regmatches(hdr_lines[3L], regexec(p3, hdr_lines[3L]))[[1L]]
    if (length(m) > 1L)
      hdr$STATS <- as.integer(strsplit(trimws(m[2L]), "\\s+")[[1L]])
  }
  if (length(hdr_lines) >= 4L) {
    m <- regmatches(hdr_lines[4L], regexec(p4, hdr_lines[4L]))[[1L]]
    if (length(m) > 1L)
      hdr$DockableTable <- toupper(m[2L]) == "TRUE"
  }
  hdr$COLUMN_NAMES <- col_names
  hdr$system       <- system
  hdr$sfreq        <- sfreq_resolved

  # ---- Skip preparation phase if requested ---------------------------------
  if (skipprep) {
    found <- FALSE
    for (j in seq(data_start, length(lines))) {
      if (trimws(lines[[j]]) == "#") {
        data_start <- j + 1L
        found <- TRUE
        break
      }
    }
    if (!found)
      cli::cli_warn("Preparation-phase marker not found; reading all samples.")
  }

  # ---- Collect data lines --------------------------------------------------
  data_lines <- lines[seq(data_start, length(lines))]
  data_lines <- data_lines[!grepl("^\\s*#", data_lines)]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0L]

  if (length(data_lines) == 0L)
    cli::cli_abort("No data lines found in: {.path {path}}")

  n_cols <- length(col_names)

  con <- textConnection(data_lines)
  on.exit(try(close(con), silent = TRUE))
  raw <- read.table(
    con,
    header     = FALSE,
    colClasses = c(rep("integer", n_cols - 1L), "character"),
    col.names  = col_names
  )
  close(con)
  on.exit()

  # ---- Parse hex marker column -> integer ----------------------------------
  markers <- strtoi(raw[[n_cols]], 16L)

  # ---- Event index vectors (1-based sample indices) ------------------------
  I <- list(
    VcgOnset      = which(bitwAnd(markers,     1L) > 0L),
    PpuOnset      = which(bitwAnd(markers,     2L) > 0L),
    TriggerResp   = which(bitwAnd(markers,     4L) > 0L),
    Measurement   = which(bitwAnd(markers,     8L) > 0L),
    ScannerStart  = which(bitwAnd(markers,    16L) > 0L),
    ScannerStop   = which(bitwAnd(markers,    32L) > 0L),
    TriggerExt    = which(bitwAnd(markers,    64L) > 0L),
    Calibration   = which(bitwAnd(markers,   128L) > 0L),
    RefTriggerVcg = which(bitwAnd(markers, 32768L) > 0L)
  )

  raw[[n_cols]] <- markers

  # ---- Channel selection ---------------------------------------------------
  if (is.null(channels) ||
      identical(channels, "all") ||
      length(channels) == 0L) {
    C <- as.matrix(raw)
  } else if (identical(channels, "none")) {
    C <- matrix(integer(0L), nrow = nrow(raw), ncol = 0L)
  } else {
    valid   <- channels[channels %in% col_names]
    col_idx <- match(valid, col_names)
    C <- as.matrix(raw[, col_idx, drop = FALSE])
  }
  storage.mode(C) <- "integer"

  # ---- Marker table --------------------------------------------------------
  idx <- which(markers > 0L)
  M   <- cbind(marker = markers[idx], index = idx)

  structure(
    list(
      C    = C,
      M    = M,
      I    = I,
      HDR  = hdr,
      path = normalizePath(path)
    ),
    class = "mrpheus_physlog"
  )
}

#' @export
print.mrpheus_physlog <- function(x, ...) {
  n_samples  <- nrow(x$C)
  duration_s <- n_samples / x$HDR$sfreq
  n_markers  <- sum(sapply(x$I, length))

  cli::cli_h1("mrpheus Philips PMU recording")
  cli::cli_inform(c(
    "i" = "Path:     {.path {x$path}}",
    "i" = "System:   {x$HDR$system} ({x$HDR$sfreq} Hz)",
    "i" = "Channels: {paste(x$HDR$COLUMN_NAMES, collapse = ', ')}",
    "i" = "Samples:  {n_samples} ({round(duration_s / 60, 1)} min)",
    "i" = "Markers:  {n_markers} events"
  ))
  invisible(x)
}
