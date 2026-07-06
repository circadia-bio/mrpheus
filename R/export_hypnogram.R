# R/export_hypnogram.R
#
# Export staged hypnogram to mrpheus_hypnogram S3 object with pre-computed
# AASM sleep architecture metrics. Ported directly from yasa.sleep_statistics()
# (src/yasa/sleepstats.py) — metric definitions and computation match YASA
# exactly so values can be compared directly between the two pipelines.
#
# Reference: Iber C (2007). The AASM Manual for the Scoring of Sleep and
# Associated Events. American Academy of Sleep Medicine.

# ── Internal: sleep architecture computation ──────────────────────────────────
#
# Direct port of yasa.sleep_statistics(hypno, sf_hyp).
# All outputs are in minutes (or % for SE, SME, and stage percentages).
#
# Stage coding (YASA convention):
#   -2 = Unscored, -1 = Artefact, 0 = Wake, 1 = N1, 2 = N2, 3 = N3, 4 = REM
#
# NOTE on latency definition (matches YASA, not AASM):
#   Lat_REM is from the START OF RECORDING, not from sleep onset.
#   AASM REM latency = Lat_REM - SOL (also returned as Lat_REM_AASM).

.sleep_architecture <- function(hypno, epoch_s = 30) {
  ep_to_min <- epoch_s / 60   # one epoch in minutes

  # TIB = total recording duration
  tib <- length(hypno) * ep_to_min

  # First and last sleep epoch (any stage > 0; artefacts < 0 are excluded)
  sleep_idx <- which(hypno > 0L)

  if (length(sleep_idx) == 0L) {
    # No sleep detected at all
    na <- NA_real_
    return(list(
      TIB = tib, SPT = 0, WASO = na, TST = 0,
      N1 = 0, N2 = 0, N3 = 0, REM = 0, NREM = 0,
      SOL = na,
      Lat_N1 = na, Lat_N2 = na, Lat_N3 = na, Lat_REM = na,
      Lat_REM_AASM = na,
      pct_N1 = na, pct_N2 = na, pct_N3 = na,
      pct_REM = na, pct_NREM = na,
      SE = 0, SME = na
    ))
  }

  first_sleep <- sleep_idx[1L]
  last_sleep  <- sleep_idx[length(sleep_idx)]

  # SPT: from first to last sleep epoch (inclusive)
  hypno_s <- hypno[first_sleep:last_sleep]
  spt     <- length(hypno_s) * ep_to_min
  waso    <- sum(hypno_s == 0L) * ep_to_min
  # TST: sleep epochs only (excludes Art/Unscored within SPT) — matches YASA >= 0.5.0
  tst     <- sum(hypno_s > 0L) * ep_to_min

  # Stage durations (full recording, not just SPT)
  n1   <- sum(hypno == 1L) * ep_to_min
  n2   <- sum(hypno == 2L) * ep_to_min
  n3   <- sum(hypno == 3L) * ep_to_min
  rem  <- sum(hypno == 4L) * ep_to_min
  nrem <- n1 + n2 + n3

  # Latencies from START OF RECORDING (0-based epoch index * ep_to_min)
  # Matches YASA: stats["SOL"] = first_sleep (0-based); stats["Lat_N1"] = min(where(hypno==1))
  sol     <- (first_sleep - 1L) * ep_to_min
  lat_n1  <- if (any(hypno == 1L)) (min(which(hypno == 1L)) - 1L) * ep_to_min else NA_real_
  lat_n2  <- if (any(hypno == 2L)) (min(which(hypno == 2L)) - 1L) * ep_to_min else NA_real_
  lat_n3  <- if (any(hypno == 3L)) (min(which(hypno == 3L)) - 1L) * ep_to_min else NA_real_
  lat_rem <- if (any(hypno == 4L)) (min(which(hypno == 4L)) - 1L) * ep_to_min else NA_real_

  # AASM REM latency = time from sleep onset to first REM epoch
  lat_rem_aasm <- if (!is.na(lat_rem)) lat_rem - sol else NA_real_

  # Stage percentages (% of TST)
  pct <- function(x) if (tst > 0) 100 * x / tst else NA_real_

  # Sleep efficiency
  se  <- 100 * tst / tib
  sme <- if (spt > 0) 100 * tst / spt else NA_real_

  list(
    TIB          = tib,
    SPT          = spt,
    WASO         = waso,
    TST          = tst,
    N1           = n1,
    N2           = n2,
    N3           = n3,
    REM          = rem,
    NREM         = nrem,
    SOL          = sol,
    Lat_N1       = lat_n1,
    Lat_N2       = lat_n2,
    Lat_N3       = lat_n3,
    Lat_REM      = lat_rem,       # from start of recording (YASA convention)
    Lat_REM_AASM = lat_rem_aasm,  # from sleep onset (AASM convention)
    pct_N1       = pct(n1),
    pct_N2       = pct(n2),
    pct_N3       = pct(n3),
    pct_REM      = pct(rem),
    pct_NREM     = pct(nrem),
    SE           = se,
    SME          = sme
  )
}

# ── Public function ───────────────────────────────────────────────────────────

#' Export a staged hypnogram for use with hypnor
#'
#' Wraps the staging tibble from [mrpheus::stage_epochs()] in the
#' `mrpheus_hypnogram` S3 class, computes AASM sleep architecture metrics, and
#' attaches them as the `sleep_architecture` attribute. This is the intended
#' input for `hypnor::new_hypnogram()` once `hypnor` is available.
#'
#' Sleep architecture metrics are a direct R port of `yasa.sleep_statistics()`
#' and produce identical values. All durations are in **minutes**; efficiency
#' and stage percentages are in **percent**.
#'
#' @section Sleep architecture metrics:
#' \describe{
#'   \item{TIB}{Time in Bed: total recording duration.}
#'   \item{SPT}{Sleep Period Time: first to last sleep epoch.}
#'   \item{WASO}{Wake After Sleep Onset: wake time within SPT.}
#'   \item{TST}{Total Sleep Time: N1 + N2 + N3 + REM within SPT
#'     (artefact epochs excluded; matches YASA >= 0.5.0).}
#'   \item{N1, N2, N3, REM, NREM}{Stage durations.}
#'   \item{SOL}{Sleep Onset Latency: time to first sleep epoch.}
#'   \item{Lat_N1, Lat_N2, Lat_N3, Lat_REM}{Latency to each stage from
#'     recording start (YASA convention).}
#'   \item{Lat_REM_AASM}{REM latency from sleep onset (AASM convention =
#'     \code{Lat_REM - SOL}).}
#'   \item{pct_N1, pct_N2, pct_N3, pct_REM, pct_NREM}{Stage durations as
#'     percentage of TST.}
#'   \item{SE}{Sleep Efficiency: TST / TIB × 100.}
#'   \item{SME}{Sleep Maintenance Efficiency: TST / SPT × 100.}
#' }
#'
#' @param staging A tibble from [mrpheus::stage_epochs()] with at minimum a
#'   `stage` column of integer stage codes (YASA convention: 0 = Wake,
#'   1 = N1, 2 = N2, 3 = N3, 4 = REM, -1 = Artefact).
#' @param epoch_s Numeric. Epoch duration in seconds. Default `30`.
#' @param start_time POSIXct or `NULL`. Recording start time. Used for
#'   clock-time axes in `hypnor` visualisations. Default `NULL`.
#' @param participant_id Character or `NULL`. Optional participant identifier
#'   forwarded to `hypnor` and `syncR`. Default `NULL`.
#'
#' @return A tibble of class `mrpheus_hypnogram` with the original staging
#'   columns plus the following attributes:
#' \describe{
#'   \item{epoch_s}{Epoch duration in seconds.}
#'   \item{start_time}{Recording start time (POSIXct or `NULL`).}
#'   \item{participant_id}{Participant identifier (character or `NULL`).}
#'   \item{source}{Always `"mrpheus"`.}
#'   \item{resolution}{Always `"AASM"`.}
#'   \item{sleep_architecture}{Named list of AASM sleep statistics (minutes /
#'     percent). See **Sleep architecture metrics** section.}
#' }
#'
#' @seealso [mrpheus::stage_epochs()]
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rec    <- read_edf("psg.edf")
#' psg    <- prepare_psg(rec) |> preprocess_psg()
#' stages <- stage_epochs(psg)
#' hyp    <- export_hypnogram(stages, participant_id = "sub-001")
#'
#' # Access sleep architecture
#' hyp |> attr("sleep_architecture")
#' }
export_hypnogram <- function(staging,
                             epoch_s        = 30,
                             start_time     = NULL,
                             participant_id = NULL) {
  if (!is.data.frame(staging) || !"stage" %in% names(staging)) {
    cli::cli_abort("`staging` must be a tibble from {.fn stage_epochs}.")
  }

  arch <- .sleep_architecture(as.integer(staging$stage), epoch_s = epoch_s)

  out <- staging
  attr(out, "epoch_s")            <- epoch_s
  attr(out, "start_time")         <- start_time
  attr(out, "participant_id")     <- participant_id
  attr(out, "source")             <- "mrpheus"
  attr(out, "resolution")         <- "AASM"
  attr(out, "sleep_architecture") <- arch

  class(out) <- c("mrpheus_hypnogram", class(out))

  cli::cli_alert_success(
    "Hypnogram ready: {nrow(out)} epochs -- \\
    TST {round(arch$TST, 1)} min, SE {round(arch$SE, 1)} %. \\
    Pass to {.code hypnor::new_hypnogram()} once {.pkg hypnor} is available."
  )
  out
}

# ── S3 methods ────────────────────────────────────────────────────────────────

#' @export
print.mrpheus_hypnogram <- function(x, ...) {
  arch <- attr(x, "sleep_architecture")
  cli::cli_h1("mrpheus hypnogram")
  cli::cli_inform(c(
    "i" = "Epochs:      {nrow(x)} x {attr(x, 'epoch_s')} s",
    "i" = "Participant: {attr(x, 'participant_id') %||% 'unset'}",
    "i" = "Source:      {attr(x, 'source')} / {attr(x, 'resolution')}"
  ))

  if (!is.null(arch)) {
    cli::cli_h2("Sleep architecture")
    cli::cli_inform(c(
      " " = "TIB:  {round(arch$TIB,  1)} min  |  SPT:  {round(arch$SPT,  1)} min",
      " " = "TST:  {round(arch$TST,  1)} min  |  WASO: {round(arch$WASO, 1)} min",
      " " = "SE:   {round(arch$SE,   1)} %    |  SME:  {round(arch$SME,  1)} %",
      " " = "SOL:  {round(arch$SOL,  1)} min  |  Lat_REM (AASM): {round(arch$Lat_REM_AASM, 1)} min",
      " " = "N1:   {round(arch$pct_N1,  1)} %  |  N2:  {round(arch$pct_N2,  1)} %  |  N3: {round(arch$pct_N3, 1)} %  |  REM: {round(arch$pct_REM, 1)} %"
    ))
  }
  NextMethod()
  invisible(x)
}
