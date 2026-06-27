#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang abort warn inform
#' @importFrom cli cli_alert_success cli_alert_warning cli_alert_info
#' @importFrom stats predict quantile sd var
#' @importFrom zoo rollapply
## usethis namespace: end
NULL

# Suppress R CMD check NOTE for dplyr NSE column names
utils::globalVariables(c("epoch", "channel", "total_power"))

#' Orpheus mosaic palette
#'
#' An 8-colour palette extracted from the Roman mosaic *Orpheus Charming the
#' Animals* (3rd century AD, Palermo Archaeological Museum). The mosaic depicts
#' Orpheus — the mythological figure whose name and story this package honours —
#' surrounded by animals, rendered in warm Mediterranean tesserae.
#'
#' The palette is intentionally earthy and muted, reflecting the natural pigments
#' of Roman mosaic work: sandy limestone backgrounds, terracotta robes, olive
#' vegetation, warm umber fauna, and the distinctive slate-teal of the birds.
#'
#' @format A named character vector of 8 hex colour codes:
#' \describe{
#'   \item{sand}{`#CDB992` — warm background tessera}
#'   \item{vermillion}{`#B83E2C` — Orpheus's robe; terracotta}
#'   \item{olive}{`#6A7840` — tree and vegetation}
#'   \item{umber}{`#7C5432` — animal fur and earth}
#'   \item{bistre}{`#3C2212` — shadows and outlines}
#'   \item{ochre}{`#B07C3A` — warm amber accent}
#'   \item{slate}{`#6C8284` — birds; dusty teal-grey}
#'   \item{ivory}{`#EAD6AA` — highlights and light tessera}
#' }
#'
#' @source Mosaic: *Orpheus Charming the Animals*, Roman, 3rd century AD.
#'   Palermo Archaeological Museum. Image via Wikimedia Commons,
#'   <https://en.wikipedia.org/wiki/Orpheus>.
#'
#' @examples
#' if (requireNamespace("scales", quietly = TRUE)) {
#'   scales::show_col(palette_orpheus)
#' }
#'
#' @export
palette_orpheus <- c(
  sand       = "#CDB992",
  vermillion = "#B83E2C",
  olive      = "#6A7840",
  umber      = "#7C5432",
  bistre     = "#3C2212",
  ochre      = "#B07C3A",
  slate      = "#6C8284",
  ivory      = "#EAD6AA"
)
