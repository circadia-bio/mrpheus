# Orpheus mosaic palette — data object
# Generated from: Roman mosaic, Orpheus Charming the Animals
# 3rd century AD, Palermo Archaeological Museum
# Source: https://en.wikipedia.org/wiki/Orpheus
#
# Run this script to regenerate palette_orpheus and update R/sysdata.rda.
# usethis::use_data(palette_orpheus, internal = FALSE, overwrite = TRUE)

palette_orpheus <- c(
  sand       = "#CDB992",  # warm background tessera
  vermillion = "#B83E2C",  # Orpheus's robe; terracotta
  olive      = "#6A7840",  # tree and vegetation
  umber      = "#7C5432",  # animal fur and earth
  bistre     = "#3C2212",  # shadows and outlines
  ochre      = "#B07C3A",  # warm amber accent
  slate      = "#6C8284",  # birds; dusty teal-grey
  ivory      = "#EAD6AA"   # highlights and light tessera
)

usethis::use_data(palette_orpheus, overwrite = TRUE)
