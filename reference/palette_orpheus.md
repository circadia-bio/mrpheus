# Orpheus mosaic palette

An 8-colour palette extracted from the Roman mosaic *Orpheus Charming
the Animals* (3rd century AD, Palermo Archaeological Museum). The mosaic
depicts Orpheus — the mythological figure whose name and story this
package honours — surrounded by animals, rendered in warm Mediterranean
tesserae.

## Usage

``` r
palette_orpheus
```

## Format

A named character vector of 8 hex colour codes:

- sand:

  `#CDB992` — warm background tessera

- vermillion:

  `#B83E2C` — Orpheus's robe; terracotta

- olive:

  `#6A7840` — tree and vegetation

- umber:

  `#7C5432` — animal fur and earth

- bistre:

  `#3C2212` — shadows and outlines

- ochre:

  `#B07C3A` — warm amber accent

- slate:

  `#6C8284` — birds; dusty teal-grey

- ivory:

  `#EAD6AA` — highlights and light tessera

## Source

Mosaic: *Orpheus Charming the Animals*, Roman, 3rd century AD. Palermo
Archaeological Museum. Image via Wikimedia Commons,
<https://en.wikipedia.org/wiki/Orpheus>.

## Details

The palette is intentionally earthy and muted, reflecting the natural
pigments of Roman mosaic work: sandy limestone backgrounds, terracotta
robes, olive vegetation, warm umber fauna, and the distinctive
slate-teal of the birds.

## Examples

``` r
if (requireNamespace("scales", quietly = TRUE)) {
  scales::show_col(palette_orpheus)
}

```
