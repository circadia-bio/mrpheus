test_that("palette_orpheus has 8 named colours", {
  expect_length(palette_orpheus, 8)
  expect_named(palette_orpheus)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", palette_orpheus)))
})
