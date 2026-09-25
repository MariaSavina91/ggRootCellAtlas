test_that("palette has one named colour per non-NA group across all maps", {
  maps <- load_maps()
  pal <- do.call(generate_common_palette, c(list("Zones"), unname(maps)))
  expect_setequal(names(pal), annotation_levels("Zones"))
  expect_false(anyNA(pal))
  expect_false(anyDuplicated(pal) > 0)
})

test_that("palette accepts a custom palette function", {
  maps <- load_maps()
  pal <- generate_common_palette("TissueTypes", maps[[1]],
                                 color_palette = function(n) rep("red", n))
  expect_true(all(pal == "red"))
})

test_that("okabe_ito_pal gives 11 distinct colours, then falls back to hues", {
  pal <- okabe_ito_pal()
  expect_length(unique(pal(11)), 11)
  expect_identical(pal(3), pal(11)[1:3])
  expect_identical(pal(20), scales::hue_pal()(20))
})

test_that("palette works on tibbles", {
  m <- dplyr::tibble(load_maps()[[2]])
  pal <- generate_common_palette("Zones", m)
  expect_setequal(names(pal), unique(stats::na.omit(m$Zones)))
})
