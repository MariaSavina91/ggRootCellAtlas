# Fraction of annotated cells in panel i that received no expression value
unmatched <- function(p, i, column = "Atlas") {
  d <- p[[i]]$data
  mean(is.na(d$Gene[!is.na(d[[column]])]))
}

test_that("accepts a plain genes x groups matrix", {
  p <- ggRootCellAtlas_expression(fake_avg_exp(), "GENE1")
  expect_s3_class(p, "patchwork")
  for (i in 1:7) expect_equal(unmatched(p, i), 0)
})

test_that("accepts Seurat v4 AverageExpression() output (named list)", {
  p <- ggRootCellAtlas_expression(list(RNA = fake_avg_exp()), "GENE1")
  for (i in 1:7) expect_equal(unmatched(p, i), 0)
})

test_that("accepts as.data.frame() of the Seurat list ('RNA.' prefix)", {
  df <- as.data.frame(list(RNA = fake_avg_exp()))
  expect_true(all(startsWith(colnames(df), "RNA.")))
  p <- ggRootCellAtlas_expression(df, "GENE1")
  for (i in 1:7) expect_equal(unmatched(p, i), 0)
})

test_that("uses the first assay of a multi-assay list, with a message", {
  input <- list(RNA = fake_avg_exp(), SCT = fake_avg_exp() * 0)
  expect_message(p <- ggRootCellAtlas_expression(input, "GENE1"), "RNA")
  expect_gt(max(p[[1]]$data$Gene, na.rm = TRUE), 0)
})

test_that("accepts Seurat v5 group names (underscore replaced by dash)", {
  p <- ggRootCellAtlas_expression(seurat_v5_names(fake_avg_exp()), "GENE1")
  for (i in 1:7) expect_equal(unmatched(p, i), 0)
})

test_that("accepts sparse matrices (Seurat v5 may return dgCMatrix)", {
  skip_if_not_installed("Matrix")
  m <- Matrix::Matrix(fake_avg_exp(), sparse = TRUE)
  p <- ggRootCellAtlas_expression(list(RNA = m), "GENE1")
  for (i in 1:7) expect_equal(unmatched(p, i), 0)
})

test_that("every annotation column can be used for matching", {
  for (col in annotation_columns) {
    m <- fake_avg_exp(col)
    inputs <- list(matrix = m, v4 = list(RNA = m), v5 = seurat_v5_names(m))
    for (nm in names(inputs)) {
      p <- tryCatch(
        ggRootCellAtlas_expression(inputs[[nm]], "GENE1", Annotation = col),
        error = function(e) e
      )
      expect_false(inherits(p, "error"), info = paste(col, nm))
      if (inherits(p, "error")) next
      for (i in 1:7) expect_equal(unmatched(p, i, col), 0, info = paste(col, nm, i))
    }
  }
})

test_that("plotted values equal the input values", {
  m <- fake_avg_exp()
  p <- ggRootCellAtlas_expression(m, "GENE2")
  d <- p[[2]]$data
  d <- d[!is.na(d$Atlas), ]
  expect_equal(unname(as.numeric(d$Gene)), unname(m["GENE2", d$Atlas]))
})

test_that("all panels share the same colour limits", {
  p <- ggRootCellAtlas_expression(fake_avg_exp(), "GENE1")
  lims <- lapply(1:7, function(i) {
    ggplot2::ggplot_build(p[[i]])$plot$scales$get_scales("fill")$get_limits()
  })
  for (l in lims) expect_equal(l, lims[[1]])
})

test_that("values outside c1..c2 are clamped, not drawn as missing (issue 002)", {
  p <- ggRootCellAtlas_expression(fake_avg_exp(), "GENE1", c1 = 1, c2 = 2,
                                  na.colour = "black")
  for (i in 1:7) {
    fill <- panel_data(p, i)$fill
    has_value <- !is.na(p[[i]]$data$Gene)
    expect_false(any(fill[has_value] == "black"), info = i)
  }
})

test_that("cells without a value use na.colour, dark grey by default", {
  m <- fake_avg_exp()
  m <- m[, !grepl("_d$", colnames(m))]  # no values for the d section
  d_panel <- 7
  missing <- is.na(ggRootCellAtlas_expression(m, "GENE1")[[d_panel]]$data$Gene)
  expect_true(any(missing))

  fill <- panel_data(ggRootCellAtlas_expression(m, "GENE1"), d_panel)$fill
  expect_true(all(hex(fill[missing]) == hex("grey30")))
  fill <- panel_data(ggRootCellAtlas_expression(m, "GENE1", na.colour = "pink"), d_panel)$fill
  expect_true(all(hex(fill[missing]) == hex("pink")))
})

# Fill scale of panel i, to map values to colours directly
fill_scale <- function(p, i = 1) {
  ggplot2::ggplot_build(p[[i]])$plot$scales$get_scales("fill")
}

test_that("default colours run from snow2 to red4", {
  sc <- fill_scale(ggRootCellAtlas_expression(fake_avg_exp(), "GENE1", c1 = 1, c2 = 2))
  expect_equal(hex(sc$map(1)), hex("snow2"))
  expect_equal(hex(sc$map(2)), hex("red4"))
})

test_that("custom colours are used for the gradient", {
  sc <- fill_scale(ggRootCellAtlas_expression(fake_avg_exp(), "GENE1", c1 = 1, c2 = 2,
                                              colours = c("white", "darkgreen")))
  expect_equal(hex(sc$map(1)), hex("white"))
  expect_equal(hex(sc$map(2)), hex("darkgreen"))
})

test_that("midpoint centres a diverging scale with symmetric default limits", {
  m <- fake_avg_exp() - 1  # values from -1 to 4
  sc <- fill_scale(ggRootCellAtlas_expression(m, "GENE1", midpoint = 0))
  spread <- max(abs(m["GENE1", ]))
  expect_equal(sc$get_limits(), c(-spread, spread))
  # the middle colour (white by default) sits exactly at the midpoint
  expect_equal(hex(sc$map(0)), hex("white"))
  expect_equal(hex(sc$map(-spread)), "#2166AC")
  expect_equal(hex(sc$map(spread)), "#B2182B")
})

test_that("midpoint works with asymmetric limits and custom colours", {
  m <- fake_avg_exp() - 1
  sc <- fill_scale(ggRootCellAtlas_expression(
    m, "GENE1", c1 = -1, c2 = 4, midpoint = 0,
    colours = c("navy", "grey90", "orange", "darkred")
  ))
  expect_equal(hex(sc$map(-1)), hex("navy"))
  expect_equal(hex(sc$map(0)), hex("grey90"))
  expect_equal(hex(sc$map(4)), hex("darkred"))
})

test_that("midpoint rejects bad settings", {
  m <- fake_avg_exp()
  expect_error(ggRootCellAtlas_expression(m, "GENE1", midpoint = 0, colours = c("blue", "red")),
               "three colours")
  expect_error(ggRootCellAtlas_expression(m, "GENE1", midpoint = 10, c1 = 0, c2 = 5),
               "between")
})

test_that("unknown gene gives an informative error", {
  expect_error(ggRootCellAtlas_expression(fake_avg_exp(), "NOPE"), "NOPE")
})

test_that("group names that match nothing give an informative error", {
  m <- matrix(1:3, nrow = 1, dimnames = list("GENE1", paste0("cluster", 0:2)))
  expect_error(ggRootCellAtlas_expression(m, "GENE1"), "match")
})

test_that("expression does not leave objects in the global environment", {
  clear_leaked_maps()
  before <- globals_snapshot()
  ggRootCellAtlas_expression(fake_avg_exp(), "GENE1")
  expect_identical(globals_snapshot(), before)
})
