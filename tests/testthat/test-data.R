test_that("all seven maps ship with the package", {
  maps <- load_maps()
  expect_length(maps, 7)
  for (m in maps) expect_s3_class(m, "data.frame")
})

test_that("root_maps() returns the bundled maps in drawing order", {
  expect_identical(root_maps(), load_maps())
})

test_that("maps have the columns the plotting functions rely on", {
  required <- c("ROI.name", annotation_columns, "ROI.id", "point", "x", "y")
  for (nm in map_names) {
    m <- load_maps()[[nm]]
    expect_true(all(required %in% names(m)), info = nm)
    expect_true(is.numeric(m$x) && is.numeric(m$y), info = nm)
    expect_false(anyNA(m$x) || anyNA(m$y), info = nm)
  }
})

test_that("every polygon has at least three points", {
  for (nm in map_names) {
    n <- table(load_maps()[[nm]]$ROI.id)
    expect_true(all(n >= 3), info = nm)
  }
})

test_that("annotations are constant within each polygon", {
  for (nm in map_names) {
    m <- load_maps()[[nm]]
    for (col in annotation_columns) {
      per_roi <- tapply(m[[col]], m$ROI.id, function(v) length(unique(v)))
      expect_true(all(per_roi == 1), info = paste(nm, col))
    }
  }
})

test_that("each cross-section is labelled with its own section", {
  # Regression: the e2 map once carried e1 labels, so the elongation zone
  # showed the same expression twice on the composite plot.
  maps <- load_maps()
  for (s in c("m1", "m2", "t", "e1", "e2", "d")) {
    m <- maps[[paste0("ggPm.At.root.crosssection.", s)]]
    expect_identical(unique(stats::na.omit(m$Sections)), s, info = s)
    suffixed <- grep("_", stats::na.omit(m$Atlas), value = TRUE)
    expect_true(all(endsWith(suffixed, paste0("_", s))), info = s)
  }
})

test_that("no two cross-sections carry identical annotation", {
  maps <- load_maps()[-1]
  for (i in seq_along(maps)[-1]) for (j in seq_len(i - 1)) {
    expect_false(identical(maps[[i]]$Atlas, maps[[j]]$Atlas),
                 info = paste(names(maps)[i], names(maps)[j]))
  }
})

test_that("group names stay distinguishable after Seurat name mangling", {
  # Names are matched loosely (case, '_', '-', '.', ' ' ignored), so two
  # groups must never collapse onto the same key.
  key <- function(x) tolower(gsub("[^[:alnum:]]", "", x))
  for (col in annotation_columns) {
    lv <- annotation_levels(col)
    expect_false(anyDuplicated(key(lv)) > 0, info = col)
  }
})
