map_names <- c(
  "ggPm.At.longroot.longitudinal",
  "ggPm.At.root.crosssection.m1",
  "ggPm.At.root.crosssection.m2",
  "ggPm.At.root.crosssection.t",
  "ggPm.At.root.crosssection.e1",
  "ggPm.At.root.crosssection.e2",
  "ggPm.At.root.crosssection.d"
)

annotation_columns <- c("SubCellTypes", "CellTypes", "TissueSubTypes",
                        "TissueTypes", "Zones", "Sections", "Atlas",
                        "Atlas_reduced")

# Load the bundled maps into a private environment (never the global one)
load_maps <- function() {
  env <- new.env()
  utils::data(list = map_names, package = "ggRootCellAtlas", envir = env)
  mget(map_names, envir = env)
}

annotation_levels <- function(column) {
  lv <- unique(unlist(lapply(load_maps(), `[[`, column)))
  sort(lv[!is.na(lv)])
}

# Fake AverageExpression() output: genes x groups, deterministic values
fake_avg_exp <- function(column = "Atlas", genes = c("GENE1", "GENE2")) {
  lv <- annotation_levels(column)
  matrix(seq(0, 5, length.out = length(genes) * length(lv)),
         nrow = length(genes), byrow = TRUE,
         dimnames = list(genes, lv))
}

# Seurat v5 AverageExpression() replaces "_" with "-" in group names
seurat_v5_names <- function(m) {
  colnames(m) <- gsub("_", "-", colnames(m))
  m
}

# Any R colour as upper-case "#RRGGBB", to compare with ggplot's fills
hex <- function(x) {
  toupper(grDevices::rgb(t(grDevices::col2rgb(x)), maxColorValue = 255))
}

panel_data <- function(plot, i) ggplot2::ggplot_build(plot[[i]])$data[[1]]

globals_snapshot <- function() sort(ls(globalenv(), all.names = TRUE))

# Earlier calls may already have leaked maps into the global environment;
# clear them so a leak-check starts from a clean state.
clear_leaked_maps <- function() {
  suppressWarnings(rm(list = map_names, envir = globalenv()))
}
