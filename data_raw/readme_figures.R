# Renders the example figures used in README.md into man/figures/.
# Expression values are toy data, made up to show the colour scales.
#
# Run from the package root: Rscript data_raw/readme_figures.R

library(ggRootCellAtlas)

save_figure <- function(plot, name, width = 5) {
  ggplot2::ggsave(file.path("man/figures", name), plot,
                  width = width, height = 5.5, dpi = 90, bg = "white")
}

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

groups <- unique(unlist(lapply(root_maps(), `[[`, "Atlas")))
groups <- groups[!is.na(groups)]
section <- sub("^.*_", "", groups)

# Toy gene rising from the meristem to the late elongation zone
level <- c(m1 = 0.5, m2 = 1, t = 2, e1 = 3.5, e2 = 5, d = 4)
toy <- matrix(ifelse(section %in% names(level), level[section], 0.2), nrow = 1,
              dimnames = list("GeneX", groups))

# Toy fold change: down in the meristem, up in the elongation zone
lfc <- toy - 2.5

save_figure(ggRootCellAtlas_annotation("TissueTypes"), "README-annotation.png", width = 6)
save_figure(ggRootCellAtlas_expression(toy, "GeneX"), "README-expression.png")
save_figure(ggRootCellAtlas_expression(lfc, "GeneX", midpoint = 0), "README-diverging.png")
