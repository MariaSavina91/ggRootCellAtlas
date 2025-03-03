#' Generate Root Cell Atlas Annotation Visualizations
#'
#' @description This function creates annotated visualizations for various root sections.
#' It uses a shared color palette to represent group annotations and combines multiple
#' plots into a single composite visualization.
#'
#' # Example usage
#' ggRootCellAtlas_annotation("Zones")
#' ggRootCellAtlas_annotation("Sections")
#' ggRootCellAtlas_annotation("SubCellTypes")
#' ggRootCellAtlas_annotation("CellTypes")
#' ggRootCellAtlas_annotation("TissueTypes")
#' ggRootCellAtlas_annotation("TissueSubTypes")
#' ggRootCellAtlas_annotation("Atlas")
#' ggRootCellAtlas_annotation("Atlas_reduce")

#'
#' @param Group_name Character. Name of the column in the datasets that represents
#'   the group annotations to be visualized.
#'
#' @return  A ggplot2 object representing the combined annotated plots, with each plot
#'   corresponding to a specific Arabidopsis thaliana root section and
#'   a shared legend for group annotations.
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr tibble
#' @importFrom stats setNames
#' @import ggPlantmap
#'
#' @export
#'


ggRootCellAtlas_annotation <- function(Group_name) {

  # load("data/ggPm.At.longroot.longitudinal.rda")
  # load("data/ggPm.At.root.crosssection.m1.rda")
  # load("data/ggPm.At.root.crosssection.m2.rda")
  # load("data/ggPm.At.root.crosssection.t.rda")
  # load("data/ggPm.At.root.crosssection.e1.rda")
  # load("data/ggPm.At.root.crosssection.e2.rda")
  # load("data/ggPm.At.root.crosssection.d.rda")

  # Generate color palette
  palette <- generate_common_palette(Group_name,
                                     ggPm.At.longroot.longitudinal,
                                     ggPm.At.root.crosssection.m1,
                                     ggPm.At.root.crosssection.m2,
                                     ggPm.At.root.crosssection.t,
                                     ggPm.At.root.crosssection.e1,
                                     ggPm.At.root.crosssection.e2,
                                     ggPm.At.root.crosssection.d
  )

  # Generate plots
  p1 <- ggPlantmap.plot(data = ggPm.At.longroot.longitudinal, eval(parse(text = Group_name)),show.legend = FALSE) +
    scale_fill_manual(values = palette, na.value = "black", name = Group_name)
  p2 <- ggPlantmap.plot(data = ggPm.At.root.crosssection.m1, eval(parse(text = Group_name)),show.legend = FALSE) +
    scale_fill_manual(values = palette, na.value = "black", name = Group_name)
  p3 <- ggPlantmap.plot(data = ggPm.At.root.crosssection.m2, eval(parse(text = Group_name)),show.legend = FALSE) +
    scale_fill_manual(values = palette, na.value = "black", name = Group_name)
  p4 <- ggPlantmap.plot(data = ggPm.At.root.crosssection.t, eval(parse(text = Group_name)),show.legend = FALSE) +
    scale_fill_manual(values = palette, na.value = "black", name = Group_name)
  p5 <- ggPlantmap.plot(data = ggPm.At.root.crosssection.e1, eval(parse(text = Group_name)),show.legend = FALSE) +
    scale_fill_manual(values = palette, na.value = "black", name = Group_name)
  p6 <- ggPlantmap.plot(data = ggPm.At.root.crosssection.e2, eval(parse(text = Group_name)),show.legend = FALSE) +
    scale_fill_manual(values = palette, na.value = "black", name = Group_name)
  p7 <- ggPlantmap.plot(data = ggPm.At.root.crosssection.d, eval(parse(text = Group_name)),show.legend = FALSE) +
    scale_fill_manual(values = palette, na.value = "black", name = Group_name)

  legend_data <- data.frame(
    x = rep(1, times = nrow(as.data.frame(palette))),
    y = rep(1, times = nrow(as.data.frame(palette))),
    label = row.names(as.data.frame(palette)),
    color = as.data.frame(palette)$palette
  )
  p8 <- ggplot(legend_data) +
    geom_point(aes(x = x, y = y, color = color), shape = 15, size = 0) +  # For plot
    scale_color_manual(values = legend_data$color, labels = legend_data$label) +
    guides(color = guide_legend(override.aes = list(size = 5))) +  # For legend
    theme_void() +
    theme(legend.position = "left")+
    labs(color = Group_name)

  # Combine plots
  combined <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8

  # Design layout
  design <- ("
  1#8
  178
  178
  168
  168
  158
  158
  148
  148
  138
  138
  128
  128
  1#8
  1#8
  ")

  # Final plot
  final_plot <- combined + plot_layout(design = design)

  return(final_plot)
}
