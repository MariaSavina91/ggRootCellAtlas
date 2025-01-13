#' Generate a Common Palette
#'
#' Creates a consistent color palette for a given grouping variable across datasets.
#'
#' @param Group_name Character. Name of the group column.
#' @param ... Data frames to be included in the palette generation.
#' @param color_palette Function. A color palette function (default is `scales::hue_pal`).
#' @return A named vector of colors.
#' @import scales
#' @importFrom dplyr filter
#' @export

generate_common_palette <- function(Group_name, ..., color_palette = scales::hue_pal()) {
  clean_data <- lapply(list(...), function(data) data[!is.na(data[,Group_name]), ])
  unique_groups <- unique(unlist(lapply(clean_data, function(data) unique(data[,Group_name]))))
  n_groups <- length(unique_groups)
  color_palette <- color_palette(n_groups)
  return(setNames(color_palette, unique_groups))
}
