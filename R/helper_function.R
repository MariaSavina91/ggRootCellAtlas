# ggPlantmap imports stringr in its NAMESPACE but does not list it in its
# DESCRIPTION, so installing ggPlantmap alone leaves it unloadable. Importing
# stringr here makes sure it gets installed alongside this package.
#' @importFrom stringr str_detect
NULL

#' Generate a Common Palette
#'
#' Creates a consistent color palette for a given grouping variable across datasets.
#'
#' @param Group_name Character. Name of the group column.
#' @param ... Data frames to be included in the palette generation.
#' @param color_palette A palette function taking the number of colours
#'   (default is `scales::hue_pal()`), or a character vector of colours:
#'   unnamed colours are assigned to the groups in order, named colours are
#'   matched to the group names.
#' @return A named vector of colors.
#' @examples
#' generate_common_palette("TissueTypes", root_maps()[[1]])
#' generate_common_palette("TissueTypes", root_maps()[[1]],
#'                         color_palette = okabe_ito_pal())
#' @importFrom stats setNames
#' @export

generate_common_palette <- function(Group_name, ..., color_palette = scales::hue_pal()) {
  groups <- unlist(lapply(list(...), function(data) data[[Group_name]]))
  unique_groups <- unique(groups[!is.na(groups)])
  n <- length(unique_groups)

  if (is.function(color_palette)) {
    return(setNames(color_palette(n), unique_groups))
  }
  if (!is.character(color_palette)) {
    stop("The palette must be a function or a character vector of colours.", call. = FALSE)
  }
  if (!is.null(names(color_palette))) {
    missing <- setdiff(unique_groups, names(color_palette))
    if (length(missing) > 0) {
      stop(sprintf("The palette has no colour for: %s.", paste(missing, collapse = ", ")),
           call. = FALSE)
    }
    return(color_palette[unique_groups])
  }
  if (length(color_palette) < n) {
    stop(sprintf("The palette has %d colours but \"%s\" has %d groups.",
                 length(color_palette), Group_name, n), call. = FALSE)
  }
  setNames(color_palette[seq_len(n)], unique_groups)
}

#' Colour-blind safe palette
#'
#' Palette function based on the Okabe-Ito colours, extended to 11 colours.
#' For more groups than that it falls back to `scales::hue_pal()`, as no
#' palette keeps that many colours distinguishable.
#'
#' @return A function that takes the number of colours and returns them.
#' @examples
#' okabe_ito_pal()(5)
#' @export

okabe_ito_pal <- function() {
  colours <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
               "#CC79A7", "#6A3D9A", "#B15928", "#999999", "#E8E8E8")
  function(n) {
    if (n <= length(colours)) colours[seq_len(n)] else scales::hue_pal()(n)
  }
}

# Maps in the order they are drawn: longitudinal section first, then the
# cross-sections from the root tip upwards.
map_names <- c(
  "ggPm.At.longroot.longitudinal",
  "ggPm.At.root.crosssection.m1",
  "ggPm.At.root.crosssection.m2",
  "ggPm.At.root.crosssection.t",
  "ggPm.At.root.crosssection.e1",
  "ggPm.At.root.crosssection.e2",
  "ggPm.At.root.crosssection.d"
)

#' Root maps bundled with the package
#'
#' Returns the longitudinal section and the six cross-sections of the
#' *Arabidopsis thaliana* root as a named list of ggPlantmap data frames.
#'
#' @return A named list of 7 data frames, in drawing order.
#' @examples
#' names(root_maps())
#' @export

root_maps <- function() {
  env <- new.env(parent = emptyenv())
  utils::data(list = map_names, package = "ggRootCellAtlas", envir = env)
  mget(map_names, envir = env)
}

check_annotation_column <- function(column, maps, arg) {
  if (!is.character(column) || length(column) != 1) {
    stop(sprintf("`%s` must be a single column name.", arg), call. = FALSE)
  }
  if (!column %in% names(maps[[1]])) {
    available <- setdiff(names(maps[[1]]), c("ROI.id", "point", "x", "y"))
    stop(sprintf("`%s` = \"%s\" is not a column of the root maps. Use one of: %s.",
                 arg, column, paste(available, collapse = ", ")), call. = FALSE)
  }
}

# Composite layout: longitudinal section (1) on the left, cross-sections
# stacked on the right from the tip (2) at the bottom to d (7) at the top.
root_layout <- "
  1#
  17
  17
  16
  16
  15
  15
  14
  14
  13
  13
  12
  12
  1#
  1#
  "
