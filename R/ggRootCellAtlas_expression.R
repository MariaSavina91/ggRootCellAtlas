#' Generate Root Cell Atlas Gene Expression Heatmaps
#'
#' @description This function creates a series of heatmaps representing the expression of a specified gene
#' or other cell feature across various root sections. It combines the heatmaps into a single composite
#' visualization using a shared color gradient.
#'
#' Group names are matched to the root map annotation ignoring case and
#' punctuation, so names altered by Seurat (e.g. `"Cortex-m1"`, `"Young.LRC"`)
#' still match `"Cortex_m1"` and `"Young LRC"`.
#'
#' @param avg_exp Average expression with genes in rows and groups in columns:
#'   a matrix, sparse matrix or data frame, or the list returned by
#'   `Seurat::AverageExpression()` (its first assay is used).
#' @param Gene A string specifying the gene (row of `avg_exp`) to visualize.
#' @param c1 Numeric value specifying the minimum value for the color gradient. Defaults to 0,
#'   or with `midpoint`, to a range symmetric around the midpoint.
#' @param c2 Numeric value specifying the maximum value for the color gradient. Defaults to the
#'   maximum expression of `Gene` across all sections. Values outside `c1`..`c2`
#'   are drawn in the colour of the nearest limit.
#' @param Annotation Character. Name of the root map column that the columns of
#'   `avg_exp` correspond to (the Idents of the Seurat object). Defaults to `"Atlas"`.
#' @param colours Colours of the gradient, from low to high. Defaults to
#'   snow2, yellow, red3, red4; with `midpoint`, to blue, white, red.
#' @param midpoint Optional value at which a diverging scale is centred, e.g. 0
#'   for fold changes or PCA scores. The middle colour of `colours` is placed at
#'   this value. Needs at least three colours.
#' @param na.colour Colour of cells without a value. Defaults to dark grey.
#'
#' @return A patchwork object combining the heatmaps, one per root section,
#'   with a shared legend.
#'
#' @examples
#' # Toy data: one value per "Atlas" group
#' groups <- unique(unlist(lapply(root_maps(), `[[`, "Atlas")))
#' groups <- groups[!is.na(groups)]
#' avg_exp <- matrix(seq_along(groups), nrow = 1,
#'                   dimnames = list("GeneID", groups))
#' ggRootCellAtlas_expression(avg_exp, "GeneID")
#' ggRootCellAtlas_expression(avg_exp, "GeneID", c1 = 10, c2 = 40)
#' ggRootCellAtlas_expression(avg_exp, "GeneID", colours = c("white", "darkgreen"))
#'
#' # Diverging scale, e.g. for log2 fold changes
#' lfc <- avg_exp - mean(avg_exp)
#' ggRootCellAtlas_expression(lfc, "GeneID", midpoint = 0)
#'
#' \dontrun{
#' # With a Seurat object
#' Seurat.object <- SetIdent(Seurat.object, value = "Atlas")
#' avg_exp <- AverageExpression(Seurat.object, assays = "RNA")
#' ggRootCellAtlas_expression(avg_exp, "GeneID")
#' }
#'
#' @import ggplot2
#' @import patchwork
#' @import ggPlantmap
#'
#' @export

ggRootCellAtlas_expression <- function(avg_exp, Gene, c1 = NA, c2 = NA, Annotation = "Atlas",
                                       colours = c("snow2", "yellow", "red3", "red4"),
                                       midpoint = NULL, na.colour = "grey30") {

  maps <- root_maps()
  check_annotation_column(Annotation, maps, "Annotation")

  # Extract Gene expression as a named vector: group -> value
  exp <- as_expression_matrix(avg_exp)
  if (!is.character(Gene) || length(Gene) != 1 || !Gene %in% rownames(exp)) {
    stop(sprintf("Gene \"%s\" is not a row of `avg_exp`.", paste(Gene, collapse = ", ")),
         call. = FALSE)
  }
  values <- exp[Gene, ]
  groups <- unique(unlist(lapply(maps, `[[`, Annotation)))
  groups <- groups[!is.na(groups)]
  names(values) <- match_group_names(names(values), groups)
  values <- values[!is.na(names(values))]
  if (length(values) == 0) {
    stop(sprintf(paste0("None of the columns of `avg_exp` match the \"%s\" groups ",
                        "of the root maps (e.g. %s)."),
                 Annotation, paste(utils::head(groups, 3), collapse = ", ")),
         call. = FALSE)
  }

  # Merge data
  maps <- lapply(maps, function(map) {
    map$Gene <- unname(values[map[[Annotation]]])
    map
  })

  # Common colour limits for all panels
  if (is.null(midpoint)) {
    if (is.na(c1)) {
      c1 <- 0
    }
    if (is.na(c2)) {
      c2 <- max(values, na.rm = TRUE)
    }
    stops <- NULL
  } else {
    if (missing(colours)) {
      colours <- c("#2166AC", "white", "#B2182B")
    }
    if (length(colours) < 3) {
      stop("A diverging scale (`midpoint`) needs at least three colours.", call. = FALSE)
    }
    spread <- max(abs(values - midpoint), na.rm = TRUE)
    if (spread == 0) spread <- 1
    if (is.na(c1)) {
      c1 <- midpoint - spread
    }
    if (is.na(c2)) {
      c2 <- midpoint + spread
    }
    if (!(c1 < midpoint && midpoint < c2)) {
      stop(sprintf("`midpoint` (%s) must lie between c1 (%s) and c2 (%s).", midpoint, c1, c2),
           call. = FALSE)
    }
    stops <- diverging_stops(length(colours), (midpoint - c1) / (c2 - c1))
  }

  # Generate plots
  plots <- lapply(maps, function(map) {
    ggPlantmap.heatmap(map, .data$Gene) +
      scale_fill_gradientn(colours = colours, values = stops,
                           limits = c(c1, c2), oob = scales::squish,
                           na.value = na.colour) +
      labs(fill = Gene)
  })

  wrap_plots(plots, design = root_layout, guides = "collect")
}

# Positions (0..1) of n colours so that the middle one sits at `mid`, with
# the others spread evenly on either side of it
diverging_stops <- function(n, mid) {
  k <- ceiling(n / 2)
  c(seq(0, mid, length.out = k), seq(mid, 1, length.out = n - k + 1)[-1])
}

# Accepts matrix / sparse matrix / data frame / Seurat AverageExpression() list
as_expression_matrix <- function(avg_exp) {
  if (is.list(avg_exp) && !is.data.frame(avg_exp)) {
    if (length(avg_exp) > 1) {
      message(sprintf("Using the first assay of `avg_exp` (\"%s\").", names(avg_exp)[1]))
    }
    avg_exp <- avg_exp[[1]]
  }
  exp <- as.matrix(avg_exp)
  if (is.null(rownames(exp)) || is.null(colnames(exp))) {
    stop("`avg_exp` needs gene names as row names and groups as column names.",
         call. = FALSE)
  }
  exp
}

# Map names as found in the data to the annotation groups, ignoring case and
# punctuation. A common "assay." prefix (as added by as.data.frame() on the
# Seurat list, e.g. "RNA.Cortex_m1") is dropped when that matches more names.
match_group_names <- function(x, groups) {
  key <- function(s) tolower(gsub("[^[:alnum:]]", "", s))
  matched <- groups[match(key(x), key(groups))]
  prefix <- unique(sub("\\..*$", "", x))
  if (length(prefix) == 1 && all(grepl(".", x, fixed = TRUE))) {
    stripped <- groups[match(key(sub("^[^.]*\\.", "", x)), key(groups))]
    if (sum(!is.na(stripped)) > sum(!is.na(matched))) matched <- stripped
  }
  matched
}
