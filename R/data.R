#' Arabidopsis root maps
#'
#' ggPlantmap maps of the *Arabidopsis thaliana* root: one longitudinal section
#' and six cross-sections. Each row is one vertex of a cell polygon.
#' Use [root_maps()] to get all seven as a list.
#'
#' The cross-sections are `m1`/`m2` (early/late meristem), `t` (transition
#' zone), `e1`/`e2` (early/late elongation zone) and `d` (differentiation zone).
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{ROI.name}{Name of the cell (region of interest).}
#'   \item{SubCellTypes, CellTypes, TissueSubTypes, TissueTypes}{Cell identity
#'     at increasing levels of grouping.}
#'   \item{Zones}{Developmental zone: `"m"`, `"t"`, `"e"`, `"d"`, `"Root cap"`
#'     or `"SCN"` (stem cell niche).}
#'   \item{Sections}{Section the cell belongs to, e.g. `"m1"` or `"e2"`.}
#'   \item{Atlas, Atlas_reduced}{Cell type combined with section (`Atlas`) or
#'     zone (`Atlas_reduced`), matching the single-cell atlas clusters.}
#'   \item{ROI.id}{Polygon identifier.}
#'   \item{point}{Order of the vertex within the polygon.}
#'   \item{x, y}{Vertex coordinates.}
#' }
#' @aliases ggPm.At.longroot.longitudinal ggPm.At.root.crosssection.m1
#'   ggPm.At.root.crosssection.m2 ggPm.At.root.crosssection.t
#'   ggPm.At.root.crosssection.e1 ggPm.At.root.crosssection.e2
#'   ggPm.At.root.crosssection.d
#' @usage data(ggPm.At.root.crosssection.m1)
#' @docType data
#' @keywords datasets
#' @name root_map_data
NULL
