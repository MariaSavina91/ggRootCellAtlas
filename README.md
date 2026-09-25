# ggRootCellAtlas

Plot cell annotations and gene expression on a map of the *Arabidopsis thaliana*
root: one longitudinal section and six cross-sections, from the meristem (m1, m2)
through the transition (t) and elongation (e1, e2) zones to the differentiation
zone (d). Built on [ggPlantmap](https://github.com/leonardojo/ggPlantmap).

## Installation

```r
# install.packages("remotes")
remotes::install_github("MariaSavina91/ggRootCellAtlas")
```

## Usage

Colour the root by an annotation column (`"SubCellTypes"`, `"CellTypes"`,
`"TissueSubTypes"`, `"TissueTypes"`, `"Zones"`, `"Sections"`, `"Atlas"` or
`"Atlas_reduced"`):

```r
library(ggRootCellAtlas)
ggRootCellAtlas_annotation("TissueTypes")
```

Show the average expression of a gene per cell group, e.g. from Seurat:

```r
Seurat.object <- SetIdent(Seurat.object, value = "Atlas")
avg_exp <- AverageExpression(Seurat.object, assays = "RNA")
ggRootCellAtlas_expression(avg_exp, "AT1G01010")

# fixed colour range; values outside it take the colour of the nearest limit
ggRootCellAtlas_expression(avg_exp, "AT1G01010", c1 = 0, c2 = 5)
```

`avg_exp` can also be a plain matrix or data frame with genes in rows and groups
in columns. Group names are matched ignoring case and punctuation, so names
changed by Seurat (`"Cortex-m1"`, `"Young.LRC"`) still match the map.

## Colours

```r
# own gradient (default: snow2 -> yellow -> red3 -> red4)
ggRootCellAtlas_expression(avg_exp, "AT1G01010", colours = c("white", "darkgreen"))

# diverging scale centred on 0 (blue -> white -> red), e.g. for fold changes
ggRootCellAtlas_expression(lfc, "AT1G01010", midpoint = 0)

# annotation colours: a palette function, or a vector (named by group or in order)
ggRootCellAtlas_annotation("TissueTypes", palette = c(
  Epidermis = "#E69F00", "Ground tissue" = "#009E73", Stele = "#CC79A7",
  "Root cap" = "#56B4E9", SCN = "#D55E00"
))
```

The default annotation palette, `okabe_ito_pal()`, is colour-blind safe for up
to 11 groups (`Zones`, `Sections`, `TissueTypes`, `TissueSubTypes`). Columns
with more groups fall back to the standard ggplot2 hues.

Cells without a value or annotation are drawn in dark grey; change this with
`na.colour`.

The maps themselves are available with `root_maps()`.
