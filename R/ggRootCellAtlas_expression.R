#' Generate Root Cell Atlas Gene Expression Heatmaps
#'
#' @description This function creates a series of heatmaps representing the expression of a specified gene
#' or oother cell fiature across various root sections. It combines the heatmaps into a single composite
#' visualization using a shared color gradient.
#'
#' # Example usage
#' Seurat.object<-SetIdent(Seurat.object,value = "AnnotationName")
#' avg_exp <- AverageExpression(Seurat.object,assays = 'RNA',slot = 'data')
#' ggRootCellAtlas_expression(avg_exp, "GeneID", 1, 10, "AnnotationName")
#'
#' ggRootCellAtlas_expression(avg_exp, Gene = "GeneID")
#'
#' @param avg_exp A data frame containing average gene expression data. Columns should represent
#'   different atlas regions, and rows should represent genes.
#' @param Gene A string specifying the gene for which expression or other feature data should be visualized.
#' @param c1 Numeric value specifying the minimum value for the color gradient. Defaults to 0.
#' @param c2 Numeric value specifying the maximum value for the color gradient. Defaults to the
#'   maximum gene expression across all datasets.
#' @param Annotation Character. Name of the column to annotate that corresponds to Idents of Seurat object.
#' Defaults to the "Atlas" Idents
#'
#' @return A `ggplot2` object representing the combined heatmaps, with each heatmap corresponding
#'   to a specific atlas dataset and a shared legend indicating the gene expression levels.
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr tibble
#' @importFrom stats setNames
#' @import ggPlantmap
#'
#' @export




ggRootCellAtlas_expression <- function(avg_exp,Gene,c1=NA,c2=NA, Annotation="Atlas") {

  # load("data/ggPm.At.longroot.longitudinal.rda")
  # load("data/ggPm.At.root.crosssection.m1.rda")
  # load("data/ggPm.At.root.crosssection.m2.rda")
  # load("data/ggPm.At.root.crosssection.t.rda")
  # load("data/ggPm.At.root.crosssection.e1.rda")
  # load("data/ggPm.At.root.crosssection.e2.rda")
  # load("data/ggPm.At.root.crosssection.d.rda")

  # Extract Gene expression
  a<-as.data.frame(avg_exp)
  colnames(a)
  for (i in 1: length(colnames(a))) {
    colnames(a)[i] <- substr(colnames(a)[i], regexpr("\\.", colnames(a)[i]) + 1, nchar(colnames(a)[i]))
    colnames(a)[i] <- gsub("G1.", "G1/", colnames(a)[i])
    colnames(a)[i] <- gsub("G2.", "G2/", colnames(a)[i])
    colnames(a)[i] <- gsub("\\.", " ", colnames(a)[i])
    colnames(a)[i] <- gsub(" m1", "_m1", colnames(a)[i])
    colnames(a)[i] <- gsub(" m2", "_m2", colnames(a)[i])
    colnames(a)[i] <- gsub(" t", "_t", colnames(a)[i])
    colnames(a)[i] <- gsub(" e1", "_e1", colnames(a)[i])
    colnames(a)[i] <- gsub(" e2", "_e2", colnames(a)[i])
    colnames(a)[i] <- gsub(" d", "_d", colnames(a)[i])
  }
  colnames(a)
  a<-as.data.frame(t(a))
  # a$Atlas<-row.names(a)

  b<-data.frame(matrix(
    vector(), nrow(a), 2, dimnames=list(c(), c(Annotation, "Gene"))),
    stringsAsFactors=F)
  b[,1]<-row.names(a)
  b$Gene<-a[,Gene]


  # Merge data
  ggPm.At.longroot.longitudinal <- ggPlantmap.merge(ggPm.At.longroot.longitudinal,b,Annotation)
  ggPm.At.root.crosssection.m1 <- ggPlantmap.merge(ggPm.At.root.crosssection.m1,b,Annotation)
  ggPm.At.root.crosssection.m2 <- ggPlantmap.merge(ggPm.At.root.crosssection.m2,b,Annotation)
  ggPm.At.root.crosssection.t <- ggPlantmap.merge(ggPm.At.root.crosssection.t,b,Annotation)
  ggPm.At.root.crosssection.e1 <- ggPlantmap.merge(ggPm.At.root.crosssection.e1,b,Annotation)
  ggPm.At.root.crosssection.e2 <- ggPlantmap.merge(ggPm.At.root.crosssection.e2,b,Annotation)
  ggPm.At.root.crosssection.d <- ggPlantmap.merge(ggPm.At.root.crosssection.d,b,Annotation)

  # Generate color palette
  common_max<-max(max(as.numeric(ggPm.At.longroot.longitudinal$Gene),na.rm = TRUE),
                  max(as.numeric(ggPm.At.root.crosssection.m1$Gene),na.rm = TRUE),
                  max(as.numeric(ggPm.At.root.crosssection.m2$Gene),na.rm = TRUE),
                  max(as.numeric(ggPm.At.root.crosssection.t$Gene),na.rm = TRUE),
                  max(as.numeric(ggPm.At.root.crosssection.e1$Gene),na.rm = TRUE),
                  max(as.numeric(ggPm.At.root.crosssection.e2$Gene),na.rm = TRUE),
                  max(as.numeric(ggPm.At.root.crosssection.d$Gene),na.rm = TRUE)
  )
  if (is.na(c1)){
    c1<-0
  }

  if(is.na(c2)){
    c2<-common_max
  }

  # Generate plots
  p1 <- ggPlantmap.heatmap(ggPm.At.longroot.longitudinal, Gene)+
    scale_fill_gradientn(colours = c("snow2","yellow","red3","red4"),limits = c(c1, c2),na.value = "black")+#0, common_max),na.value = "black")+
    labs(fill = Gene)
  p2 <- ggPlantmap.heatmap(ggPm.At.root.crosssection.m1, Gene)+
    scale_fill_gradientn(colours = c("snow2","yellow","red3","red4"),limits = c(c1, c2),na.value = "black")+#0, common_max),na.value = "black")+
    labs(fill = Gene)
  p3 <- ggPlantmap.heatmap(ggPm.At.root.crosssection.m2, Gene)+
    scale_fill_gradientn(colours = c("snow2","yellow","red3","red4"),limits = c(c1, c2),na.value = "black")+#0, common_max),na.value = "black")+
    labs(fill = Gene)
  p4 <- ggPlantmap.heatmap(ggPm.At.root.crosssection.t, Gene)+
    scale_fill_gradientn(colours = c("snow2","yellow","red3","red4"),limits = c(c1, c2),na.value = "black")+#0, common_max),na.value = "black")+
    labs(fill = Gene)
  p5 <- ggPlantmap.heatmap(ggPm.At.root.crosssection.e1, Gene)+
    scale_fill_gradientn(colours = c("snow2","yellow","red3","red4"),limits = c(c1, c2),na.value = "black")+#0, common_max),na.value = "black")+
    labs(fill = Gene)
  p6 <- ggPlantmap.heatmap(ggPm.At.root.crosssection.e2, Gene)+
    scale_fill_gradientn(colours = c("snow2","yellow","red3","red4"),limits = c(c1, c2),na.value = "black")+#0, common_max),na.value = "black")+
    labs(fill = Gene)
  p7 <- ggPlantmap.heatmap(ggPm.At.root.crosssection.d, Gene)+
    scale_fill_gradientn(colours = c("snow2","yellow","red3","red4"),limits = c(c1, c2),na.value = "black")+#0, common_max),na.value = "black")+
    labs(fill = Gene)

  # Combine plots
  combined <- p1 + p2 + p3 + p4 + p5 + p6 + p7

  # Design layout
  design <- ("
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
  ")

  # Final plot
  final_plot <- combined + plot_layout(design = design, guides = "collect")+ theme(legend.position = "right")#, legend.title = paste(Group_name))

  return(final_plot)
}

