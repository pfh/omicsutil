
library(Seurat)
library(purrr)

MyDimPlot <- function(seu, ...) {
    DimPlot(seu, ..., shuffle=TRUE) + coord_fixed() + theme_void()
}

MyFeaturePlot <- function(seu, genes, ..., ncol=3) {
    plots <- map(genes,\(g) FeaturePlot(seu, g, ...) + coord_fixed() + theme_void())
    cowplot::plot_grid(plotlist=plots, ncol=ncol)
}