library(ggplot2)
library(vegan)
library(gridExtra)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dplyr)

setwd("Ciliate_symbiosis/symbiont")
data <- read.table("StrainFrequencies.txt", header=T, row.names=1)
mat <- as.matrix(data)
pop <- read.table("sample_metadata.txt", header=T)
col = colorRamp2(c(0, 1), hcl_palette = "Blues", reverse=TRUE)

ha <- HeatmapAnnotation(Location = pop$Location, annotation_height = unit(4, "mm"), col = list(Location = c("FRESH25B" = "#8E6CAD", "FRESH22" = "#A057A0", "FRESH24" = "#9B457C", "FRESH26" = "#615496", "SALT210" = "#3A7FA5", "SALT16B" = "#3AB1DB", "FARO3" = "#F4B664", "BLAMEL4B" = "#A1C666", "BI27A" = "#39B54A", "JUMA2M" = "#28815B")))

bray_dist <- function(x) { vegdist(x, method = "bray") }

pdf("Strain_heatmap.pdf", height=2.5, width=8)
map <- Heatmap(mat, border = "black", top_annotation = ha, km = 1, name = "Strain proportion", col = col, show_row_names = TRUE, cluster_rows = FALSE, clustering_distance_columns = bray_dist, clustering_method_columns = "average", show_column_names = FALSE, show_column_dend = TRUE, column_dend_reorder = TRUE)
draw(map, annotation_legend_side = "right", merge_legend = TRUE)
dev.off()
