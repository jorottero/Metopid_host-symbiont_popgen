library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)

setwd("Ciliate_symbiosis/symbiont")

data <- t(read.table("kegg-metabolism-module_pathwise_completeness-MATRIX_modified.txt", header=T, row.names=1))
data2 <- t(read.table("kegg-metabolism-module_stepwise_completeness-MATRIX_modified.txt", header=T, row.names=1))
data3 <- t(read.table("kegg-metabolism-enzyme_hits-MATRIX_modified.txt", header=T, row.names=1))

pop <- read.table("metadata.txt", header=T)

dist <- vegdist(data, method="euclidean", binary=F, na.rm=TRUE)
dist[is.na(dist)] = 0
pcoa <- pcoa(dist, correction = "cailliez")
dist2 <- vegdist(data2, method="euclidean", binary=F, na.rm=TRUE)
dist2[is.na(dist2)] = 0
pcoa2 <- pcoa(dist2, correction = "cailliez")
dist3 <- vegdist(data3, method="jaccard", binary=T, na.rm=TRUE)
dist3[is.na(dist3)] = 0
pcoa3 <- pcoa(dist3, correction = "cailliez")

Axis1 <- pcoa$vectors[,1]
Axis2 <- pcoa$vectors[,2]
var1 <- round(pcoa$values[1,2], digits=4)*100
var2 <- round(pcoa$values[2,2], digits=4)*100
x_axis <- paste("Axis1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis2"," ","[",var2,"%]",sep="")
Axis1b <- pcoa2$vectors[,1]
Axis2b <- pcoa2$vectors[,2]
var1b <- round(pcoa2$values[1,2], digits=4)*100
var2b <- round(pcoa2$values[2,2], digits=4)*100
x_axisb <- paste("Axis1"," ","[",var1b,"%]",sep="")
y_axisb <- paste("Axis2"," ","[",var2b,"%]",sep="")
Axis1c <- pcoa3$vectors.cor[,1]
Axis2c <- pcoa3$vectors.cor[,2]
var1c <- round(pcoa3$values[1,3], digits=4)*100
var2c <- round(pcoa3$values[2,3], digits=4)*100
x_axisc <- paste("Axis1"," ","[",var1c,"%]",sep="")
y_axisc <- paste("Axis2"," ","[",var2c,"%]",sep="")

plot1 <- cbind(data, Axis1, Axis2, pop)
plot2 <- cbind(data2, Axis1b, Axis2b, pop)
plot3 <- cbind(data3, Axis1c, Axis2c, pop)

plot1$Location <- factor(plot1$Location, levels = c("FRESH22", "FRESH24", "FRESH25B", "FRESH26", "BLAMEL4B", "BI27A", "JUMA2M", "SALT210", "SALT16B", "FARO3"))
plot2$Location <- factor(plot2$Location, levels = c("FRESH22", "FRESH24", "FRESH25B", "FRESH26", "BLAMEL4B", "BI27A", "JUMA2M", "SALT210", "SALT16B", "FARO3"))
plot3$Location <- factor(plot3$Location, levels = c("FRESH22", "FRESH24", "FRESH25B", "FRESH26", "BLAMEL4B", "BI27A", "JUMA2M", "SALT210", "SALT16B", "FARO3"))

col <- c("FRESH25B" = "#8E6CAD", "FRESH22" = "#A057A0", "FRESH24" = "#9B457C", "FRESH26" = "#615496", "SALT210" = "#3A7FA5", "SALT16B" = "#3AB1DB", "FARO3" = "#F4B664", "BLAMEL4B" = "#A1C666", "BI27A" = "#39B54A", "JUMA2M" = "#28815B")

pdf("KEGG_module_pathwise_completeness_PCoA.12.pdf", width=10, height=5)
p <- ggplot(plot1, mapping = aes(x = Axis1, y = Axis2, color=Location)) + geom_point(position=position_jitter(width=0.008, height=0.008), aes(fill = factor(Location)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, aesthetics = c("colour", "fill"), name="Location") + labs(x = x_axis, y = y_axis)
p
dev.off()

pdf("KEGG_module_stepwise_completeness_PCoA.12.pdf", width=10, height=5)
p <- ggplot(plot2, mapping = aes(x = Axis1b, y = Axis2b, color=Location)) + geom_point(position=position_jitter(width=0.008, height=0.008), aes(fill = factor(Location)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, aesthetics = c("colour", "fill"), name="Location") + labs(x = x_axisb, y = y_axisb)
p
dev.off()

pdf("KEGG_enzyme_hits_presence_PCoA.12.pdf", width=10, height=5)
p <- ggplot(plot3, mapping = aes(x = Axis1c, y = Axis2c, color=Location)) + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Location)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, aesthetics = c("colour", "fill"), name="Location") + labs(x = x_axisc, y = y_axisc)
p
dev.off()

df <- data.frame(pop)

sink("PERMANOVA-KEGG.txt")
# Adonis test
adonis2(dist ~ Location, data = df, add = "cailliez", permutations = 999, by = NULL)
adonis2(dist2 ~ Location, data = df, add = "cailliez", permutations = 999, by = NULL)
adonis2(dist3 ~ Location, data = df, add = "cailliez", permutations = 999, by = NULL)

# Dispersion test
beta1 <- betadisper(dist, df$Location)
permutest(beta1)
beta2 <- betadisper(dist2, df$Location)
permutest(beta2)
beta3 <- betadisper(dist3, df$Location)
permutest(beta3)
sink()