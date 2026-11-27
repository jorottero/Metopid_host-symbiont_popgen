library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)

setwd("Ciliate_symbiosis/host")

#IBS: Identical by descent/state, average genotype distance between individuals, ignores allele frequencies, used in MDS
data <- as.matrix(read.table("F25B-A1A10.ibsMat"))
data[is.na(data)] = 0 
#COV: Covariance matrix, weighting scheme based on allele frequencies, used in PCA
data2 <- as.matrix(read.table("F25B-A1A10.covMat"))
data2[is.na(data2)] = 0
pop <- read.table("F25B-A1A10.txt", header=T)
pcoa <- pcoa(data, correction = "cailliez")
pca <- prcomp(data2, scale=TRUE, center=TRUE)

Axis1 <- pcoa$vectors.cor[,1]
Axis2 <- pcoa$vectors.cor[,2]
var1 <- round(pcoa$values[1,3], digits=4)*100
var2 <- round(pcoa$values[2,3], digits=4)*100
x_axis <- paste("Axis1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis2"," ","[",var2,"%]",sep="")
PC1 <- pca$rotation[,1]
PC2 <- pca$rotation[,2]
eigs <- pca$sdev^2
pc1_var <- round(eigs[1]/sum(eigs), digits=4)*100
pc2_var <- round(eigs[2]/sum(eigs), digits=4)*100
pc1_axis <- paste("PC1"," ","[",pc1_var,"%]",sep="")
pc2_axis <- paste("PC2"," ","[",pc2_var,"%]",sep="")
plot1 <- cbind(data, Axis1, Axis2, pop)
plot2 <- cbind(data2, PC1, PC2, pop)

plot1$Location <- factor(plot1$Location, levels = c("FRESH22", "FRESH24", "FRESH25B", "FRESH26", "BLAMEL4B", "BI27A", "JUMA2M", "SALT210", "SALT16B", "SALT15A", "FARO3"))
plot2$Location <- factor(plot2$Location, levels = c("FRESH22", "FRESH24", "FRESH25B", "FRESH26", "BLAMEL4B", "BI27A", "JUMA2M", "SALT210", "SALT16B",  "SALT15A", "FARO3"))

col <- c("FRESH25B" = "#8E6CAD", "FRESH22" = "#A057A0", "FRESH24" = "#9B457C", "FRESH26" = "#615496", "SALT210" = "#3A7FA5", "SALT16B" = "#3AB1DB", "SALT15A" = "#CAE1FF", "FARO3" = "#F4B664", "BLAMEL4B" = "#A1C666", "BI27A" = "#39B54A", "JUMA2M" = "#28815B")

pdf("F25B-A1A10_PCoA.12.pdf")
p <- ggplot(plot1, aes(Axis1, Axis2, color=Location)) + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Location)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, aesthetics = c("colour", "fill"), name="Location") + labs(x = x_axis, y = y_axis)
p
dev.off()

pdf("F25B-A1A10_PCA.12.pdf")
p <- ggplot(plot2, aes(PC1, PC2, color=Location)) + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Location)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, aesthetics = c("colour", "fill"), name="Location") + labs(x = pc1_axis, y = pc2_axis)
p
dev.off()

dist <- as.dist(data)
dist2 <- dist(data2, method="euclidean")
df <- data.frame(pop)

sink("PERMANOVA.txt")
# Adonis test
adonis2(dist ~ Location, data = df, add = "cailliez", permutations = 999, by = NULL)
adonis2(dist2 ~ Location, data = df, add = "cailliez", permutations = 999, by = NULL)

# Dispersion test
beta1 <- betadisper(dist, df$Location)
permutest(beta1)
beta2 <- betadisper(dist2, df$Location)
permutest(beta2)
sink()
