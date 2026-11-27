library(vegan)
library(ggplot2)
library(ape)
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)
library(ggrepel)

setwd("Ciliate_symbiosis/symbiont")

data <- read.table("Methanocorpusculum.Freebayes.FINAL.CT.FORMAT", header=T, row.names=1)
com <- t(apply(data, 2, function(i) i/sum(i, na.rm=TRUE)))
data2 <- t(read.table("Methanocorpusculum.Freebayes.FINAL.GT.FORMAT", header=T, row.names=1))
pop <- read.table("Methanocorpusculum.clst", header=T)

dist <- vegdist(com, method="bray", binary=F, na.rm=TRUE)
matrix <- df2genind(data2, ploidy = 1, ncode = 1)
dist2 <- vegdist(matrix, method="euclidean", binary=F, na.rm=TRUE)
dist2[is.na(dist2)] = 0
pcoa <- pcoa(dist, correction = "cailliez")
pcoa2 <- pcoa(dist2, correction = "cailliez")

Axis1 <- pcoa$vectors[,1]
Axis2 <- pcoa$vectors[,2]
var1 <- round(pcoa$values[1,2], digits=4)*100
var2 <- round(pcoa$values[2,2], digits=4)*100
x_axis <- paste("Axis1"," ","[",var1,"%]",sep="")
y_axis <- paste("Axis2"," ","[",var2,"%]",sep="")
Axis1b <- pcoa2$vectors.cor[,1]
Axis2b <- pcoa2$vectors.cor[,2]
var1b <- round(pcoa2$values[1,3], digits=4)*100
var2b <- round(pcoa2$values[2,3], digits=4)*100
x_axisb <- paste("Axis1"," ","[",var1b,"%]",sep="")
y_axisb <- paste("Axis2"," ","[",var2b,"%]",sep="")

plot1 <- cbind(com, Axis1, Axis2, pop)
plot2 <- cbind(matrix, Axis1b, Axis2b, pop)

plot1$Location <- factor(plot1$Location, levels = c("FRESH22", "FRESH24", "FRESH25B", "FRESH26", "BLAMEL4B", "BI27A", "JUMA2M", "SALT210", "SALT16B", "FARO3"))
plot2$Location <- factor(plot2$Location, levels = c("FRESH22", "FRESH24", "FRESH25B", "FRESH26", "BLAMEL4B", "BI27A", "JUMA2M", "SALT210", "SALT16B", "FARO3"))

col <- c("FRESH25B" = "#8E6CAD", "FRESH22" = "#A057A0", "FRESH24" = "#9B457C", "FRESH26" = "#615496", "SALT210" = "#3A7FA5", "SALT16B" = "#3AB1DB", "FARO3" = "#F4B664", "BLAMEL4B" = "#A1C666", "BI27A" = "#39B54A", "JUMA2M" = "#28815B")

# Colored by Location
pdf("Methanocorpusculum_PCoA.12.pdf")
p <- ggplot(plot1, aes(Axis1, Axis2, color= Location)) + geom_point(position=position_jitter(width=0.008, height=0.008), aes(fill = factor(Location)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, name="Location") + labs(x = x_axis, y = y_axis)
p
dev.off()

pdf("Methanocorpusculum_PCoA.geno.12.pdf")
p <- ggplot(plot2, aes(Axis1b, Axis2b, color= Location)) + geom_point(position=position_jitter(width=0.008, height=0.008), aes(fill = factor(Location)), shape=21, color="gray32", size=6) + theme_classic() + theme(text = element_text(size = 15)) + scale_fill_manual(values = col, name="Location") + labs(x = x_axisb, y = y_axisb)
p
dev.off()

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