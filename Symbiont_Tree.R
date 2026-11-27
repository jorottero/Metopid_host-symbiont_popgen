library(phangorn)
library(ggtree)
library(ggplot2)
library(ape)

setwd("Ciliate_symbiosis/symbiont")

data <- read.table("tree_metadata.txt", header=T)
data$Group <- factor(data$Group, levels = c("FRESH22", "FRESH24", "FRESH25B", "FRESH26", "BLAMEL4B", "BI27A", "JUMA2M", "JUMA2P", "SALT210", "SALT16B", "SALT15A", "FARO3", "Outgroup"))

col <- c("FRESH25B" = "#8E6CAD", "FRESH22" = "#A057A0", "FRESH24" = "#9B457C", "FRESH26" = "#615496", "SALT210" = "#3A7FA5", "SALT16B" = "#3AB1DB", "SALT15A" = "#CAE1FF", "FARO3" = "#F4B664", "BLAMEL4B" = "#A1C666", "BI27A" = "#39B54A", "JUMA2M" = "#28815B", "JUMA2P" = "#2F4F4F", "Outgroup" = "black")

tree <- read.tree("RAxML_bipartitions.core_genes.tre")
rooted_tree <- midpoint(tree)

pdf("RAxML_bipartitions.core_genes.pdf", width=10, height=10)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0,0.14) + geom_tiplab(align=TRUE, hjust=-.15, size=3) + annotate("point", x=0, y=37, shape=21, fill="darkgray", color="black", size=3) + annotate("text", x=0.003, y=37, label = "> 75% node support", hjust = "left") + geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 75), size=3, shape=21, fill="darkgray", color="black") + geom_treescale(0,40)
t2 <- t %<+% data + geom_tippoint(aes(color=factor(Group)), shape=20, size=1) + theme(legend.position = c(0.1,0.8)) + guides(color = guide_legend(override.aes = list(label = "\u25A1", size = 0))) + aes(color=factor(Group)) + scale_color_manual(values = col, name="Group", na.value="black")
t2
dev.off()

t1 <- drop.tip(tree, "Methanocorpusculum_labreanum", rooted=FALSE, trim.internal=FALSE)
t2 <- drop.tip(t1, "WGA_72", rooted=FALSE, trim.internal=FALSE)
t3 <- drop.tip(t2, "WGA_73", rooted=FALSE, trim.internal=FALSE)
t4 <- drop.tip(t3, "WGA_74", rooted=FALSE, trim.internal=FALSE)
t5 <- drop.tip(t4, "WGA_75", rooted=FALSE, trim.internal=FALSE)
t6 <- drop.tip(t5, "WGA_76", rooted=FALSE, trim.internal=FALSE)
t7 <- drop.tip(t6, "WGA_77", rooted=FALSE, trim.internal=FALSE)

t8 <- midpoint(t7)
write.tree(t8, file="RAxML_bipartitions.core_genes.midpoint.reduced.tre")

