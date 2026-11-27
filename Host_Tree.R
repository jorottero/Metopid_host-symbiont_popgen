library(phangorn)
library(ggtree)
library(ggplot2)
library(ape)

setwd("Ciliate_symbiosis/host")

data <- read.table("tree_metadata.txt", header=T)
data$Group <- factor(data$Group, levels = c("FRESH22", "FRESH24", "FRESH25B", "FRESH26", "BLAMEL4B", "BI27A", "JUMA2M", "SALT210", "SALT16B", "SALT15A", "FARO3", "Outgroup"))

col <- c("FRESH25B" = "#8E6CAD", "FRESH22" = "#A057A0", "FRESH24" = "#9B457C", "FRESH26" = "#615496", "SALT210" = "#3A7FA5", "SALT16B" = "#3AB1DB", "SALT15A" = "#CAE1FF", "FARO3" = "#F4B664", "BLAMEL4B" = "#A1C666", "BI27A" = "#39B54A", "JUMA2M" = "#28815B", "Outgroup" = "black")

tree <- read.tree("Astral-Pro3.tre")
rooted_tree <- midpoint(tree)
write.tree(rooted_tree, file="Astral-Pro3.midpoint.tre")

pdf("Astral-Pro3.pdf")
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0,0.05) + geom_tiplab(align=TRUE, hjust=-.15, size=3) + annotate("point", x=0, y=20, shape=21, fill="darkgray", color="black", size=3) + annotate("text", x=0.002, y=20, label = "> 75% node support", hjust = "left") + geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 0.75), size=3, shape=21, fill="darkgray", color="black") + geom_treescale(0,25)
t2 <- t %<+% data + geom_tippoint(aes(color=factor(Group)), shape=20, size=1) + theme(legend.position = c(0.1,0.7)) + guides(color = guide_legend(override.aes = list(label = "\u25A1", size = 0))) + aes(color=factor(Group)) + scale_color_manual(values = col, name="Group", na.value="black")
t2
dev.off()

