library(ggplot2)
library(dplyr)
library(forcats)
library(MetBrewer)
library(tidyr)

setwd("Ciliate_symbiosis/symbiont")
df <- read.table("Functional_enrichment.txt", header=T)

data <- separate_longer_delim(df, cols = Associated_groups, delim = ",")

colors_signac <- met.brewer("Signac", type = "discrete")
colors_pissaro <- met.brewer("Pissaro", type = "discrete")
colors_isfahan <- met.brewer("Isfahan1", type = "discrete", direction = -1)

combined_base <- c(colors_signac, colors_isfahan)

stepped_ramp <- colorRampPalette(combined_base)
cols <- stepped_ramp(29)

#cols <- met.brewer(name = "Signac", n = 31, type = "continuous")

data_sorted <- data %>%
  arrange(tolower(Category), Function) %>%
  mutate(Function = fct_inorder(Function))
  
pdf("Functional_enrichment_categories_dotplot.pdf", width=21.5, height=15)
ggplot(data_sorted, aes(x = Associated_groups, y = Function, size = Enrichment_score, color = Category)) + geom_point(alpha=0.8) + scale_size_continuous(range = c(1,5)) + theme_classic() + theme(text = element_text(size = 15)) + labs(x = "Group", y = "Function") + scale_x_discrete(limits=c("FRESH-SALT210", "BLAMEL4B", "SALT16B", "FARO3", "BI27A", "JUMA2M"))+ scale_y_discrete(limits = rev) + labs(color = "Category", size = "Enrichment score") + scale_color_manual(values = cols) + guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()
