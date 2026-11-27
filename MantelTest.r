library(vegan)
library(ggplot2)
library(ncf)
library(tidyverse)
library(ape)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)
library(geosphere)

setwd("Ciliate_symbiosis/Congruence_testing")

host <- as.matrix(read.table("Host_patristic_distance.txt", header = TRUE, row.names=1))
sym <- as.matrix(read.table("Symbiont_patristic_distance.txt", header = TRUE, row.names=1))

sink("MantelTestPatristic.txt")
mantel.test(host, sym, resamp = 9999, method="spearman", quiet = FALSE)
sink()
