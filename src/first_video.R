setwd("C:/Users/mahdi.nohtani/Desktop/Bioinformatics_tutorial")

library(GEOquery)
library(Biobase)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
series <-"GSE39132"
platform <- "GPL570"


gset <- getGEO(series, AnnotGPL = TRUE, GSEMatrix = TRUE, destdir = "./data")

gset <- gset[[1]]
class(gset)
ex <- exprs(gset)
gr <- c(rep("norm", 5), rep("tumour", 5))


min(ex)
max(ex)
head(ex)
pdf("./results/boxplot.pdf")
boxplot(ex)
dev.off()
x <- cor(ex)
pdf("./results/corheatmap.pdf")
pheatmap(x, labels_col = gr, labels_row = gr)
dev.off()
