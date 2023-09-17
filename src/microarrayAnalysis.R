setwd("E:/learning/bioinf_advance1")
setRepositories()
1 2
### then Esc to exit
library("RColorBrewer")
library (limma)
library(Biobase)
library(reshape2)
library(ggplot2)
library(gplots)
library(plyr)
library(pheatmap)
library (GEOquery)
series <- "GSE9476"
platform <- "GPL96"
  gset <- getGEO( series, GSEMatrix = TRUE , AnnotGPL = TRUE , destdir = "Data/")
  if ( length(gset) > 1) idx <- grep( platform , attr( gset, "names")) else idx <- 1
  gset <- gset[[idx]]
#### if there is only one platform you can just type "gset <- gset[[1]]  
  class(gset)
  ### it must return ExpressionSet
  #### now grouping of data
  gr <- c ( "CD34", rep("BM", 10), rep("CD34", 7), rep("AML", 26), rep("PB", 10), rep("CD34", 10))
length(gr)  
###64 which means the work was done correctly
ex <- exprs(gset)
###if tha max was very high means it is not in log2 scale and you can turn it to log2 scale
max(ex)
min(ex)
### if it is not in log 2 scale, 
ex <- log2 (ex + 1)
####now draw a boxplot to see they are normalised or not
### to adjust the plot you can add width=64 to make it bigger
pdf("Results/boxplot.pdf", width = 64)
boxplot(ex)
dev.off()
#### this dataset is normalised, but if it was not we could normalised it with the following command
ex<- normalizeQuantiles(ex)
#### now we want to do quality control with correlation heatmap
pdf("Results/corHeatmap.pdf", width = 15 , height = 15)
pheatmap(cor(ex))
dev.off()
#### this heatmap is just a visual representation of the correlation matrix which you can see it
#excor <- cor (ex)
#excor[1:5, 1:5]
###the problem here is that we don't know what is the name of coloumns and rows
#### pheatmap command have some options like labels_row which we can use
#####you can also omit borders and add different colours, like greenred or bluered, bluered is better
######as some people maybe colourblind

pdf("Results/corHeatmap.pdf", width = 15 , height = 15)
pheatmap(cor(ex), labels_row = gr , labels_col = gr, border_color = NA, color = bluered(256))
dev.off()
###################################
########################################
# Principal component analysis
########################################
##################################
pc <- prcomp(ex)
pdf("Results/PCA.pdf")
plot(pc)
dev.off()
# if you want to see the PC
names(pc)
# sdev  rotation center scale x , x: matrix after being in new PCA, so they are genes in PCA
##ROTATION is samples in reduced space or PCA
###if dim(pc$x)
dim(pc$x)
#but if you see colnames(pc$x)
colnames(pc$x)
# you see that it is not samples, it is the PC1 to PC64, so instead of seeing all of them you just see
##1 and 2
plot(pc$x[ , 1:2])
#every dot is a gene, pc1 shows the most variation and after that pc2
##the genes in either side of pc1 are high or 0 expressed in every sample
### so we need to modify it, as we need genes that show variation
#### so we subtract each gene to its avarage expression
##### therefore, from now on we just have differences
######but we do not do it on our expression vector, we create a alaki matrix and do this
ex.scale <- t(scale(t(ex), scale = FALSE ))
PC <- prcomp(ex.scale)
pdf("Results/pc_scaled.pdf")
plot(PC)
plot(PC$x[ ,1:2])
dev.off()
########
##########
#PCA of samples in pc$rotation or pc$r
pcr <- data.frame(pc$r[ , 1:3], Group= gr)
pdf("Results/PCA_samples.pdf")
### in aes() for x and y you just can give your coloumns, so if the data that you need is not in coloumns
### you need to change your data
ggplot( pcr, aes( x= PC1 , y= PC2, color= Group)) + geom_point()
dev.off()
## as class pc$rotation is matrix and if you add string to matrix it will turn numbers to strings we
## change it to data.frame
## as no other variable strats with r we can tell pc$r instead of pc$rotation
## The other reason we changed matrix to data.frame is that ggplot only works with data.frame

###################################
########################################
#Differential Expression analysis
##################################
#########################################
gr <- factor(gr)
gset$description <- factor(gr)
## design matrix is a simple matrix that shows each sample belongs to which group, head(design) and see
### it for yourself
design <- model.matrix( ~ description + 0 , gset)
colnames(design) <- levels(gr)
## based on design matrix it fits a linear model to your data
fit <- lmFit(gset , design)
count.matrix <- makeContrasts(AML-CD34, levels = design)
fit2<- contrasts.fit(fit, count.matrix)
fit2<- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = Inf)
tT <- subset(tT, select = c("Gene.symbol", "Gene.ID" , "adj.P.Val", "logFC" ))
write.table(tT, "Results/AML-CD34.txt", row.names = FALSE, sep = "\t", quote = FALSE)
################ UP and DOWN genes
aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(aml.up$Gene.symbol)
aml.up.genes <- sub("///.*", "", aml.up.genes)
write.table(aml.up.genes, file = "Results/AML-CD34_Up.txt", quote = FALSE, row.names = FALSE , col.names = FALSE)
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(aml.down$Gene.symbol)
aml.down.genes <- sub("///.*", "", aml.down.genes)
write.table(aml.down.genes, file = "Results/AML-CD34_Down.txt", quote = FALSE, row.names = FALSE , col.names = FALSE)
