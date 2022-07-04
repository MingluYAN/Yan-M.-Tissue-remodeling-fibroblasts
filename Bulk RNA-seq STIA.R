# Bulk RNA-seq analysis, count files were processed using CLC genomics Workbench v.12.
rm(list = ls())
source("http://Bioconductor.org/biocLite.R")
library(ggplot2)
library(DESeq2)
library(pheatmap)
setwd("/Users/apple/Desktop/")

# read count files
aSTIA <- read.delim("a4countsSTIA.txt", header = TRUE, row.names = 1)
# analyze 
countmatrix <- as.matrix(aSTIA[1:6])
head(countmatrix)
table <- data.frame(name = c("6ko", "3ko", "1ko", "5con", "4con", "2con"), condition = c("2", "2", "2", "1", "1", "1"))
dds2 <- DESeqDataSetFromMatrix (countmatrix, colData = table, design = ~condition)
dds2 <- DESeq(dds2)
dds2
colData (dds2)

res = results (dds2, contrast = c("condition", "2", "1" ))
res = res [order(res$pvalue),]
head(res)
write.csv (res, file = "STIAcountsDeseq.csv")
res = res [order(res$pvalue),]

# Fig.5c
RANKL <-read.csv ("Osteoclatogenesis.csv",row.names = 1)
RANKL <- as.matrix(RANKL)
pdf(file="RANKL.pdf",width = 15,height = 15)
pheatmap(RANKL, scale = "row", cluster_cols=FALSE, cellwidth = 15, cellheight = 15, fontsize = 10)
dev.off()

MMP <-read.csv ("MMP.csv",row.names = 1)
MMP <- as.matrix(MMP)
pdf(file="MMP.pdf",width = 15,height = 15)
pheatmap(MMP, scale = "row", cluster_cols=FALSE, cellwidth = 15, cellheight = 15, fontsize = 10)
dev.off()

ECM <-read.csv ("ECM.csv",row.names = 1)
ECM <- as.matrix(ECM)
pdf(file="ECM.pdf",width = 15,height = 15)
pheatmap(ECM, scale = "row", cluster_cols=FALSE, cellwidth = 15, cellheight = 15, fontsize = 10)
dev.off()

cytokine <-read.csv ("Cytokines and chemokines.csv",row.names = 1)
cyto <- as.matrix(cytokine)
pdf(file="cyto.pdf",width = 15,height = 15)
pheatmap(cyto, scale = "row", cluster_cols=FALSE, cellwidth = 15, cellheight = 15, fontsize = 10)
dev.off()

surface <-read.csv ("Surface molecules.csv",row.names = 1)
surface <- as.matrix(surface)
pdf(file="surface.pdf",width = 15,height = 15)
pheatmap(surface, scale = "row", cluster_cols=FALSE, cellwidth = 15, cellheight = 15, fontsize = 10)
dev.off()


Proliferation <-read.csv ("Proliferation.csv",row.names = 1)
Proliferation <- as.matrix(Proliferation)
pdf(file="proliferation.pdf",width = 15,height = 15)
pheatmap(Proliferation, scale = "row", cluster_cols=FALSE, cellwidth = 15, cellheight = 15, fontsize = 10)
dev.off()

