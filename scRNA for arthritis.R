# Yan et al. Pathological tissue remodeling fibroblasts 
# CIA scRNA-seq data 
library(Seurat)
library(ggplot2)
library(RColorBrewer)
setwd("/Users/Yan/Desktop/")

# load data
data <- Read10X(data.dir = "/Users/Yan/Desktop/")
Seuratdata <- CreateSeuratObject(counts = data, project = "CtrlaggreCIA", min.cells = 3, min.features = 200)
Seuratdata

#QC
Seuratdata[["percent.mt"]] <- PercentageFeatureSet(Seuratdata, pattern = "^mt-")
VlnPlot(Seuratdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Seuratdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seuratdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Filter
Seuratdata <- subset(Seuratdata, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
Seuratdata

# Normalization
Seuratdata <- NormalizeData(Seuratdata, normalization.method = "LogNormalize", scale.factor = 10000)

Seuratdata <- FindVariableFeatures(Seuratdata, selection.method = "vst", nfeatures = 2000)
top15 <- head(VariableFeatures(Seuratdata), 15)
plot1 <- VariableFeaturePlot(Seuratdata)
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE)
plot1 + plot2

all.genes <- rownames(Seuratdata)
Seuratdata <- ScaleData(Seuratdata, features = all.genes)

#PCA
Seuratdata <- RunPCA(Seuratdata, features = VariableFeatures(object = Seuratdata))
print(Seuratdata[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Seuratdata, dims = 1:2, reduction = "pca")
DimPlot(Seuratdata, reduction = "pca")
DimHeatmap(Seuratdata, dims = 1:15, cells = 500, balanced = TRUE)
Seuratdata <- JackStraw(Seuratdata, num.replicate = 100)
Seuratdata <- ScoreJackStraw(Seuratdata, dims = 1:20)
JackStrawPlot(Seuratdata, dims = 1:20)
ElbowPlot(Seuratdata)

#analysis
Seuratdata <- FindNeighbors(Seuratdata, dims = 1:10)
Seuratdata <- FindClusters(Seuratdata, resolution = 0.04)
Seuratdata <- RunUMAP(Seuratdata, dims = 1:10)

DimPlot(Seuratdata, reduction = "umap", label = TRUE)

# cell type annotation
FeaturePlot(Seuratdata,"Ptprc",pt.size = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) 
FeaturePlot(Seuratdata,"Pdpn",pt.size = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Seuratdata,"Fap",pt.size = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) 
FeaturePlot(Seuratdata,"Cdh5",pt.size = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Seuratdata,"Mcam",pt.size = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) 
FeaturePlot(Seuratdata,"Col6a1",pt.size = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) 

# annotation cell type, merge CD45+ cells
cluster.merged <- RenameIdents(object = Seuratdata,'0' = '0', '1' = '0','4' = '0')

pdf("Fig1a.pdf",useDingbats = F, height = 5, width = 6)
DimPlot(cluster.merged, reduction = "umap", label = F, pt.size = 5, cols = c('0'='#F68282','1'='#31C53F','2'='#00A9FF'))
dev.off()

# reanalysis of fibroblast population
Seuratdata <- FindNeighbors(Seuratdata, dims = 1:10)
Seuratdata <- FindClusters(Seuratdata, resolution = 0.3)
Seuratdata <- RunUMAP(Seuratdata, dims = 1:10)

DimPlot(Seuratdata, reduction = "umap", label = TRUE)
VlnPlot(Seuratdata, features = c("Ptprc","Cd3g","Cd19","Cd79a","Cd14","Itgam","Pdpn","Fap","Thy1","Pecam1","Mcam","Hapln1"))

# subset fibroblast clusters c2 and c8 (containing few Hapln1+ chondrocytes)

Mes <- WhichCells(object = Seuratdata, idents = c("2","8"))
Mes<- SetIdent(object = Seuratdata, cells = Mes, value = "mesenchy")
Mes <- subset(Mes, idents = "mesenchy")
Mes <- ScaleData(Mes, verbose = FALSE)
Mes <- RunPCA(Mes, npcs = 20, verbose = FALSE)
Mes <- RunUMAP(Mes, reduction = "pca", dims = 1:10)
Mes <- FindNeighbors(Mes, reduction = "pca", dims = 1:10)
Mes <- FindClusters(Mes, resolution = 0.3)
DimPlot(Mes, reduction = "umap",pt.size = 3)
VlnPlot(Mes, features = c("Pdpn","Thy1","Fap","Prg4","Hapln1"),pt.size = 0)

# deplete contaminated Hapln1+ chondrocyte cluster1
Fib <- WhichCells(object = Mes, idents = c("0","2","3"))
Fib <- SetIdent(object = Mes, cells = Fib, value = "Fibroblasts")
Fib <- subset(Fib, idents = c("Fibroblasts"))
Fib <- ScaleData(Fib, verbose = FALSE)
Fib <- RunPCA(Fib, npcs = 30, verbose = FALSE)
Fib <- RunUMAP(Fib, reduction = "pca", dims = 1:20)
Fib <- FindNeighbors(Fib, dims = 1:20)
Fib <- FindClusters(Fib, resolution = 0.5)

pdf("Fig.1b.pdf", useDingbats = F, height = 5, width = 6)
DimPlot(Fib, reduction = "umap",pt.size = 4, cols = c('0'='#ABA300','1'='#31C53F','2'='#00A9FF','3'='#C77CFF'))
dev.off()

# Fig.1b dot plot
Fib <- RenameIdents(object=Fib,"3"="mFib1")
Fib <- RenameIdents(object=Fib,"1"="mFib2")
Fib <- RenameIdents(object=Fib,"0"="mFib3")
Fib <- RenameIdents(object=Fib,"2"="mFib4")

markers <- c("Pdpn","Fap","Thy1","Prg4","Tnfsf11","Mmp9","Mmp3","Mmp13","Mmp19")

pdf("Fig.1b.dot.pdf", useDingbats = F, height = 5, width = 7)
DotPlot(Fib, features =markers, dot.scale = 10) + RotatedAxis()
dev.off()

# Extended Data Fig.3d dot plot for Ets1
pdf("Extended Data Fig3d.pdf", useDingbats = F, height = 5.5, width = 3.8)
DotPlot(Fib, features =c ("Ets1"), dot.scale = 5) + RotatedAxis()
dev.off()

# Extended Data Fig.3e MSigDB analysis
install_github("YosefLab/VISION")
library(VISION)

# download signatures “c3.tft.tft_legacy.v7.5.1.symbols.gmt” http://www.gsea-msigdb.org/gsea/downloads.jsp
signatures <- c("/Users/Yan/Desktop/c3.tft.tft_legacy.v7.5.1.symbols.gmt")
vision.obj <- Vision(Fib,signatures = signatures)
vision.obj <- analyze(vision.obj)
viewResults(vision.obj)

# Extended Data Fig.4
FeaturePlot(Seuratdata,"Col6a1",pt.size = 3) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) 
VlnPlot(Fib, features = c("Thy1","Col6a1"))

# Reanalysis of scRNA-seq data dbGap Study Accession: phs001529.v1.p1. Stephenson, W. et al. “Single-cell RNA-seq of rheumatoid arthritis synovial tissue using low-cost microfluidic instrumentation”  Nat. Commun. 9, 791, 2018
 
#data import
counts <- read.csv('RA_5Knees_Expression_Matrix.csv', header=TRUE, sep=",")
rownames(counts)<-counts[,1]
counts<-counts[,-1]

#create seurat object
Data1 <- CreateSeuratObject(counts = counts)

#QC 
Data1 <- PercentageFeatureSet(Data1, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(object = Data1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# QC filter
Data1 <- subset(x = Data1, subset = percent.mt <10 & nFeature_RNA <7500 & nCount_RNA <40000)
# Normalization
Data1 <- NormalizeData(Data1, verbose = FALSE)
Data1$dataset <- "data1"

Data1 <- FindVariableFeatures(Data1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Data1), 10)
plot1 <- VariableFeaturePlot(Data1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# scale data
all.genes <- rownames(Data1)
Data1 <- ScaleData(Data1, features = all.genes)
# run PCA to variable features
Data1 <- RunPCA(Data1, features = VariableFeatures(object = Data1))
# check dims and PCA;
ElbowPlot(Data1)

Data1 <- FindNeighbors(Data1, dims = 1:15)
Data1 <- FindClusters(Data1, resolution = 0.1)
Data1 <- RunUMAP(object = Data1, dims = 1:15)

pdf("Extended Data Fig.1b.pdf", useDingbats = F, height = 5, width = 7)
DimPlot(object = Data1, label = T, reduction = "umap", label.size = 3) + theme(aspect.ratio = 0.8) + NoLegend()
dev.off()
#Extended Data Fig.1b
VlnPlot(Data1, features =c("PDPN","FAP"),log=T, pt.size = 0)+ RotatedAxis() 
        
# subset fibroblasts cluster 1 and cluster 3
Fib <- WhichCells(object = Data1, idents = c("1","3"))
Fib<- SetIdent(object = Data1, cells = Fib, value = "Fib")
Fib1 <- subset(Fib, idents = "Fib")
Fib1@assays #7939 cells
        
# Reanalysis of scRNA-seq data ImmPort Accession: SDY998. Zhang, F.  et al. “Defining inflammatory cell states in rheumatoid arthritis joint synovial tissues by integrating single-cell transcriptomics and mass cytometry” Nat. Immunol. 20, 928-942, 2019. 
# data import
umi.counts <- read.table('celseq_matrix_ru1_molecules.tsv.gz', sep = "\t", header=TRUE)
rownames(umi.counts)<-umi.counts[,1]
umi.counts<-umi.counts[,-1]
umi.counts[is.na(umi.counts)]<-0
Meta <- read.table("celseq_meta.tsv.725591.gz", sep = "\t", header = TRUE)
rownames(Meta)<-Meta[,1]

#create seurat object
Data2 <- CreateSeuratObject(counts = umi.counts, min.cells = 3, min.features = 500, meta.data = Meta[,4:5])
Data2$orig.ident <- Data2$type
Idents(object = Data2) <- Data2$type
        
#QC 
Data2[["percent.mt"]] <- PercentageFeatureSet(object = Data2, pattern = "^MT-")
VlnPlot(object = Data2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        
# filter
Data2 <- subset(x = Data2, subset = percent.mt <30 & nFeature_RNA < 9000 & nCount_RNA <60000)
Data2
Data2$dataset <- "data2"
        
# Analysis 
Data2 <- FindVariableFeatures(Data2, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Data2), 10)
plot1 <- VariableFeaturePlot(Data2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
        
# scale data
all.genes <- rownames(Data2)
Data2 <- ScaleData(Data2, features = all.genes)
        
# run PCA 
Data2 <- RunPCA(Data2, features = VariableFeatures(object = Data2))
        
# check dims and PCA; dims = 1: 15
ElbowPlot(Data2)
        
Data2 <- FindNeighbors(Data2, dims = 1:15)
Data2 <- FindClusters(Data2, resolution = 0.1)
Data2 <- RunUMAP(object = Data2, dims = 1:15)
        
Data2@assays #9553 cells
        
pdf("Extended Data Fig.1c.pdf", useDingbats = F, height = 5, width = 7)
DimPlot(object = Data2, label = T, reduction = "umap", label.size = 3) + theme(aspect.ratio = 0.8) + NoLegend()
dev.off()

pdf("Extended Data Fig.1c.Vln.pdf", useDingbats = F, height = 15, width = 10)
VlnPlot(Data2, features =c("FAP","PDPN"),log=T, pt.size = 0)+ RotatedAxis()
dev.off()
        
# subset fibroblast cluster 0 and 5 
Fib2 <- WhichCells(object = Data2, idents = c("0","5"))
Fib2<- SetIdent(object = Data2, cells = Fib2, value = "Fib2")
Fib2 <- subset(Fib2, idents = "Fib2")
Fib2@assays  #3249 fibroblasts
        
# Merge Fib1 and Fib2
merge <- merge(x= Fib1, y = Fib2)
        
#QC
Idents(object = merge) <- merge$dataset
VlnPlot(object = merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(object = merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        
# filter
merge <- subset(x = merge, subset = percent.mt <15 & nFeature_RNA <9000 & nCount_RNA <60000)
        
# SCT
library(sctransform)
merge.list <- SplitObject(object = merge, split.by = "dataset")
for (i in 1:length(x = merge.list)) {merge.list[[i]] <- SCTransform(object = merge.list[[i]], verbose = FALSE,vars.to.regress = "percent.mt")}
merge.features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 3000)
merge.list <- PrepSCTIntegration(object.list = merge.list, anchor.features = merge.features, verbose = FALSE)
merge.anchors <- FindIntegrationAnchors(object.list = merge.list, normalization.method = "SCT", anchor.features = merge.features, verbose = FALSE)
merge.integration <- IntegrateData(anchorset = merge.anchors, normalization.method = "SCT", verbose = FALSE)
        
# PCA, UMAP
merge.integration <- RunPCA(object = merge.integration, verbose = FALSE)
merge.integration <- RunUMAP(object = merge.integration, dims = 1:5)
# FindClusters
merge.integration <- FindClusters(merge.integration, resolution = 0.5)
#cluster merge
cluster.merge <- RenameIdents(object = merge.integration,'3' = '3', '1' = '3')
# rename cluster
cluster.merge <- RenameIdents(object=cluster.merge,"4"="F1")
cluster.merge <- RenameIdents(object=cluster.merge,"5"="F2")
cluster.merge <- RenameIdents(object=cluster.merge,"8"="F3")
cluster.merge <- RenameIdents(object=cluster.merge,"0"="F4")
cluster.merge <- RenameIdents(object=cluster.merge,"6"="F5")
cluster.merge <- RenameIdents(object=cluster.merge,"2"="F6")
cluster.merge <- RenameIdents(object=cluster.merge,"3"="F7")
cluster.merge <- RenameIdents(object=cluster.merge,"7"="F8")
        
# UMAP Extended Data Fig.1d
pdf("Extended Data Fig.1d.pdf", useDingbats = F, height = 7, width = 10)
DimPlot(object = cluster.merge, label = T, reduction = "umap", label.size = 3,pt.size = 2) + theme(aspect.ratio = 0.8) + NoLegend()
dev.off()

# Vlnplot 
pdf("TNFSF11.pdf", useDingbats = F, height = 7, width = 10)
VlnPlot(cluster.merge,features = "TNFSF11",pt.size = 0)

pdf("MMP13.pdf",useDingbats = F, height = 7, width = 10)
VlnPlot(cluster.merge,features = "MMP13",pt.size = 0)

pdf("MMP19.pdf",useDingbats = F, height = 7, width = 10)
VlnPlot(cluster.merge,features = "MMP19",pt.size = 0)

#Extended Data Fig.3c
pdf("Extended Data Fig3c.pdf", useDingbats = F, height = 5.5, width = 3.8)
DotPlot(cluster.merge,features = "ETS1",dot.scale = 5) + RotatedAxis()
dev.off()
