# Reanalysis of scRNA-seq data GSE136103 Ramachandran, P. et al.  “Resolving the fibrotic niche of human liver cirrhosis at single-cell level. Nature 2019 Nov;575(7783):512-518”. 
# load data, A: cirrhotic1 CD45-A; B: cirrhotic1 CD45-B; L2: cirrhotic2 CD45-; L3: cirrhotic3 CD45-; 

A <- Read10X(data.dir = "/Users/Yan/Desktop/")
data1A <- CreateSeuratObject(counts = A, project = "1A", min.cells = 3, min.features = 200)
B <- Read10X(data.dir = "/Users/Yan/Desktop/")
data1B <- CreateSeuratObject(counts = B, project = "1B", min.cells = 3, min.features = 200)
L2 <- Read10X(data.dir = "/Users/Yan/Desktop/")
data2 <- CreateSeuratObject(counts = L2, project = "2", min.cells = 3, min.features = 200)
L3 <- Read10X(data.dir = "/Users/Yan/Desktop/")
data3 <- CreateSeuratObject(counts = L3, project = "3", min.cells = 3, min.features = 200)

View(data1A@meta.data)
View(data1B@meta.data)
View(data2@meta.data)
View(data3@meta.data)

data1A$group <- "1A"
data1B$group <- "1B"
data2$group <- "L2"
data3$group <- "L3"

scRNAlist <- list()
scRNAlist[[1]] = NormalizeData(data1A)
scRNAlist[[1]] = FindVariableFeatures(scRNAlist[[1]],selection.method = "vst")
scRNAlist[[2]] = NormalizeData(data1B)
scRNAlist[[2]] = FindVariableFeatures(scRNAlist[[2]],selection.method = "vst")
scRNAlist[[3]] = NormalizeData(data2)
scRNAlist[[3]] = FindVariableFeatures(scRNAlist[[3]],selection.method = "vst")
scRNAlist[[4]] = NormalizeData(data3)
scRNAlist[[4]] = FindVariableFeatures(scRNAlist[[4]],selection.method = "vst")


scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = 20000)
scRNA <- IntegrateData(anchorset = scRNA.anchors)

Liver <- scRNA
proj_name <- data.frame(proj_name=rep("Liver",ncol(Liver)))
rownames(proj_name) <- row.names(Liver@meta.data)
Liver <- AddMetaData(Liver,proj_name)

DefaultAssay(Liver) <- "RNA"

VlnPlot(Liver, features = c("nFeature_RNA", "nCount_RNA","percent.mt"))
Liver1 <- subset(Liver, subset = nCount_RNA < 6000)
Liver1 <- FindVariableFeatures(Liver1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Liver1), 10)
all.genes = rownames(Liver1)
Liver1 = ScaleData(Liver1, features = all.genes)
Liver1 <- RunPCA(Liver1, features = VariableFeatures(object = Liver1))
Liver1 <- JackStraw(Liver1, num.replicate = 100)
Liver1  <- ScoreJackStraw(Liver1 ,dims = 1:20)
JackStrawPlot(Liver1 ,dims = 1:20)
ElbowPlot(Liver1)

Liver1  <- FindNeighbors(Liver1,dims = 1:20)
Liver1@assays #21294 features for 6898 cells
Liver1  <- FindClusters(Liver1 , resolution = 0.5)
Liver1  <- RunUMAP(Liver2 , dims = 1:20)

# Extended Data Fig.10a
DimPlot(Liver1 , reduction = "umap", label = T)

# Extended Data Fig.10b

pdf("Liver PTPRC.pdf", useDingbats = F, height = 5, width = 6)
FeaturePlot(Liver1, features = "PTPRC")
dev.off()

pdf("Liver PECAM1.pdf", useDingbats = F, height = 5, width = 6)
FeaturePlot(Liver1, features = "PECAM1")
dev.off()

pdf("Liver EPCAM.pdf", useDingbats = F, height = 5, width = 6)
FeaturePlot(Liver1, features = "EPCAM")
dev.off()

pdf("Liver ALB.pdf", useDingbats = F, height = 5, width = 6)
FeaturePlot(Liver1, features = "ALB")
dev.off()

pdf("Liver PDGFRB.pdf", useDingbats = F, height = 5, width = 6)
FeaturePlot(Liver1, features = "PDGFRB")
dev.off()

pdf("Liver ACTA2.pdf", useDingbats = F, height = 5, width = 6)
FeaturePlot(Liver1, features = "ACTA2")
dev.off()

pdf("Liver TAGLN.pdf", useDingbats = F, height = 5, width = 6)
FeaturePlot(Liver1, features = "TAGLN")
dev.off()


pdf("Liver RGS5.pdf", useDingbats = F, height = 5, width = 6)
FeaturePlot(Liver1, features = "RGS5")
dev.off()

# Subset fibroblasts cluster 10 and 13
Fib.liver <- subset(Liver1, idents = c("10","13"))

#Extended Data Fig.10c

pdf("Fib.MYH11.pdf", useDingbats = F, height = 3, width = 5)
FeaturePlot(Fib.liver,"MYH11",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

pdf("Fib.RGS5.pdf", useDingbats = F, height = 3, width = 5)
FeaturePlot(Fib.liver,"RGS5",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

pdf("Fib.LUM.pdf", useDingbats = F, height = 3, width = 5)
FeaturePlot(Fib.liver,"LUM",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

pdf("Fib.COL3A1.pdf", useDingbats = F, height = 3, width = 5)
FeaturePlot(Fib.liver,"COL3A1",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

pdf("Fib.PDGFRA.pdf", useDingbats = F, height = 3, width = 5)
FeaturePlot(Fib.liver,"PDGFRA",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

pdf("Fib.COL1A1.pdf", useDingbats = F, height = 3, width = 5)
FeaturePlot(Fib.liver,"COL1A1",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

#Extended Data Fig.10d
FeaturePlot(Fliver, features = c("ETS1", "SPI1"), pt.size = 1.2, blend = T, max.cutoff = 1,cols = c("black", "red", "blue"))

#Extended Data Fig.10i
FeatureScatter(object = Fliver,feature1 = 'ETS1', feature2 = 'SPI1',plot.cor = TRUE, shuffle =T, jitter = TRUE, pt.size = 2.5)

# Reanalysis of scRNA-seq data GSE135893 Habermann, A. C. et al.  “Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial and mesenchymal lineages in pulmonary fibrosis”. Sci Adv 2020 Jul;6(28):eaba1972.
# load data
ild <- Read10X(data.dir = "/Users/Yan/Desktop/", gene.column = 1)
ild <- CreateSeuratObject(ild, project = "scRNA lung", min.features = 200)
ild <- PercentageFeatureSet(object = ild, pattern = "^MT-", col.name = "percent.mt")

# QC
pdf("QC.pdf", useDingbats = F, height = 15, width = 50)
VlnPlot(ild, features = c("nFeatureRNA","nCountRNA", "percent.mt"), ncol = 3)
dev.off()

# Filter out genes with less than 1000 features and more than 25% MT
ild <- subset(ild, subset = nFeature_RNA > 1000&percent.mt < 25)

# Run SCTansform
ild <- SCTransform(object = ild, verbose = T)

# Chose a PC for analysis
ild <- FindVariableFeatures(ild, verbose = T, nfeatures = 3000)
ild <- ScaleData(ild, features = row.names(ild@assays$SCT@data))
ild <- RunPCA(ild)
ElbowPlot(ild)

# RunUMAP, FindClusters, FindNeighbors dims =20, res =0.01

ild <- RunUMAP(object = ild, dims = 1:20, verbose = F)
ild <- FindNeighbors(object = ild, dims = 1:20, verbose = F)
ild <- FindClusters(object = ild, resolution = 0.01, verbose = F)

DimPlot(ild)
ild@assays #33694 features for 119438 cells

#cell type annotation, Extended Data Fig10f

FeaturePlot(ild,"EPCAM",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(ild,"PTPRC",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(ild,"PECAM1",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(ild,"MYH11",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(ild,"LUM",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(ild,"COL1A1",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

# Immune cells: PTPRC
# Epithelial cells: EPCAM
# Endothelial cells: PECAM1+ PTPRC-
# Mesenchymal cells: no expression of those markers


Idents(ild, cells = WhichCells(ild, idents = c(0, 3))) <- "Immune"
Idents(ild, cells = WhichCells(ild, idents = c(1, 2))) <- "Epithelial"
Idents(ild, cells = WhichCells(ild, idents = c(4))) <- "Endothelial"
Idents(ild, cells = WhichCells(ild, idents = c(5))) <- "Mesenchymal"

# Extended dATA Fig.10e
DimPlot (ild)

#subset mesenchymal cells including fibroblasts and smooth muscle cells
meso <- subset(ild, idents = "Mesenchymal")


# determine PCA
meso <- FindVariableFeatures(meso, verbose = T, nfeatures = 3000)
meso <- ScaleData(meso,features = row.names(meso@assays$SCT@data))
meso <- RunPCA(meso)

# run UMAP
meso <- RunUMAP(object = meso, dims = 1:14, verbose = F)
meso <- FindNeighbors(meso, dims = 1:20)
meso <- FindClusters(meso, resolution = 0.5)

DimPlot(object = meso, label = T, reduction = "umap", label.size = 5) + theme(aspect.ratio = 0.8) + NoLegend()

# check the duplicate cells and smooth muscle cells
VlnPlot(meso, features = c("PLIN2","HAS1","ACTA2","COL1A1","LUM","WT1","PTPRC","PECAM1","EPCAM"))

# cell type annotation
Idents(meso, cells = WhichCells(meso, idents = c(2, 3))) <- "Smooth Muscle Cells"
Idents(meso, cells = WhichCells(meso, idents = c(6))) <- "Mesothelial Cells"
Idents(meso, cells = WhichCells(meso, idents = c(8,10,12,13,14))) <- "Doublets"
Idents(meso, cells = WhichCells(meso, idents = c(7))) <- "HAS1 High Fibrboblasts"
Idents(meso, cells = WhichCells(meso, idents = c(5))) <- "Fibrboblasts"
Idents(meso, cells = WhichCells(meso, idents = c(0))) <- "PLIN2+ Fibroblasts"
Idents(meso, cells = WhichCells(meso, idents = c(1, 4, 9, 11,15))) <- "Myofibroblasts"

# subset fibroblast populations
FIB <- subset(meso, idents = c("HAS1 High Fibrboblasts","PLIN2+ Fibroblasts","Myofibroblasts","Fibrboblasts"))
FIB@assays #33694 features for 3590 cells

#
FIB1 <- RunUMAP(object = FIB, dims = 1:10, verbose = F)
FIB1 <- FindNeighbors(FIB1, dims = 1:10)
FIB1 <- FindClusters(FIB1, resolution = 0.5)

# Extended Fig.10g
pdf("Extended Fig.10g.pdf", useDingbats = F, height = 3, width = 5)
DimPlot(object = FIB1, label = T, reduction = "umap", label.size = 5) + theme(aspect.ratio = 0.8) + NoLegend()
dev.off()

FeaturePlot(FIB1,"PLIN2",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(FIB1,"PDGFRA",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(FIB1,"ACTA2",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(FIB1,"HAS1",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(FIB1,"LUM",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(FIB1,"COL1A1",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(FIB1,"ETS1",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(FIB1,"SPI1",pt.size = 1.5) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

#Extended Data Fig.10h
FeaturePlot(FIB1, features = c("ETS1", "SPI1"), pt.size = 0.8, blend = T, max.cutoff = 1,cols = c("black", "red", "blue"))
#Extended Data Fig.10j
FeatureScatter(object = FIB1,feature1 = 'ETS1', feature2 = 'SPI1',plot.cor = TRUE, shuffle =T, jitter = TRUE, pt.size = 0.5)

# Reanalysis of scRNA-seq data GSE140023 Conway, B. R. et al.  “Kidney Single-Cell Atlas Reveals Myeloid Heterogeneity in Progression and Regression of Kidney Disease”. J Am Soc Nephrol 2020 Dec;31(12):2833-2854. 
# load data
sham <- Read10X(data.dir = "/Users/Yan/Desktop/")
sham <- CreateSeuratObject(counts = sham, project = "sham", min.cells = 3, min.features = 200)
uuo2 <- Read10X(data.dir = "/Users/Yan/Desktop/")
uuo2 <- CreateSeuratObject(counts = uuo2, project = "uuo2", min.cells = 3, min.features = 200)
uuo7 <- Read10X(data.dir = "/Users/Yan/Desktop/")
uuo7 <- CreateSeuratObject(counts = uuo7, project = "uuo7", min.cells = 3, min.features = 200)
ruuo <- Read10X(data.dir = "/Users/Yan/Desktop/")
ruuo <- CreateSeuratObject(counts = ruuo, project = "ruuo", min.cells = 3, min.features = 200)

View(sham@meta.data)
View(uuo2@meta.data)
View(uuo7@meta.data)
View(ruuo@meta.data)

sham$group <- "sham"
uuo2$group <- "uuo2"
uuo7$group <- "uuo7"
ruuo$group <- "ruuo"

scRNAlist <- list()
scRNAlist[[1]] = NormalizeData(sham)
scRNAlist[[1]] = FindVariableFeatures(scRNAlist[[1]],selection.method = "vst")
scRNAlist[[2]] = NormalizeData(uuo2)
scRNAlist[[2]] = FindVariableFeatures(scRNAlist[[2]],selection.method = "vst")
scRNAlist[[3]] = NormalizeData(uuo7)
scRNAlist[[3]] = FindVariableFeatures(scRNAlist[[3]],selection.method = "vst")
scRNAlist[[4]] = NormalizeData(ruuo)
scRNAlist[[4]] = FindVariableFeatures(scRNAlist[[4]],selection.method = "vst")


scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = 2000)
scRNA <- IntegrateData(anchorset = scRNA.anchors)

UUO <- scRNA
UUO@active.ident
proj_name <- data.frame(proj_name=rep("UUO",ncol(UUO)))

rownames(proj_name) <- row.names(UUO@meta.data)
UUO <- AddMetaData(UUO,proj_name)

DefaultAssay(UUO) <- "RNA"

VlnPlot(UUO, features = c("nFeature_RNA", "nCount_RNA","percent.mt"))
UUO <- subset(UUO, subset = nCount_RNA < 6000)
UUO <- FindVariableFeatures(UUO, selection.method = "vst", nfeatures = 2000)
UUO <- head(VariableFeatures(UUO), 10)
all.genes = rownames(UUO)
UUO = ScaleData(UUO, features = all.genes)
UUO<- RunPCA(UUO, features = VariableFeatures(object = UUO))
UUO <- JackStraw(UUO, num.replicate = 100)
UUO  <- ScoreJackStraw(UUO ,dims = 1:20)
JackStrawPlot(UUO ,dims = 1:20)
ElbowPlot(UUO)

UUO1  <- FindNeighbors(UUO,dims = 1:20)
UUO1  <- FindClusters(UUO1 , resolution = 1.5)
UUO1  <- RunUMAP(UUO1 , dims = 1:20)
DefaultAssay(object = UUO1) <- "RNA"

# Extended Data Fig.10k
DimPlot(UUO1 , reduction = "umap", label = T,pt.size = 2.5)

# cell type annotation, Extended Data Fig.10l

pdf("Kid.prprc.pdf", useDingbats = F, height = 2.5, width = 10)
VlnPlot(UUO1, features = "Ptprc",pt.size = 0) 
dev.off()

pdf("Kid.slc34a1.pdf", useDingbats = F, height = 2.5, width = 10)
VlnPlot(UUO1, features = "Slc34a1",pt.size = 0) 
dev.off()

pdf("Kid.Tfcp2l1.pdf", useDingbats = F, height = 2.5, width = 10)
VlnPlot(UUO1, features = "Tfcp2l1",pt.size = 0) 
dev.off()

pdf("Kid.Emcn.pdf", useDingbats = F, height = 2.5, width = 10)
VlnPlot(UUO1, features = "Emcn",pt.size = 0) 
dev.off()

pdf("Kid.Pdgfra.pdf", useDingbats = F, height = 2.5, width = 10)
VlnPlot(UUO1, features = "Pdgfra",pt.size = 0) 
dev.off()

pdf("Kid.Pdgfrb.pdf", useDingbats = F, height = 2.5, width = 10)
VlnPlot(UUO1, features = "Pdgfrb",pt.size = 0) 
dev.off()

pdf("Kid.Col1a1.pdf", useDingbats = F, height = 2.5, width = 10)
VlnPlot(UUO1, features = "Col1a1",pt.size = 0) 
dev.off()

# subset cluster 24 fibroblasts
FUUO <- subset(UUO1, idents = c("24"))

# check duplicated cells and contaminated cells
DefaultAssay(FUUO) <- "integrated"
FUUO <- ScaleData(FUUO, verbose = FALSE)
FUUO <- RunPCA(FUUO, npcs = 30, verbose = FALSE)
FUUO <- RunUMAP(object = FUUO, dims = 1:20, verbose = F)
FUUO <- FindNeighbors(FUUO, dims = 1:20)
FUUO <- FindClusters(FUUO, resolution = 0.5)
DimPlot(FUUO)
VlnPlot(FUUO, features = c("Ptprc","Epcam","Cdh5","Pdgfra","Pdgfrb")) 

# cluster 3 > Epcam+ epithelial; cluster 4 Cdh5+ endo
# deplete contaminated cells retain fibroblasts
Fibkd <- subset(FUUO, idents = c("0","1","2"))

Fibkd@assays # 201 cells for 17934 features
Fibkd <- ScaleData(Fibkd, verbose = FALSE)
Fibkd <- RunPCA(Fibkd, npcs = 30, verbose = FALSE)

Fibkd <- RunUMAP(object = Fibkd, dims = 1:20, verbose = F)
Fibkd <- FindNeighbors(Fibkd, dims = 1:20)
Fibkd <- FindClusters(Fibkd, resolution = 0.5)

pdf("Extended Data Fig.10m.pdf", useDingbats = F, height = 3, width = 5)
DimPlot(Fibkd, pt.size = 3, cols = c('0'='#ABA300','1'='#31C53F','2'='#00A9FF'))
dev.off()

FeaturePlot(Fibkd,"Col3a1",pt.size = 1.5, min.cutoff = 0) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Fibkd,"Acta2",pt.size = 1.5, min.cutoff = 0) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Fibkd,"Postn",pt.size = 1.5, min.cutoff = 0) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Fibkd,"Sfrp2",pt.size = 1.5, min.cutoff = 0) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

# Extended Data Fig.10n
FeaturePlot(Fibkd, features = c("Ets1", "Spi1"), pt.size = 1.2, blend = T, max.cutoff = 1,cols = c("black", "red", "blue"))
# Extended Data Fig.10o
FeatureScatter(object = Fibkd,feature1 = 'Ets1', feature2 = 'Spi1',plot.cor = TRUE, shuffle =T, jitter = TRUE, pt.size = 2.5)
