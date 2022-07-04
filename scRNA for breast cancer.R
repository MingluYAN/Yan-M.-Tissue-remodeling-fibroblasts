# Reanalysis of scRNA-seq data GSE176078 Wu, S.Z. et al. “A single-cell and spatially resolved atlas of human breast cancers”. Nat Genet 2021 Sep;53(9):1334-1347.

# load data
Bcar <- Read10X(data.dir = "/Users/Yan/Desktop/", gene.column = 1)
Bcar <- CreateSeuratObject(Bcar, project = "breast cancer", min.features = 200)

# QC
Bcar <- PercentageFeatureSet(object = Bcar, pattern = "^MT-", col.name = "percent.mt")
pdf("Bcar.QC.pdf", useDingbats = F, height = 15, width = 50)
VlnPlot(Bcar, features = c("nFeatureRNA","nCountRNA", "percent.mt"), ncol = 3)
dev.off()

# Run SCTansform
Bcar <- SCTransform(object = Bcar, verbose = T)

# Chose a PC for analysis
Bcar<- FindVariableFeatures(Bcar, verbose = T, nfeatures = 3000)
Bcar <- ScaleData(Bcar, features = row.names(Bcar@assays$SCT@data))
Bcar <- RunPCA(Bcar)

ElbowPlot(Bcar)

# RunUMAP, FindClusters, FindNeighbors dims = 1: 20 

Bcar <- RunUMAP(object = Bcar, dims = 1:20, verbose = F)
Bcar <- FindNeighbors(object = Bcar, dims = 1:20, verbose = F)
Bcar <- FindClusters(object = Bcar, resolution = 1.5, verbose = F)

# Fig.6c
DimPlot(Bcar, reduction = "umap", label = T,pt.size = 0.5)

# cell type annotation, Fig.6d
FeaturePlot(Bcar,"EPCAM",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Bcar,"PTPRC",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Bcar,"MCAM",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Bcar,"RGS5",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Bcar,"PFN2",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Bcar,"PECAM1",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Bcar,"PDGFRA",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Bcar,"PDGFRB",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Bcar,"COL1A1",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(Bcar,"LRRC15",pt.size = 1) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

# subset CAFs populations
CAFs <- subset(Bcar, idents = c("48","44","16","27","42")) # 6362 CAFs
CAFs @assay

# myCAF and iCAF signatures
total <- c("TNFSF11","ETS1","LRRC15","INHBA","ACTA2","TPM1","TAGLN","CTHRC1","CLEC3B","C7","IL33","CXCL2","C3")
mylev <- c("16","44","42","27","48")
factor(Idents(CAFs), levels= mylev)
Idents(CAFs) <- factor(Idents(CAFs), levels= mylev)
#Fig.6e
DotPlot(CAFs, features = total, dot.scale = 30) + RotatedAxis() 
