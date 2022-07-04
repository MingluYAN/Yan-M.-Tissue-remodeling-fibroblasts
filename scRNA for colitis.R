# Reanalysis of scRNA-seq data GSE172261 “Colon stroma mediates an inflammation-driven fibroblastic response controlling matrix remodeling and healing”. PLoS Biol 2022 Jan;20(1):e3001532.
# load data
DSS <- Read10X(data.dir = "/Users/Yan/Desktop/", gene.column = 1)
DSS <- CreateSeuratObject(DSS, project = "scRNA DSS", min.features = 300, min.cells = 3)
#QC
DSS <- PercentageFeatureSet(object = DSS, pattern = "^mt-", col.name = "percent.mt")
pdf("DSSQC.pdf", useDingbats = F, height = 15, width = 50)
VlnPlot(DSS, features = c("nFeature_RNA","v", "percent.mt"), ncol = 3)
dev.off()

# Filter out genes with less than 300 features and more than nCount RNA less than 500
DSS <- subset(DSS, subset =  nFeature_RNA > 300 & nCount_RNA > 500)
# Run SCTansform
DSS <- SCTransform(object = DSS, verbose = T)
# Chose a PC for analysis
DSS <- FindVariableFeatures(DSS, verbose = T, nfeatures = 2000)
DSS <- ScaleData(DSS, features = row.names(DSS@assays$SCT@data))
DSS <- RunPCA(DSS)

# RunUMAP, FindClusters, FindNeighbors dims =20
DSS <- RunUMAP(object = DSS, dims = 1:20, verbose = F)
DSS <- FindNeighbors(object = DSS, dims = 1:20, verbose = F)
DSS <- FindClusters(object = DSS, resolution = 1, verbose = F)

pdf("Extended Data Fig.8a.pdf")
DimPlot(DSS, label = T)
dev.off()
        
saveRDS(DSS,file = "220410DSS.RDS") 
DSS <- readRDS("220410DSS.RDS")
DSS@assays #16878 features for 34197 cells
        
Subtype <- c("Pecam1","Kdr","Plvap","Acta2","Cnn1","Mylk","Lyve1","Prox1","Rgs5","Gfap","Fabp7","S100b","Kit","Ano1","Smoc2","Dcn","Fn1","Col1a1","Col3a1","Loxl1","Loxl2","Mmp2")
my <- c("6","19","23","4","13","18","20","24","22","21","25","17","0","1","2","3","5","7","8","9","10","11","12","14","15","16")
factor(Idents(DSS), levels= my)
Idents(DSS) <- factor(Idents(DSS), levels= my)
DotPlot(DSS, features = Subtype, dot.scale = 10) + RotatedAxis() # Extended Data Fig.8b
        
# subset fibroblasts
FibDSS <- subset(DSS, idents = c("0","1","2","3","5","7","8","9","10","11","12","14","15","16"))
FibDSS@assays #24873 cellls
        
# analysis
FibDSS <- RunUMAP(object = FibDSS, dims = 1:10, verbose = F)
FibDSS <- FindNeighbors(FibDSS, dims = 1:10)
FibDSS <- FindClusters(FibDSS, resolution = 0.5) 
DimPlot(FibDSS,label = F)
        
# check duplicaed cells
VlnPlot(FibDSS, features = c("Pecam1","Epcam","Rgs5","Plvap")) 
        
# remove cluster 12 Rgs5+ perivascular cells
FibDSS2 <- subset(FibDSS, idents = c("0","1","2","3","4","5","6","7","8","9","10","11"))
FibDSS2@assays #24566 fibroblasts
        
#analysis fibroblasts
FibDSS2 <- RunUMAP(object = FibDSS2, dims = 1:10, verbose = F)
FibDSS2 <- FindNeighbors(FibDSS2, dims = 1:10)
FibDSS2 <- FindClusters(FibDSS2, resolution = 0.15) 
        
pdf("Fig.7a.pdf", useDingbats = F, height = 5, width = 7)
DimPlot(FibDSS2, reduction = "umap",pt.size = 1, cols = c('0'='#ABA300','1'='#31C53F','2'='#00A9FF','3'='#C77CFF','4'= "orange"))
dev.off()
        
Sig <- c("Ets1","Wnt5a","Bmp2","Bmp7","Bmp5","Pdgfra","Acta2","Procr","F3","Adamdec1","Il6","Cxcl1","Il33","Grem1","Thy1","Rspo1","Wnt2","Wnt2b","Il11")
mylev <- c("2","3","0","4","1")
factor(Idents(FibDSS2), levels= mylev)
Idents(FibDSS2) <- factor(Idents(FibDSS2), levels= mylev)
# Fig.7b
DotPlot(FibDSS2, features = Sig, dot.scale = 10) + RotatedAxis()