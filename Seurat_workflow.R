library(Seurat)
library(ggplot2)
library(patchwork)

HBEC_10x <- Read10X(data.dir="path/to/10Xdata/filtered_feature_bc_matrix")
data <- CreateSeuratObject(counts=HBEC_10x, min.cells=3, min.features=200, project="HBEC")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")

VlnPlot(data, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
p1 <- VlnPlot(data, features="nFeature_RNA") + geom_hline(aes(yintercept=3500), color="blue", linetype='longdash') + geom_hline(aes(yintercept=5000), color="blue", linetype='longdash') + guides(fill="none", color="none", linetype="none")
p2 <- VlnPlot(data, features="nCount_RNA") + geom_hline(aes(yintercept=15000), color="blue", linetype='longdash') + geom_hline(aes(yintercept=32000), color="blue", linetype='longdash') + guides(fill="none", color="none", linetype="none")
p3 <- VlnPlot(data, features="percent.mt") + geom_hline(aes(yintercept=5), color="blue", linetype='longdash') + guides(fill="none", color="none", linetype="none")
pdf("HBEC_scRNA_QC.pdf")
wrap_plots(p1, p2, p3, ncol = 3) 
dev.off()

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

data_filt <- subset(data, subset=nFeature_RNA > 3500 & nFeature_RNA < 5000 & nCount_RNA > 15000 & nCount_RNA < 32000 & percent.mt < 5)
data_filt <- NormalizeData(data_filt, normalization.method="LogNormalize", scale.factor=10000)
data_filt <- NormalizeData(data_filt)

data_filt_df <- as.data.frame(data_filt@assays$RNA$counts)

write.csv(data_filt_df, file="HBEC_10X_filtered.csv", quote=F)

data_filt <- FindVariableFeatures(data_filt, selection.method="vst", nfeatures = 2000)
HVGs <- VariableFeatures(data_filt)
data_filt_df_HVG <- data_filt_df[HVGs,]

write.table(as.data.frame(HVGs), file="HBEC_10X_HVG.txt", quote=F, row.names=F, col.names=F)

