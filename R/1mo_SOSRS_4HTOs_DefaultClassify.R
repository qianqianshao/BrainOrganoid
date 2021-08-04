# 12.30.2020 by Qianyi
R 

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)

home="/scratch/junzli_root/junzli/qzm/Dropseq_analysis/"
source(paste0(home,"data_DGE/Rcode_multiplot.R"))

dgefile="/gpfs/accounts/junzli_root/junzli/qzm/Dropseq_analysis/HTO/"
setwd(dgefile)
names=c("3GEX","HTO")
names=c("3GEX","HTO")
run=c("NovaA-320","NovaA-320")
sample=c("2025-AT","2025-AT")
dataset=c("2025-AT-1-3GEX_CGTGACAT-TTTAGACC","2025-AT-1-HTO")
subject=c("Human1","Human1")
organ=c("Brain","Brain")
n=length(dataset)
datainfo=data.frame(run,dataset,names,subject,organ)


# Read in data

# Load in the UMI matrix
i=1
#	dir=paste0("/scratch/junzli_root/junzli/qzm/Data/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[i],"/outs/filtered_feature_bc_matrix/")
  dir=paste0("/nfs/turbo/umms-junzli/qzm/Data/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[i],"/outs/filtered_feature_bc_matrix/")
	data <- Read10X(data.dir = dir)
dim(data) # [1] 36601  8964
	dgedata=as.matrix(data)
	### Filter for cells (>500 genes and <10% MT)
nGeneperCell <- colSums(dgedata>0)
nUMIperCell <- colSums(dgedata)
print(c(mean(nUMIperCell),mean(nGeneperCell)))
dgedata.tmp=dgedata[,nGeneperCell>500]
mito.genes <- grep("^MT-|^mt-", rownames(dgedata.tmp), value = T) 
percent.mito <- colSums(dgedata.tmp[mito.genes, ])/colSums(dgedata.tmp)
dgedata.tmp=dgedata[,nGeneperCell>500][,percent.mito<0.1] 
print(c(ncol(dgedata),ncol(dgedata[,nGeneperCell>500]),ncol(dgedata.tmp)))
print(c(ncol(dgedata),ncol(dgedata[,nGeneperCell>500]),ncol(dgedata.tmp)))
### Filter for genes (Non-0 genes)
nCellperGene <- rowSums(dgedata.tmp>0)
dgedata2=dgedata.tmp[which(nCellperGene >= 1),]
print(c(nrow(dgedata),nrow(dgedata2)))
print(summary(rowSums(dgedata2)))

# Load hashtag count matrix 
i=2
	dir2=paste0("2025-AT-1-HTO2/umi_count/")
	data2 <- Read10X(data.dir = dir2,gene.column=1)
dim(data2) #[1]    5 19980
data2[,1:2]
#                     AAGAACAAGTGATAAC CGATGCGAGCATGATA
#tag1-GTCAACTCTTTAGCG               26               26
#tag2-AGTAAGTTCAGCGTA               18              191
#tag3-TTCCGCCTCTCTTTG              270              139
#tag4-TGATGGCCTATTGGG               74              141
#unmapped                           23               23
data2=as.matrix(data2)
data2=data2[1:4,]
rownames(data2)=paste0("HTO-",c("A","B","C","D"))
data2[,1:2]
#      AAGAACAAGTGATAAC CGATGCGAGCATGATA
#HTO-A               26               26
#HTO-B               18              191
#HTO-C              270              139
#HTO-D               74              141

nHTOperCell <- colSums(data2)
nHTOperCell <- sort(nHTOperCell,decreasing=T)
nHTOperCell2 = nHTOperCell
cells2=names(nHTOperCell2)
max(nHTOperCell2[which(!(cells2 %in% cells1))]) # 670
png(file="plot/hto_counts_20kcell.png",res=300,height=1000,width=2000)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(1:length(nHTOperCell),nHTOperCell,log="xy",pch=16,cex=0.8,col=rgb(0,0,0,0.5),xlab="Cell Barcode",ylab="Total HTO counts")
plot(1:length(nHTOperCell),nHTOperCell,log="y",pch=16,cex=0.8,col=rgb(0,0,0,0.5),xlab="Cell Barcode",ylab="Total HTO counts")
dev.off()

nHTOperCell <- sort(nHTOperCell,decreasing=T)
cells=names(nHTOperCell)
nAperCell <- data2[1,cells]
nBperCell <- data2[2,cells]
nCperCell <- data2[3,cells]
nDperCell <- data2[4,cells]
which(names(nAperCell) != cells)
png(file="plot/hto_4counts_20kcell.png",res=300,height=1000,width=1000)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(1:length(nHTOperCell),nDperCell,log="xy",pch=16,cex=0.5,col=rgb(1,0,0,0.8),xlab="Cell Barcode",ylab="4 HTO counts")
points(1:length(nHTOperCell),nCperCell,pch=16,cex=0.5,col=rgb(0,1,0,0.6))
points(1:length(nHTOperCell),nBperCell,pch=16,cex=0.5,col=rgb(0,0,1,0.4))
points(1:length(nHTOperCell),nAperCell,pch=16,cex=0.5,col=rgb(1,0,1,0.3))
dev.off()


# Select cell barcodes detected by both RNA and HTO 
colnames(data2)=paste0(colnames(data2),"-1")
dge.umis=dgedata2
dge.htos=data2
dim(dge.umis) # [1] 27420  7799
dim(dge.htos) # [1]     4 19980
joint.bcs <- intersect(colnames(dge.umis), colnames(dge.htos))
length(joint.bcs) # 6009 # used this


# Subset RNA and HTO counts by joint cell barcodes
dge.umis <- dge.umis[, joint.bcs]
dge.htos <- as.matrix(dge.htos[, joint.bcs])
write.table(dge.htos,paste0("plot/",sample[2],"_rawHTO.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(dge.umis,paste0("plot/",sample[2],"_rawRNA.txt"),col.names=T,row.names=T,sep="\t",quote=F)
# uploaded to GEO

nHTOperCell <- colSums(dge.htos)
nHTOperCell <- sort(nHTOperCell,decreasing=T)
png(file="plot/hto_counts_6kcell.png",res=300,height=1000,width=2000)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(1:length(nHTOperCell),sort(nHTOperCell,decreasing=T),log="xy",pch=16,cex=0.8,col=rgb(0,0,0,0.5),xlab="Cell Barcode",ylab="Total HTO counts")
plot(1:length(nHTOperCell),sort(nHTOperCell,decreasing=T),log="y",pch=16,cex=0.8,col=rgb(0,0,0,0.5),xlab="Cell Barcode",ylab="Total HTO counts")
dev.off()
cells=names(nHTOperCell)
nAperCell <- dge.htos[1,cells]
nBperCell <- dge.htos[2,cells]
nCperCell <- dge.htos[3,cells]
nDperCell <- dge.htos[4,cells]
which(names(nAperCell) != cells)
png(file="plot/hto_4counts_6kcell.png",res=300,height=1000,width=1000)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(1:length(nHTOperCell),nDperCell,log="xy",pch=16,cex=0.5,col=rgb(1,0,0,0.8),xlab="Cell Barcode",ylab="4 HTO counts")
points(1:length(nHTOperCell),nCperCell,pch=16,cex=0.5,col=rgb(0,1,0,0.6))
points(1:length(nHTOperCell),nBperCell,pch=16,cex=0.5,col=rgb(0,0,1,0.4))
points(1:length(nHTOperCell),nAperCell,pch=16,cex=0.5,col=rgb(1,0,1,0.3))
dev.off()


# Confirm that the HTO have the correct names
rownames(dge.htos)
## [1] "HTO_A" "HTO_B" "HTO_C" "HTO_D" "HTO_E" "HTO_F" "HTO_G" "HTO_H"



### Setup Seurat object and add in the HTO data

# Setup Seurat object
dge.hashtag <- CreateSeuratObject(counts = dge.umis)

# Normalize RNA data with log normalization
dge.hashtag <- NormalizeData(dge.hashtag)
# Find and scale variable features
dge.hashtag <- FindVariableFeatures(dge.hashtag, selection.method = "mean.var.plot")
dge.hashtag <- ScaleData(dge.hashtag, features = VariableFeatures(dge.hashtag))

dge.hashtag[["percent.mt"]] <- PercentageFeatureSet(dge.hashtag, pattern = "^MT-")
print(c(mean(dge.hashtag$nFeature_RNA),mean(dge.hashtag$nCount_RNA),mean(dge.hashtag$nCount_HTO),mean(dge.hashtag$percent.mt)))

### Adding HTO data as an independent assay

# Add HTO data as a new assay independent from RNA
dge.hashtag[["HTO"]] <- CreateAssayObject(counts = dge.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
dge.hashtag <- NormalizeData(dge.hashtag, assay = "HTO", normalization.method = "CLR")

# details for normalization of HTO counts
setwd("D:/HTO_BrainOrganoids/fig_4HTO_202012")
library(pbapply)
library(future)
normalization.method = "CLR"
object=read.table("hto_counts_6kcell.txt")
object = as.matrix(dge.hashtag@assays$HTO@counts)
  scale.factor = 1e4;
  margin = 1;
  block.size = NULL;
  verbose = TRUE;
CustomNormalize <- function(data, custom_function, margin, verbose = TRUE) {
  if (is.data.frame(x = data)) {
    data <- as.matrix(x = data)
  }
  if (!inherits(x = data, what = 'dgCMatrix')) {
    data <- as(object = data, Class = "dgCMatrix")
  }
  myapply <- ifelse(test = verbose, yes = pbapply, no = apply)
  # margin <- switch(
  #   EXPR = across,
  #   'cells' = 2,
  #   'features' = 1,
  #   stop("'across' must be either 'cells' or 'features'")
  # )
  if (verbose) {
    message("Normalizing across ", c('features', 'cells')[margin])
  }
  norm.data <- myapply(
    X = data,
    MARGIN = margin,
    FUN = custom_function)
  if (margin == 1) {
    norm.data = t(x = norm.data)
  }
  colnames(x = norm.data) <- colnames(x = data)
  rownames(x = norm.data) <- rownames(x = data)
  return(norm.data)
}
    bb=switch(
      EXPR = normalization.method,
      'CLR' = CustomNormalize(
        data = object,
        custom_function = function(x) {
          return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
        },
        margin = margin,
        verbose = verbose
        # across = across
      )
    )
x=object
c=log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x))))
a = as.matrix(dge.hashtag@assays$HTO@data)
# first convert the counts to a sparse matrix if stored in other formats
data <- as.matrix(x = object)
data <- as(object = data, Class = "dgCMatrix")
# perform CLR normalization:
norm.data <- myapply(
    X = data,
    MARGIN = 1,
    FUN = function(x) {
          return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
        }) 
b <- t(norm.data)
which(b != bb) # integer(0)
bb[,1:4]
a[,1:4]
bb[,6008:6009]
a[,6008:6009]
write.table(a,"plot/hto_normalized_6kcell.txt",row.names=T,col.names=T,quote=F,sep="\t")

# Demultiplex cells based on HTO enrichment
# Here we use the Seurat function HTODemux() to assign single cells back to their sample origins.
# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using the
# default settings
dge.hashtag <- HTODemux(dge.hashtag, assay = "HTO", positive.quantile = 0.99)
dge.hashtag$ID <- Idents(dge.hashtag)
# HTO 20k cells - overlapped 6k cells
Cutoff for HTO-A : 37 reads
Cutoff for HTO-B : 44 reads
Cutoff for HTO-C : 117 reads
Cutoff for HTO-D : 303 reads
# HTO 30k cells - overlapped 6.8k cells
Cutoff for HTO-A : 40 reads
Cutoff for HTO-B : 49 reads
Cutoff for HTO-C : 148 reads
Cutoff for HTO-D : 293 reads


# check distribution of the HTO tags before deciding the quantile
a = as.matrix(dge.hashtag@assays$HTO@data)
png(file="plot/hto_4counts_6kcell_normalized.png",res=300,height=1000,width=2000)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(1.5, .5, 0))
plot(density(a[1,]),col=1,main="",ylab="Density")
legend("topright",rownames(a)[1],lty=1,col=1,text.col=1)
plot(density(a[2,]),col=2,main="",ylab="Density")
legend("topright",rownames(a)[2],lty=1,col=2,text.col=2)
plot(density(a[3,]),col=3,main="",ylab="Density")
legend("topright",rownames(a)[3],lty=1,col=3,text.col=3)
plot(density(a[4,]),col=4,main="",ylab="Density")
legend("topright",rownames(a)[4],lty=1,col=4,text.col=4)
dev.off()



# Visualize demultiplexing results
# Output from running HTODemux() is saved in the object metadata. We can visualize how many cells are classified as singlets, doublets and negative/ambiguous cells.

# Global classification results
table(dge.hashtag$HTO_classification.global)
# HTO 20k cells - overlapped 6k cells
 Doublet Negative  Singlet 
    1903     1054     3052 

class=dge.hashtag@meta.data$HTO_classification
names(class)=rownames(dge.hashtag@meta.data)
write.table(class,"plot/hto_defaultclassify_6kcell.txt",row.names=T,col.names=F,quote=F,sep="\t")


# Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
Idents(dge.hashtag) <- "ID"
p=list()
p[[1]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")
p[[2]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-A", feature2 = "hto_HTO-C")
p[[3]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-A", feature2 = "hto_HTO-D")
p[[4]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-B", feature2 = "hto_HTO-C")
p[[5]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-B", feature2 = "hto_HTO-D")
p[[6]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-C", feature2 = "hto_HTO-D")
png(file="plot/HTO_pairs.png",res=300,height=2000,width=3600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(p,cols=3)
dev.off()

# Visualize enrichment for selected HTOs with ridge plots

# Group cells based on the max HTO signal
Idents(dge.hashtag) <- "HTO_maxID"
pdf(file="plot/HTO_maxID.pdf",height=6,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
RidgePlot(dge.hashtag, assay = "HTO", features = rownames(dge.hashtag[["HTO"]]), ncol = 2)
dev.off()


# Compare number of UMIs for singlets, doublets and negative cells

Idents(dge.hashtag) <- "HTO_classification.global"
pdf(file="plot/HTO_classification.pdf",height=3,width=4)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
VlnPlot(dge.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()

# Generate a two dimensional tSNE embedding for HTOs.Here we are grouping cells by singlets and doublets for simplicity.

# First, we will remove negative cells from the object
dge.hashtag.subset <- subset(dge.hashtag, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(dge.hashtag.subset) <- "HTO"
dge.hashtag.subset <- ScaleData(dge.hashtag.subset, features = rownames(dge.hashtag.subset), 
    verbose = FALSE)
dge.hashtag.subset <- RunPCA(dge.hashtag.subset, features = rownames(dge.hashtag.subset), approx = FALSE)
dge.hashtag.subset <- RunTSNE(dge.hashtag.subset, dims = 1:4, perplexity = 100)
png(file="plot/HTO_tSNE.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.hashtag.subset)
dev.off()

# You can also visualize the more detailed classification result by running Idents(object) <-
# 'HTO_classification' before plotting. Here, you can see that each of the small clouds on the
# tSNE plot corresponds to one of the 28 possible doublet combinations.
# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
pdf(file="plot/HTO_classification_heatmap.pdf",height=4,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
HTOHeatmap(dge.hashtag, assay = "HTO", ncells = 5000)
dev.off()

# Cluster and visualize cells using the usual scRNA-seq workflow, and examine for the potential presence of batch effects.

# Extract the singlets
dge.singlet <- subset(dge.hashtag, idents = "Singlet")

# Select the top 1000 most variable features
dge.singlet <- FindVariableFeatures(dge.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
dge.singlet <- ScaleData(dge.singlet, features = VariableFeatures(dge.singlet))

# Run PCA
dge.singlet <- RunPCA(dge.singlet, features = VariableFeatures(dge.singlet))
# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
dge.singlet <- FindNeighbors(dge.singlet, reduction = "pca", dims = 1:10)
dge.singlet <- FindClusters(dge.singlet, resolution = 0.6, verbose = FALSE)
dge.singlet <- RunTSNE(dge.singlet, reduction = "pca", dims = 1:10)

# Projecting singlet identities on TSNE visualization
png(file="plot/HTO_singlet_tSNE.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, group.by = "HTO_classification")
dev.off()


save(dge.hashtag,file="HTO20kcell_6kcell_hashtag.Robj")
save(dge.singlet,file="HTO20kcell_6kcell_singlet.Robj")

