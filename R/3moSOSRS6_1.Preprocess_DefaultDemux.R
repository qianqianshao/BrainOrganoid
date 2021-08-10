# scRNA-seq analysis for 6 SOSRSs at 3-month old with HTOs on 2.4.2021 by Qianyi
# 6 replicates sequenced in 2 samples, each with 3 SOSRS (3 HTO tags)
# kept cell barcodes present in both raw RNA counts and HTO counts matrices
# used default HTO Demultiplexing from Seurat3

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
run=c("NovaZ-9","NovaZ-9")
sample=c("2514-AT","2514-AT")
gexset=c("2514-AT-1","2514-AT-2")
htoset=list(
  1:3,
  4:6
)
dataset=list(
  c("2514-AT-1-3GEX_GCTACAAA-AGGGCACG","2514-AT-1-HTO"),
  c("2514-AT-2-3GEX_CGCTGAAA-GCAGACAC","2514-AT-2-HTO")
)
subject=c("Human2","Human3")
organ=c("Brain","Brain")
n=length(dataset)
datainfo=data.frame(run,dataset,names,subject,organ)


gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}


# Read in data

# Load in the UMI matrix
i=1
datalist=list()
for(i in 1:2){
#	dir=paste0("/scratch/junzli_root/junzli/qzm/Data/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[[i]][1],"/outs/filtered_feature_bc_matrix/")
  dir=paste0("/nfs/turbo/umms-junzli/qzm/Data/",run[i],"/Client/",sample[i],"/",run[i],"/Sample_",dataset[[i]][1],"/outs/filtered_feature_bc_matrix/")
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
### Filter for genes (Non-0 genes)
nCellperGene <- rowSums(dgedata.tmp>0)
dgedata2=dgedata.tmp[which(nCellperGene >= 1),]
print(c(nrow(dgedata),nrow(dgedata2)))
print(summary(rowSums(dgedata2)))
colnames(dgedata2)=paste(gexset[i],colnames(dgedata2),sep="_")
datalist[[i]]=dgedata2
}

# Load in the HTO count matrix
i=1
htolist=list()
for(i in 1:2){
	dir1=paste0(dataset[[i]][2],"/umi_count/")
	data1 <- Read10X(data.dir = dir1,gene.column=1)
print(dim(data1)) #[1]    4 200000
data1[,1:2]
data1=as.matrix(data1)
data1=data1[1:3,]
rownames(data1)=paste0("HTO-",htoset[[i]])
colnames(data1)=paste(gexset[i],colnames(data1),sep="_")
colnames(data1)=paste0(colnames(data1),"-1")
htolist[[i]]=data1
print(data1[,1:2])
#      2514-AT-1_GCCTGTTTCAACCCGG 2514-AT-1_TCACTCGAGAGTGTGC
#HTO-1                          0                          2
#HTO-2                          1                          2
#HTO-3                          2                          3
#      2514-AT-2_GGCAGTCTCCCGAATA 2514-AT-2_TCCTCCCCAAGACAAT
#HTO-4                         22                         28
#HTO-5                         37                         26
#HTO-6                         47                          7


nHTOperCell <- colSums(data1)
nHTOperCell <- sort(nHTOperCell,decreasing=T)
nHTOperCell1 = nHTOperCell
cells1=names(nHTOperCell1)
png(file=paste0("plot/hto_",i,"_counts_200kcells.png"),res=300,height=1000,width=2000)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(1:length(nHTOperCell),nHTOperCell,log="xy",pch=16,cex=0.8,col=rgb(0,0,0,0.5),xlab="Cell Barcode",ylab="Total HTO counts")
plot(1:length(nHTOperCell),nHTOperCell,log="y",pch=16,cex=0.8,col=rgb(0,0,0,0.5),xlab="Cell Barcode",ylab="Total HTO counts")
dev.off()

nHTOperCell <- sort(nHTOperCell,decreasing=T)
cells=names(nHTOperCell)
nAperCell <- data1[1,cells]
nBperCell <- data1[2,cells]
nCperCell <- data1[3,cells]
which(names(nAperCell) != cells)
png(file=paste0("plot/hto_",i,"3counts_100kcell.png"),res=300,height=1000,width=1000)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(1:length(nHTOperCell),nCperCell,log="xy",pch=16,cex=0.5,col=rgb(1,0,0,0.8),xlab="Cell Barcode",ylab="HTO counts")
points(1:length(nHTOperCell),nBperCell,pch=16,cex=0.5,col=rgb(0,0,1,0.4))
points(1:length(nHTOperCell),nAperCell,pch=16,cex=0.5,col=rgb(0,1,0,0.3))
dev.off()
}

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(dge.umis), colnames(dge.htos))
length(joint.bcs) # 16916
i=1
umis=list()
htos=list()
for(i in 1:2){
dge.umis=datalist[[i]]
dge.htos=htolist[[i]]
print(dim(dge.umis)) # [1] 28049  6374
print(dim(dge.htos)) # [1]     3 200000
joint.bcs <- intersect(colnames(dge.umis), colnames(dge.htos))
print(length(joint.bcs)) # 6373


# Subset RNA and HTO counts by joint cell barcodes
dge.umis <- dge.umis[, joint.bcs]
dge.htos <- as.matrix(dge.htos[, joint.bcs])
write.table(dge.htos,paste0("plot/hto_",i,"_counts_6kcell.txt"),row.names=T,col.names=T,quote=F,sep="\t")
umis[[i]]=dge.umis
htos[[i]]=dge.htos

if(i==1){
  rownames(dge.htos)=c("HTO-A","HTO-B","HTO-C")
} 
if(i==2){
  rownames(dge.htos)=c("HTO-D","HTO-E","HTO-F")
}
write.table(dge.htos,paste0("plot/",gexset[i],"_rawHTO.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(dge.umis,paste0("plot/",gexset[i],"_rawRNA.txt"),col.names=T,row.names=T,sep="\t",quote=F)



nHTOperCell <- colSums(dge.htos)
nHTOperCell <- sort(nHTOperCell,decreasing=T)
png(paste0(file="plot/hto_",i,"_counts_6kcell.png"),res=300,height=1000,width=2000)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(1:length(nHTOperCell),sort(nHTOperCell,decreasing=T),log="xy",pch=16,cex=0.8,col=rgb(0,0,0,0.5),xlab="Cell Barcode",ylab="Total HTO counts")
plot(1:length(nHTOperCell),sort(nHTOperCell,decreasing=T),log="y",pch=16,cex=0.8,col=rgb(0,0,0,0.5),xlab="Cell Barcode",ylab="Total HTO counts")
dev.off()
cells=names(nHTOperCell)
nAperCell <- dge.htos[1,cells]
nBperCell <- dge.htos[2,cells]
nCperCell <- dge.htos[3,cells]
which(names(nAperCell) != cells)
png(paste0(file="plot/hto_",i,"_3counts_6kcell.png"),res=300,height=1000,width=1000)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(1:length(nHTOperCell),nCperCell,log="xy",pch=16,cex=0.5,col=rgb(1,0,0,0.8),xlab="Cell Barcode",ylab="HTO counts")
points(1:length(nHTOperCell),nBperCell,pch=16,cex=0.5,col=rgb(0,0,1,0.4))
points(1:length(nHTOperCell),nAperCell,pch=16,cex=0.5,col=rgb(0,1,0,0.3))
dev.off()

}


# Confirm that the HTO have the correct names
rownames(dge.htos)

### Setup Seurat object and add in the HTO data
i=1
dgelist=list()
for(i in 1:2){
dge.umis=umis[[i]]
dge.htos=htos[[i]]

# Setup Seurat object
dge.hashtag <- CreateSeuratObject(counts = dge.umis)
  dge.hashtag[["percent.mt"]] <- PercentageFeatureSet(dge.hashtag, pattern = "^MT-")
print(c(mean(dge.hashtag$nFeature_RNA),mean(dge.hashtag$nCount_RNA),mean(dge.hashtag$nCount_HTO),mean(dge.hashtag$percent.mt)))

# Normalize RNA data with log normalization
dge.hashtag <- NormalizeData(dge.hashtag)
# Find and scale variable features
dge.hashtag <- FindVariableFeatures(dge.hashtag, selection.method = "mean.var.plot")
dge.hashtag <- ScaleData(dge.hashtag, features = VariableFeatures(dge.hashtag))


### Adding HTO data as an independent assay

# Add HTO data as a new assay independent from RNA
dge.hashtag[["HTO"]] <- CreateAssayObject(counts = dge.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
dge.hashtag <- NormalizeData(dge.hashtag, assay = "HTO", normalization.method = "CLR")

dge.hashtag
dgelist[[i]]=dge.hashtag

}

# Demultiplex cells based on HTO enrichment
# Here we use the Seurat function HTODemux() to assign single cells back to their sample origins.
# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using the
# default settings
for(i in 1:2){
dge.hashtag <- dgelist[[i]]
dge.hashtag <- HTODemux(dge.hashtag, assay = "HTO", positive.quantile = 0.99)
dge.hashtag$ID <- Idents(dge.hashtag)
save(dge.hashtag,file=paste0(gexset[i],"_hashtag.Robj"))

# HTO 200k cells - overlapped 6k cells
Cutoff for HTO-1 : 38 reads
Cutoff for HTO-2 : 25 reads
Cutoff for HTO-3 : 50 reads

Cutoff for HTO-4 : 33 reads
Cutoff for HTO-5 : 41 reads
Cutoff for HTO-6 : 22 reads


# check distribution of the HTO tags before deciding the quantile
a = as.matrix(dge.hashtag@assays$HTO@data)
png(file=paste0("plot/hto_",i,"3counts_6kcell_normalized.png"),res=300,height=1000,width=2500)
par(mfrow=c(2,3),mar=c(4,4,1,1),mgp=c(1.5, .5, 0))
plot(density(a[1,])$y,pch=16,cex=0.8,col=1,main="",ylab="Density")
legend("topright",rownames(a)[1],lty=1,col=1,text.col=1)
plot(density(a[2,])$y,pch=16,cex=0.8,col=2,main="",ylab="Density")
legend("topright",rownames(a)[2],lty=1,col=2,text.col=2)
plot(density(a[3,])$y,pch=16,cex=0.8,col=3,main="",ylab="Density")
legend("topright",rownames(a)[3],lty=1,col=3,text.col=3)
plot(density(a[1,]),col=1,main="",ylab="Density")
legend("topright",rownames(a)[1],lty=1,col=1,text.col=1)
plot(density(a[2,]),col=2,main="",ylab="Density")
legend("topright",rownames(a)[2],lty=1,col=2,text.col=2)
plot(density(a[3,]),col=3,main="",ylab="Density")
legend("topright",rownames(a)[3],lty=1,col=3,text.col=3)
dev.off()


hto=t(as.matrix(dge.hashtag@assays$HTO@counts))
a=t(hto)
png(paste0(file="plot/hto_",i,"_3counts_6kcell_density.png"),res=300,height=1000,width=2500)
par(mfrow=c(2,3),mar=c(4,4,1,1),mgp=c(1.5, .5, 0))
plot(density(a[1,])$y,pch=16,cex=0.8,col=1,main="",ylab="Density")
legend("topright",rownames(a)[1],lty=1,col=1,text.col=1)
plot(density(a[2,])$y,pch=16,cex=0.8,col=2,main="",ylab="Density")
legend("topright",rownames(a)[2],lty=1,col=2,text.col=2)
plot(density(a[3,])$y,pch=16,cex=0.8,col=3,main="",ylab="Density")
legend("topright",rownames(a)[3],lty=1,col=3,text.col=3)
plot(density(a[1,]),col=1,main="",ylab="Density")
legend("topright",rownames(a)[1],lty=1,col=1,text.col=1)
plot(density(a[2,]),col=2,main="",ylab="Density")
legend("topright",rownames(a)[2],lty=1,col=2,text.col=2)
plot(density(a[3,]),col=3,main="",ylab="Density")
legend("topright",rownames(a)[3],lty=1,col=3,text.col=3)
dev.off()

# Visualize demultiplexing results
# Output from running HTODemux() is saved in the object metadata. We can visualize how many cells are classified as singlets, doublets and negative/ambiguous cells.

# Global classification results
table(dge.hashtag$HTO_classification.global)
## 
##  Doublet Negative  Singlet 
##     2598      346    13972

class=dge.hashtag@meta.data$HTO_classification
names(class)=rownames(dge.hashtag@meta.data)
write.table(class,"plot/hto_",i,"_defaultclassify_6kcell.txt",row.names=T,col.names=F,quote=F,sep="\t")

# Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
Idents(dge.hashtag) <- "ID"
table(dge.hashtag$ID)
 Doublet    HTO-A    HTO-B    HTO-C    HTO-D Negative 
    1903      538      806      805      903     1054

p=list()
p[[1]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = paste0("hto_HTO-",htoset[[i]][1]), feature2 = paste0("hto_HTO-",htoset[[i]][2]))
p[[2]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = paste0("hto_HTO-",htoset[[i]][1]), feature2 = paste0("hto_HTO-",htoset[[i]][3]))
p[[3]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = paste0("hto_HTO-",htoset[[i]][2]), feature2 = paste0("hto_HTO-",htoset[[i]][3]))
png(file=paste0("plot/hto_",i,"_pairs.png"),res=300,height=1000,width=3600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(p,cols=3)
dev.off()

# zoomed-in view, order the singlets to be on top
Idents(dge.hashtag) <- "ID"
ids=as.numeric(dge.hashtag$ID)
names(ids)=names(dge.hashtag$ID)
table(dge.hashtag$HTO_classification.global)
order=names(sort(dge.hashtag$HTO_classification.global))

feature1=paste0("hto_HTO-",htoset[[i]][1])
feature2=paste0("hto_HTO-",htoset[[i]][2])
feature3=paste0("hto_HTO-",htoset[[i]][3])
data = FetchData(object = dge.hashtag, vars = c(feature1, feature2, feature3))
data=data[order,]
ids=ids[order]
ids1=ids2=ids3=ids
ids1[which(ids>3)]<-1
ids2[which(ids==3 | ids==5)]<-1
ids2[which(ids==4)]<-3
ids3[which(ids==2 | ids==5)]<-1
ids3[which(ids==3)]<-2
ids3[which(ids==4)]<-3
png(file=paste0("plot/hto_",i,"FitClassify_pairs_zoom2.png"),res=300,height=1000,width=3000)
par(mfrow=c(1,3),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(data[,1],data[,2],cex.lab=1.5,xlim=c(0,2),ylim=c(0,2),xlab=colnames(data)[1],ylab=colnames(data)[2],pch=16,col=ids1,cex=.5)
plot(data[,1],data[,3],cex.lab=1.5,xlim=c(0,2),ylim=c(0,2),xlab=colnames(data)[1],ylab=colnames(data)[3],pch=16,col=ids2,cex=.5)
plot(data[,2],data[,3],cex.lab=1.5,xlim=c(0,2),ylim=c(0,2),xlab=colnames(data)[2],ylab=colnames(data)[3],pch=16,col=ids3,cex=.5)
dev.off()

 Visualize enrichment for selected HTOs with ridge plots

# Group cells based on the max HTO signal
Idents(dge.hashtag) <- "HTO_maxID"
pdf(file=paste0("plot/hto_",i,"_maxID.pdf"),height=3,width=9)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
RidgePlot(dge.hashtag, assay = "HTO", features = rownames(dge.hashtag[["HTO"]]), ncol = 3)
dev.off()


# Compare number of UMIs for singlets, doublets and negative cells

Idents(dge.hashtag) <- "HTO_classification.global"
pdf(file=paste0("plot/hto_",i,"_classification.pdf"),height=4,width=5)
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
dge.hashtag.subset <- RunTSNE(dge.hashtag.subset, dims = 1:3, perplexity = 100,check_duplicates = FALSE)

# remove negative and visulize classification using default classification
png(file=paste0("plot/hto_",i,"Classify_PCA2.png"),res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.hashtag.subset,reduction="pca",group="hash.ID")
dev.off()
png(file=paste0("plot/hto_",i,"Classify_tSNE.png"),res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.hashtag.subset)
dev.off()
png(file=paste0("plot/hto_",i,"Classify_tSNE2.png"),res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.hashtag.subset,group="hash.ID")
dev.off()

# You can also visualize the more detailed classification result by running Idents(object) <-
# 'HTO_classification' before plotting. Here, you can see that each of the small clouds on the
# tSNE plot corresponds to one of the 28 possible doublet combinations.
# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
pdf(file=paste0("plot/hto_",i,"Classification_heatmap.pdf"),height=4,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
HTOHeatmap(dge.hashtag, classification="hash.ID",global.classification="HTO_classification.global",assay = "HTO", ncells = nrow(dge.hashtag@meta.data))
dev.off()

}
