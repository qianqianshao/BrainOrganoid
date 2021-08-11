# scRNA-seq analysis for 6 SOSRSs at 3-month old with HTOs on 5.12.2021 by Qianyi
# 3 replicates sequenced in 2 samples of 1 run, each with 3 SOSRS (3 HTO tags)
# sample 1: HTO-1,2,3
# sample 2: HTO-4,5,6 (same as HTO-1,2,3, respectively)
# CCA-integrate the 6 replicates and do downstream clustering

R

library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(gplots)
library(patchwork)

home="/scratch/junzli_root/junzli/qzm/Dropseq_analysis/"
source(paste0(home,"data_DGE/Rcode_multiplot.R"))
redblue100<-rgb(read.table(paste0(home,'data_DGE/redblue100.txt'),sep='\t',row.names=1,header=T))
gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette1=gg_color_hue(6)
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(7,"Set2"),brewer.pal(12,"Paired")[c(10,12)])
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])

dgefile="/gpfs/accounts/junzli_root/junzli/qzm/Dropseq_analysis/HTO/"
setwd(dgefile)
names=c("3GEX","HTO")
run=c("2931-AT","2931-AT")
sample=c("2931-AT","2931-AT")
gexset=c("2931-AT-1","2931-AT-2")
htoset=list(
  1:3,
  4:6
)
dataset=list(
  c("2931-AT-1-GEX_GTCCCATC-GTCACGTT","2931-AT-1-HTO"),
  c("2931-AT-2-GEX_CCGGAGGA-AACATCCG","2931-AT-2-HTO")
)
subject=c("Human4","Human4")
organ=c("Brain","Brain")
n=length(dataset)
datainfo=data.frame(run,dataset,names,subject,organ)




# load directly-merged singlet
  load(file=paste0(sample[2],"_singlet.Robj"))
dge=dge.singlet
table(dge$label)
HTO-1 HTO-2 HTO-3 HTO-4 HTO-5 HTO-6 
 1972  1709  1395  1003   942   770 

### removed non-detected genes 
nUMIperGene=apply(dge@assays$RNA@data,1,sum) 
genes.use=names(nUMIperGene)[which(nUMIperGene>0)]
length(nUMIperGene)  # 29488
length(genes.use)    # 28811

### separate 6 replicates - 6 HTO tags
datalists=SplitObject(dge,split.by="label")

datalists
datalists=datalists[c(3,1,2,5,6,4)]
datalists


datalist <- lapply(X = datalists, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})



dgefile="/gpfs/accounts/junzli_root/junzli/qzm/Dropseq_analysis/HTO/plot/CCA6_"


### Integration 
k.filter <- min(200, min(sapply(datalist, ncol)))
k.filter # revised k.filter based on the number of cells in the smallest dataset
anchors <- FindIntegrationAnchors(object.list = datalist,k.filter=k.filter, dims = 1:30) # anchor.features =length(genes.use),
save(anchors,file=paste0("plot/CCA6-anchors.Robj"))

load(file=paste0("plot/CCA6-anchors.Robj"))
integrated <- IntegrateData(anchorset = anchors, features.to.integrate = genes.use, dims = 1:30)

save(integrated,file=paste0(sample[2],"_singlet_CCA6allgenes.Robj"))


### Run PCA using integrated data
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
dge <- integrated
dge <- ProjectDim(dge)

save(integrated,file=paste0(sample[2],"_singlet_CCA6allgenes.Robj"))


###### Determine top PCs
numPCs=c(8,11);
i=2
pdf(paste("plot/CCA6-allgenes_PCA_Variablel_variation_",numPCs[i],".pdf",sep=""),height=4,width=8)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(Stdev(dge,reduction="pca")[1:30],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,lty=2)
### density plot of Eigenvalue
eigenvalue=Stdev(dge,reduction="pca")[numPCs[i]]
print(eigenvalue)
plot(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(Stdev(dge,reduction="pca")),col="black")
lines(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.2,0.55,col="red",paste(numPCs[i],"PCs"))
dev.off()

###### Using top PCs
DefaultAssay(dge) <- "integrated"
### tSNE
  dge <- RunTSNE(dge, dims = 1:numPCs[i])
  dge <- RunUMAP(dge, dims = 1:numPCs[i])
### Louvain-Jaccard Clustering
dge <- FindNeighbors(dge,dims=1:numPCs[i])
dge <- FindClusters(dge, resolution = seq(0.1,1.2,by=0.1))

save(dge,file=paste0(sample[2],"_singlet_CCA6allgenes_",numPCs[i],"PCs.Robj"))


print(c( length(unique(dge$integrated_snn_res.0.1)),length(unique(dge$integrated_snn_res.0.2)),length(unique(dge$integrated_snn_res.0.3)),length(unique(dge$integrated_snn_res.0.4)),length(unique(dge$integrated_snn_res.0.5)),length(unique(dge$integrated_snn_res.0.6)),length(unique(dge$integrated_snn_res.0.7)),length(unique(dge$integrated_snn_res.0.8)),length(unique(dge$integrated_snn_res.0.9)),length(unique(dge$integrated_snn_res.1)),length(unique(dge$integrated_snn_res.1.1)),length(unique(dge$integrated_snn_res.1.2)) ))


###### order clusters for each individual organ
res=c(paste0("integrated_snn_res.0.",c(1)),paste0("integrated_snn_res.0.",c(3)));
resi=2


## order cell by cluster ID and randomly shuffle cells within each batch
Idents(dge)<-dge[[res[resi]]]
print(c(organs[resi],length(unique(Idents(dge)))))
levels=levels(Idents(dge))

centroid=log(AverageExpression(dge)$RNA+1)

### Reordering cluster centroid using dissimilarity matrix
library(seriation)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(dataset)
bb=1
tmp=centroid[,levels]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
# note: dist calculate distance between each row
length(da) # 91
da
pdf(file=paste0("plot/Centroid_norm_Seriation_",res[resi],".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO")))
 hmap(da) # default method="OLO"
dev.off()

### get order of seriation
 do=seriate(da,method="OLO")
levelss=get_order(seriate(da,method="OLO"))

levels=levelss

### Reordered clusters for all cells
cells.use=colnames(dge)
ident=factor(Idents(dge),levels=levels-1)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i]-1)]
   if(length(tmp)>0){
      tmpname=names(tmp)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident.ordered=factor(as.numeric(cells),ordered=TRUE)
names(cells.ident.ordered)=cells.use

### re-order based on previous ordering
if(resi==1){
levels=levelss[c(1:4,7:5)]
}
if(resi==2){
levels=rev(levelss[c(10:13,9:8,1:7)])
}

### Reordered clusters for all cells
cells.use=colnames(dge)
Idents(dge)<-dge[[res[resi]]]
ident=factor(Idents(dge),levels=levels-1)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i]-1)]
   if(length(tmp)>0){
      tmpname=names(tmp)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident.ordered=factor(as.numeric(cells),ordered=TRUE)
names(cells.ident.ordered)=cells.use


### save ordered cluster ID in dge object
which(unique(cells.ident.ordered)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"ordered")

dge[[ordered]]=cells.ident.ordered
Idents(dge)<-dge[[ordered]]
# save the dge file

save(dge,file=paste0(sample[2],"_singlet_CCA6allgenes_",numPCs[resi],"PCs.Robj"))



res=paste0("integrated_snn_res.0.",c(1,3),"ordered")
i=resi=1
dge$run2CCA6=dge$integrated_snn_res.0.1ordered
dge.singlet=dgeall=dge

i=resi=2
Idents(dge) <- dge[[res[resi]]]
dge$run2CCA6cc13=dge$integrated_snn_res.0.3ordered

id=as.character(dge$label)
names(id)=rownames(dge@meta.data)
id[which(id=="HTO-1")]<- "HTO-1/4"
id[which(id=="HTO-2")]<- "HTO-2/5"
id[which(id=="HTO-3")]<- "HTO-3/6"
id[which(id=="HTO-4")]<- "HTO-1/4"
id[which(id=="HTO-5")]<- "HTO-2/5"
id[which(id=="HTO-6")]<- "HTO-3/6"
id=as.factor(id)
table(id)
dge$sample=id

png(file="plot/CCA6hto_FitClassify_singlet_tSNE.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="tsne",group.by = "label")
dev.off()

png(file="plot/CCA6hto_FitClassify_singlet_UMAP.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="umap",group.by = "label")
dev.off()


png(file="plot/CCA6hto_FitClassify_singlet_cluster.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet,reduction="tsne",label=T,cols=myBrewerPalette)
dev.off()

png(file="plot/CCA6hto_FitClassify_singlet_cluster_UMAP.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="umap",label=T,cols=myBrewerPalette)
dev.off()


### Contribution of each batch to each cluster
dge=dge.singlet
label="label"
label="sample"
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@meta.data[,label],Idents(dge))
  print(ncellsbatchcluster)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  print(percentcellsbatch)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  print(percentcellscluster)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0("plot/ncellspercluster_FitClassifiedSinglet_",label,".txt"),quote=F,row.names=T,col.names=T,sep="\t")


pdf(file=paste0(dgefile,"PerCellAttributes_ViolinPlot.pdf"),height=4,width=6)
Idents(dge)<-dge[[res[resi]]]
VlnPlot(dge, features="nFeature_RNA", ncol = 1,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="nCount_RNA", log=T, ncol = 1,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
dev.off()


### compare cluster centroids with the first batch of data
table(dge$run2directmerge,dge$run2CCA6)

dge5=dge
centroid5=log(AverageExpression(dge5)$RNA+1)
dim(centroid5)

load("HTO20kcell_6kcell_singlet.Robj")
dge.singlet
dge1=dge.singlet

load("2514-AT_singlet.Robj")
dge.singlet
dge2=dge.singlet

centroid1=read.table("plot/2025-AT_singlet_7clustercentroids.txt")
centroid2=read.table("plot/2514-AT_singlet_CCA6_7clustercentroids.txt")
dim(centroid1)
dim(centroid2)

# use union or intersect of HVG
hvg1=VariableFeatures(dge1)
hvg2=VariableFeatures(dge2)
hvg5=VariableFeatures(dge5)

union=unique(c(hvg5,hvg1))
union=union[which(union %in% intersect(rownames(dge5),rownames(dge1)))]
inter=intersect(hvg5,hvg1)
length(union) # 1515
length(inter) # 619
# rank correlation among cluster centroids between two runs
cor(centroid5[union,],centroid1[union,],method="sp")
cor(centroid5[inter,],centroid1[inter,],method="sp")

union=unique(c(hvg5,hvg2))
union=union[which(union %in% intersect(rownames(dge5),rownames(dge2)))]
inter=intersect(hvg5,hvg2)
length(union) # 2525
length(inter) # 629
# rank correlation among cluster centroids between two runs
cor(centroid5[union,],centroid2[union,],method="sp")
cor(centroid5[inter,],centroid2[inter,],method="sp")


# find differentially-expressed markers
DefaultAssay(dge) <- "RNA"        # should use this
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
write.table(markers,paste0("plot/hto_FitClassify_singlet_CCA6_",res[resi],"_mindiff0.2_logfc2fold_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
  print(table(markers$cluster))
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(1.5),min.diff.pct=0.1)
write.table(markers,paste0("plot/hto_FitClassify_singlet_CCA6_",res[resi],"_mindiff0.1_logfc1.5fold_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
  print(table(markers$cluster))

print(table(Idents(dge)))
  print(i)
  print(table(markers$cluster))


plotlist=plotlist2=list()
print(table(Idents(dge)))
  print(i)
  print(table(markers$cluster))
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
size=sqrt(length(top2$gene))
p <- FeaturePlot(dge, top2$gene, min.cutoff = "q9", cols = c("lightblue", 
    "red"), pt.size = .5,ncol=ceiling(size), combine = FALSE)
p2 <- FeaturePlot(dge, top2$gene, reduction="tsne",min.cutoff = "q9", cols = c("lightblue", 
    "red"), pt.size = .5,ncol=ceiling(size), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
  p2[[j]] <- p2[[j]] + NoLegend() + NoAxes()
}
plotlist[[i]]=cowplot::plot_grid(plotlist = p,ncol=6)
plotlist2[[i]]=cowplot::plot_grid(plotlist = p2,ncol=6)
}

png(paste0("plot/markerstop.png"),res=300,height=600*5,width=500*6)
plotlist
#plotlist2
dev.off()

# visualize known markers 3.12.2021
knownmarkers=c("EBF1","DMRTA2","RELN","NR2F2","EIF4E","LHX1","LHX5")
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 7

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

pdf("plot/clusters_ordered0_knownmarkers1_Violin.pdf",height=4,width=10)
VlnPlot(dge,gene,ncol=4,pt.size=-1)
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
png("plot/clusters_ordered0_knownmarkers1_Feature.png",res=300,height=1800,width=1600)
print(plotlist)
dev.off()

knownmarkers=c("HOPX","EOMES","FEZF2","CALB2","DLX5","LHX5","LHX6","CRYAB","TBR1","SATB2","PROX1","RORB","GAD1","GAD2")
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 14

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

pdf("plot/clusters_ordered0_knownmarkers2_Violin.pdf",height=6,width=12.5)
VlnPlot(dge,gene,ncol=5,pt.size=-1)
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
png("plot/clusters_ordered0_knownmarkers2_Feature.png",res=300,height=3600,width=3200)
print(plotlist)
dev.off()


knownmarkers=c("VIM","NES","FABP7","EMX2","SFRP1","RSPO2","TTR","STMN2","GSX2","CENPF","MKI67","BCL11B","EBF2")
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 13

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

pdf("plot/clusters_ordered0_knownmarkers3_Violin.pdf",height=6,width=12.5)
VlnPlot(dge,gene,ncol=5,pt.size=-1)
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
png("plot/clusters_ordered0_knownmarkers3_Feature.png",res=300,height=3600,width=3200)
print(plotlist)
dev.off()

knownmarkers=c("AQP4","GFAP","S100B","ALDH1L1","OLIG2","MBP","PAX6","FOXG1")
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 8

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

pdf("plot/clusters_ordered0_knownmarkers4_Violin.pdf",height=4,width=10)
VlnPlot(dge,gene,ncol=4,pt.size=-1)
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
png("plot/clusters_ordered0_knownmarkers4_Feature.png",res=300,height=1800,width=1600)
print(plotlist)
dev.off()



knownmarkers=c("VIM","FABP7","HOPX","EMX2","SFRP1","RSPO2","TTR","RELN","STMN2","EOMES","TBR1","GAD2","BCL11B","FEZF2","SATB2","MKI67","FOXG1","PAX6")
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 18

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

pdf("plot/clusters_ordered0_knownmarkers5_Violin.pdf",height=8,width=12.5)
VlnPlot(dge,gene,ncol=5,pt.size=-1)
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
png("plot/clusters_ordered0_knownmarkers5_Feature.png",res=300,height=3600,width=3600)
print(plotlist)
dev.off()


# 5/21/2021
knownmarkers=c("RELN","LHX1","EBF3","POU3F2","POU3F3","CUX1","CUX2","SATB2","RORB","NECAB1","BCL11B","FEZF2","TLE4","TBR1","EOMES","PAX6","FABP7","SFRP1","HOPX","TTR","RSPO2","AQP4","GFAP","S100B","OLIG2","MBP","MKI67","CENPF","VIM","STMN2","FOXG1","GAD2")
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 32

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

DefaultAssay(dge) <- "RNA"
pdf("plot/clusters_ordered0_knownmarkers6_Violin.pdf",height=10,width=17.5)
VlnPlot(dge,gene,ncol=7,pt.size=-1,cols=myBrewerPalette)
dev.off()
DefaultAssay(dge) <- "integrated"
pdf("plot/clusters_ordered0_knownmarkers_DotPlot_RNA.pdf",height=5,width=15)
DotPlot(dge, features = gene,cols=c("white","blue")) + RotatedAxis()
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
png("plot/clusters_ordered0_knownmarkers6_Feature.png",res=300,height=4500,width=4320)
print(plotlist)
dev.off()


# 6/30/2021
knownmarkers="ASH1L"
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 1

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

DefaultAssay(dge) <- "RNA"
pdf("plot/clusters_ordered0_ASH1L_Violin.pdf",height=2.5,width=4)
VlnPlot(dge,gene,pt.size=-1,cols=myBrewerPalette)
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
png("plot/clusters_ordered0_ASH1L_Feature.png",res=300,height=900,width=720)
print(plotlist)
dev.off()



### Visualize individual batch and subject
dgeall=dge
for(label in c("sample","label",res[resi])){
plotlist=list()
Idents(dgeall)<-dgeall[[label]]
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlim=ylim=list()
pos=c("topleft","topright","bottomright","bottomright")
# dge4 used top 9 PCs for tSNE
xlim=list(range(Embeddings(dge,"pca")[,1]),range(Embeddings(dge,"pca")[,1]),range(Embeddings(dge,"umap")[,1]),range(Embeddings(dge,"tsne")[,1]))
ylim=list(range(Embeddings(dge,"pca")[,2]),range(Embeddings(dge,"pca")[,3]),range(Embeddings(dge,"umap")[,2]),range(Embeddings(dge,"tsne")[,2]))
### plot PCs and tSNE for each batch using the other batches as background
dge=dgeall
sets=levels(Idents(dgeall))
if(length(sets)>1){
  cols=gg_color_hue(length(sets))
plot2set=plot3set=plot4set=plottset=NULL
dim=list(Embeddings(dge,"pca")[,1:2],Embeddings(dge,"pca")[,c(1,3)],Embeddings(dge,"umap"),Embeddings(dge,"tsne"))
size=length(sets)
if(size>4){
  pdf(paste0("plot/",label,"_indiv.pdf"),height=2.3*round(sqrt(size)),width=2.3*ceiling(sqrt(size)))
  par(mfrow=c(round(sqrt(size)),ceiling(sqrt(size))),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
} else {
  pdf(paste0("plot/",label,"_indiv.pdf"),height=2.5,width=2.5*size)
  par(mfrow=c(1,size),mar=c(2.5,2.5,0.5,0.5),mgp=c(1.2, 0.5, 0))
}
for(j in 1:4){
#  j=1
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(Idents(dgeall))
names(ident)=names(Idents(dgeall))
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.4,col=c("grey70",alpha(cols[seti],0.8))[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}
if(size<(round(sqrt(size))*ceiling(sqrt(size)))){
  for(new in (size+1):(round(sqrt(size))*ceiling(sqrt(size)))){
    plot.new()
  }
}
}
dev.off()
}
}
}


## Rank cor

### using either union of HVG or top 50 markers for each cluster
res=c(paste0("integrated_snn_res.0.",c(1),"ordered"))
reps=rownames(dge@assays$HTO@data)

hvg.union=VariableFeatures(dge.singlet)
length(hvg.union) # 2000

top50markers=NULL
markers=read.table(paste0("plot/hto_FitClassify_singlet_CCA6_",res,"_mindiff0.1_logfc1.5fold_1.2021.txt"),stringsAsFactors=FALSE)
markers %>% group_by(cluster) %>% top_n(50, avg_logFC)  -> top50
top50markers=unique(top50$gene)
length(top50markers) #267

genelist=list(hvg.union,top50markers)
genelabels=c("HVG","Top50Markers")


### order cells by batch first, then by clusters of each batch
blockident=paste(dgeall$label,dgeall$integrated_snn_res.0.1ordered,sep="_")

### Clusters ordered first by batches, then by res
batch=reps
nbatch=length(batch)
ncluster=NULL
for(i in batch){
  ncluster=c(ncluster,length(unique(Idents(dgeall))))
}
ncluster 
clusters=list()
for(i in 1:nbatch){
  clusters[[i]]=rep(1:ncluster[i])
}
levels2=NULL
for(bb in 1:nbatch){
    cluster=clusters[[bb]]
    levels2=c(levels2,paste(batch[bb],cluster,sep="_"))
}
levels2=levels2[which(levels2 %in% unique(blockident))]
levels=levels2

#names(blockident)=paste0("X",names(blockident))
#names(blockident)=gsub("-",".",names(blockident))
ident=factor(blockident,levels=levels)

### Calculate correlation for each normalized centroid using HVG
### for each cluster, calculate average normalized expression of each gene
dge=dgeall
Idents(dge)<-ident
dge$indivclusters<-ident

dgeall=dge
dge.singlet=dgeall
save(dge,file=paste0(sample[2],"_singlet_CCA6allgenes.Robj"))


print(cbind(table(Idents(dge))))

centroid=log(AverageExpression(dge)$RNA+1)
write.table(centroid,paste0("plot/HTO_singlet_CCA6_rep_Centroid.txt"),row.names=T,col.names=T,quote=F,sep="\t")

#centroid1=read.table(paste0("plot/",dataset[1],"_",res[1],"_Centroid.txt"),header=T,row.names=1)
#centroid[1:5,1:5]
#centroid1[1:5,1:5]

pdf(file=paste("plot/organs_parts_Centroid_RankedCorrelation_HVG.pdf",sep=""),height=7.5,width=7)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))

for(g in 1:length(genelist)){
genelabel=genelabels[g]
genes=genelist[[g]]
print(length(genes))

cc=cor(as.matrix(centroid)[genes,],method="spearman")
dim(cc)
min(cc) # 0.036

data.use=cc[levels,levels]
### save the cluster centroid rank correlation using HVG
write.table(data.use,paste0("plot/HTO_singlet_rep_Centroid_rho_",genelabel,".txt"),row.names=T,col.names=T,quote=F,sep="\t")
### note: once calculated, next time can directly read table without calculating again

### load cluster centroid rank correlation using HVG
data.use=read.table(paste0("plot/HTO_singlet_rep_Centroid_rho_",genelabel,".txt"),header=T,row.names=1)
colnames(data.use)=rownames(data.use)
levels=rownames(data.use)
#batch=names[ft]

### labeling
colsep.use=cumsum(table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))])
col.lab=rep("",length(levels))
col.lab[round(cumsum(table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))])-table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))]/2)+0]=unique(gsub("_.*","",levels))
row.lab=gsub(".*_","",levels)

sidecol=do.call(rbind,strsplit(levels,"_"))
batchtmp=batch[which(batch %in% unique(sidecol[,1]))]
for(rep in 1:length(unique(sidecol[,1]))){
a=batchtmp[rep]
sidecol[which(sidecol[,1]==a),1]<-rep
}

rlab=matrix(0,2,length(levels))
rlab[1,]=rep(c("white","black"),6)[as.numeric(sidecol[,1])]
for(i in 1:nrow(sidecol)){
  rlab[2,i]=myBrewerPalette[as.numeric(sidecol[i,2])]
}
clab=cbind(rlab[2,],rlab[1,])
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")

col.use=redblue100

library(gplots)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=.8,cexRow=.8,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))


}
dev.off()
