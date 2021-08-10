# scRNA-seq analysis for 6 SOSRSs at 3-month old with HTOs by Qianyi
# 6 replicates sequenced in 2 samples, each with 3 SOSRS (3 HTO tags)
# Directly-merge the 2 samples and do downstream clustering

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


# load singlet 
dgelist2=list()
for(i in 1:2){
  load(file=paste0(gexset[i],"_singlet.Robj"))
  dgelist2[[i]]=dge.singlet
}

# merge the 6 batches sequenced in 2 samples together for downstream analysis
dge1=dgelist2[[1]]
dge2=dgelist2[[2]]
dge.singlet=merge(dge1,dge2)
dge.singlet
print(c(mean(dge1$nFeature_RNA),mean(dge1$nCount_RNA),mean(dge1$nCount_HTO),mean(dge1$percent.mt)))
print(c(mean(dge2$nFeature_RNA),mean(dge2$nCount_RNA),mean(dge2$nCount_HTO),mean(dge2$percent.mt)))
print(c(mean(dge.singlet$nFeature_RNA),mean(dge.singlet$nCount_RNA),mean(dge.singlet$nCount_HTO),mean(dge.singlet$percent.mt)))

save(dge.singlet,file=paste0(sample[2],"_singlet.Robj"))

# Select the top 1000 most variable features
dge.singlet <- FindVariableFeatures(dge.singlet, selection.method = "mean.var.plot")
# Scaling RNA data, we only scale the variable features here for efficiency
dge.singlet <- ScaleData(dge.singlet, features = VariableFeatures(dge.singlet))
# Run PCA
dge.singlet <- RunPCA(dge.singlet, features = VariableFeatures(dge.singlet))

dge=dge.singlet
numPCs=12
i=1
pdf(paste("plot/dge_PCA_Variablel_variation.pdf",sep=""),height=3,width=6)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(Stdev(dge,reduction="pca")[1:40],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,lty=2)
eigenvalue=Stdev(dge,reduction="pca")[numPCs[i]]
print(eigenvalue)
plot(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(Stdev(dge,reduction="pca")),col="black")
lines(density(Stdev(dge,reduction="pca")),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.2,0.55,col="red",paste(numPCs[i],"PCs"))
dev.off()


# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
dge.singlet <- FindNeighbors(dge.singlet, reduction = "pca", dims = 1:numPCs[i])
dge.singlet <- FindClusters(dge.singlet, resolution = seq(0.1,1,0.1), verbose = FALSE)
dge.singlet <- RunTSNE(dge.singlet, reduction = "pca", dims = 1:numPCs[i])
dge.singlet <- RunUMAP(dge.singlet, reduction = "pca", dims = 1:numPCs[i])

# Projecting singlet identities on TSNE visualization
png(file="plot/hto_FitClassify_singlet_tSNE.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="tsne",group.by = "label")
dev.off()

png(file="plot/hto_FitClassify_singlet_UMAP.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="umap",group.by = "label")
dev.off()

# Clustering result
print(c( length(unique(dge.singlet$RNA_snn_res.0.1)),length(unique(dge.singlet$RNA_snn_res.0.2)),length(unique(dge.singlet$RNA_snn_res.0.3)),length(unique(dge.singlet$RNA_snn_res.0.4)),length(unique(dge.singlet$RNA_snn_res.0.5)),length(unique(dge.singlet$RNA_snn_res.0.6)),length(unique(dge.singlet$RNA_snn_res.0.7)),length(unique(dge.singlet$RNA_snn_res.0.8)),length(unique(dge.singlet$RNA_snn_res.0.9)),length(unique(dge.singlet$RNA_snn_res.1)) ))


png(file="plot/hto_FitClassify_singlet_cluster.png",res=300,height=1000,width=1100)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, group.by = "RNA_snn_res.0.2",label=T)
dev.off()

png(file="plot/hto_FitClassify_singlet_cluster_tSNE.png",res=300,height=1000,width=1100)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="tsne",group.by = "RNA_snn_res.0.2",label=T)
dev.off()

###### order clusters for each dataset
dge=dge.singlet
res=paste0("RNA_snn_res.0.",c(2))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()
resi=1
Idents(dge)<-dge[[res[resi]]]
print(c(dataset[resi],length(unique(Idents(dge)))))
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
pdf(file=paste0("plot/Centroid_norm_Seriation_",dataset[resi],"_",res[resi],".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO")))
 hmap(da) # default method="OLO"
dev.off()

### get order of seriation
 do=seriate(da,method="OLO")
levelss[[resi]]=get_order(seriate(da,method="OLO"))
# levelss=levels[get_order(seriate(da,method="OLO"))]
levelss[[resi]]
levels=levelss[[resi]]

### Reordered clusters for all cells
cells.use=colnames(dge)
# random shuffling cells within ordered clusters
ident=factor(Idents(dge),levels=levels-1)

cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i]-1)]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident.ordered=factor(as.numeric(cells),ordered=TRUE)
names(cells.ident.ordered)=cells.use

### save ordered cluster ID in dge object
which(unique(cells.ident)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"ordered")

dge[[ordered]]=cells.ident.ordered
Idents(dge)<-dge[[ordered]]
# save the dge file
dge.singlet=dge

save(dge.singlet,file=paste0(sample[2],"_singlet.Robj"))
# visualize Clustering result
print(c( length(unique(dge.singlet$RNA_snn_res.0.1)),length(unique(dge.singlet$RNA_snn_res.0.2)),length(unique(dge.singlet$RNA_snn_res.0.3)),length(unique(dge.singlet$RNA_snn_res.0.4)),length(unique(dge.singlet$RNA_snn_res.0.5)),length(unique(dge.singlet$RNA_snn_res.0.6)),length(unique(dge.singlet$RNA_snn_res.0.7)),length(unique(dge.singlet$RNA_snn_res.0.8)),length(unique(dge.singlet$RNA_snn_res.0.9)),length(unique(dge.singlet$RNA_snn_res.1)) ))

png(file="plot/hto_FitClassify_singlet_cluster.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet,reduction="tsne",label=T)
dev.off()

png(file="plot/hto_FitClassify_singlet_cluster_UMAP.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="umap",label=T)
dev.off()

### Contribution of each batch to each cluster
dge=dge.singlet
label="label"
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@meta.data[,"label"],Idents(dge))
  print(ncellsbatchcluster)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  print(percentcellsbatch)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  print(percentcellscluster)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0("plot/ncellspercluster_FitClassifiedSinglet.txt"),quote=F,row.names=T,col.names=T,sep="\t")

### compare cluster centroids with the first batch of data
dge2=dge
load("HTO20kcell_6kcell_singlet.Robj")
dge1=dge.singlet
centroid1=log(AverageExpression(dge1)$RNA+1)
centroid2=log(AverageExpression(dge2)$RNA+1)
write.table(centroid1,paste0("plot/2025-AT_singlet_7clustercentroids.txt"),quote=F,row.names=T,col.names=T,sep="\t")
write.table(centroid2,paste0("plot/2514-AT_singlet_9clustercentroids.txt"),quote=F,row.names=T,col.names=T,sep="\t")
dim(centroid1)
dim(centroid2)
# use union or intersect of HVG
hvg1=VariableFeatures(dge1)
hvg2=VariableFeatures(dge2)
union=unique(c(hvg1,hvg2))
inter=intersect(hvg1,hvg2)
length(hvg1)  # 1158
length(hvg2)  # 1047
length(union) # 1668
length(inter) # 537
length(which(union %in% rownames(dge1)))
length(which(union %in% rownames(dge2)))
# rank correlation among cluster centroids between two runs
cor(centroid2[union,],centroid1[union,],method="sp")
cor(centroid2[inter,],centroid1[inter,],method="sp")

# find markers
res="RNA_snn_res.0.2ordered"
#res="RNA_snn_res.0.6ordered"
Idents(dge)<-dge[[res]]
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
write.table(markers,paste0("plot/hto_FitClassify_singlet_",res,"_mindiff0.2_logfc2fold_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(1.5),min.diff.pct=0.1)
write.table(markers,paste0("plot/hto_FitClassify_singlet_",res,"_mindiff0.1_logfc1.5fold_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
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

pdf(paste0("plot/markerstop.pdf"),height=2*3,width=1.8*6)
plotlist
plotlist2
dev.off()

# visualize known markers 1.12.2021
knownmarkers=c("PAX6","EOMES","TBR1","FOXG1","DLX2","NKX2-1","GSX2","BCL11B","SATB2","HOPX","PHETA1","VIM","FAM107A","PTPRZ1","NES","RELN","EMX2","EMX1","GLI3","SATB2")
# TBR2:    EOMES
# FAM109A: PHETA1
# note NKX2-1 is not TTF1; TTF1 is the official gene symbol of another gene
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 20

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

pdf("plot/clusters_ordered0_knownmarkers1_Violin.pdf",height=8,width=12.5)
VlnPlot(dge,gene,ncol=5,pt.size=-1)
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
png("plot/clusters_ordered0_knownmarkers1_Feature.png",res=300,height=3000,width=3000)
print(plotlist)
dev.off()


### Visualize individual batch and subject
dgeall=dge
for(label in c("label","RNA_snn_res.0.2ordered")){
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
#dim=list(dge@dr$pca@cell.embeddings[,1:2])
#i=1
#xlim=ylim=list()
#xlim[[1]]=range(dgeall@reductions$pca@cell.embeddings[,1])
#ylim[[1]]=range(dgeall@reductions$pca@cell.embeddings[,2])
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


