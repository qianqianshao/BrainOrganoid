# scRNA-seq analysis for 4 SOSRSs at 1-month old with HTO on 1.2.2021 by Qianyi
# cluster for demultiplexed 4 SOSRSs at 1-month old


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
dge.singlet <- FindClusters(dge.singlet, resolution = seq(0.1,1,0.1), verbose = FALSE)
dge.singlet <- RunTSNE(dge.singlet, reduction = "pca", dims = 1:10)
dge.singlet <- RunUMAP(dge.singlet, reduction = "pca", dims = 1:10)


dge=dge.singlet
numPCs=10;i=1
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

# Projecting singlet identities on TSNE visualization
png(file="plot/HTO_Classify2_singlet_tSNE.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="tsne",group.by = "label")
dev.off()

png(file="plot/HTO_Classify2_singlet_UMAP.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="umap",group.by = "label")
dev.off()

# Clustering result
print(c( length(unique(dge.singlet$RNA_snn_res.0.1)),length(unique(dge.singlet$RNA_snn_res.0.2)),length(unique(dge.singlet$RNA_snn_res.0.3)),length(unique(dge.singlet$RNA_snn_res.0.4)),length(unique(dge.singlet$RNA_snn_res.0.5)),length(unique(dge.singlet$RNA_snn_res.0.6)),length(unique(dge.singlet$RNA_snn_res.0.7)),length(unique(dge.singlet$RNA_snn_res.0.8)),length(unique(dge.singlet$RNA_snn_res.0.9)),length(unique(dge.singlet$RNA_snn_res.1)) ))
print(c(mean(dge.singlet$nFeature_RNA),mean(dge.singlet$nCount_RNA),mean(dge.singlet$nCount_HTO),mean(dge.singlet$percent.mt)))


png(file="plot/HTO_Classify2_singlet_cluster.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, group.by = "RNA_snn_res.0.2",label=T)
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


# visualize Clustering result
print(c( length(unique(dge.singlet$RNA_snn_res.0.1)),length(unique(dge.singlet$RNA_snn_res.0.2)),length(unique(dge.singlet$RNA_snn_res.0.3)),length(unique(dge.singlet$RNA_snn_res.0.4)),length(unique(dge.singlet$RNA_snn_res.0.5)),length(unique(dge.singlet$RNA_snn_res.0.6)),length(unique(dge.singlet$RNA_snn_res.0.7)),length(unique(dge.singlet$RNA_snn_res.0.8)),length(unique(dge.singlet$RNA_snn_res.0.9)),length(unique(dge.singlet$RNA_snn_res.1)) ))

png(file="plot/HTO_Classify2_singlet_cluster.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet,reduction="tsne",label=T)
dev.off()

png(file="plot/HTO_Classify2_singlet_cluster_UMAP.png",res=300,height=1000,width=1200)
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
  write.table(ncellscluster,paste0("plot/ncellspercluster_JunClassifiedSinglet.txt"),quote=F,row.names=T,col.names=T,sep="\t")


# find markers
res="RNA_snn_res.0.2ordered"
Idents(dge)<-dge[[res]]
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(2),min.diff.pct=0.2)
write.table(markers,paste0("plot/HTO20kcell_6kcell_Classify2_singlet_",res,"_mindiff0.2_logfc2fold_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")
markers=FindAllMarkers(dge,only.pos=TRUE,logfc.threshold = log(1.5),min.diff.pct=0.1)
write.table(markers,paste0("plot/HTO20kcell_6kcell_Classify2_singlet_",res,"_mindiff0.1_logfc1.5fold_1.2021.txt"),col.names=T,row.names=T,quote=F,sep="\t")

markers=read.table(paste0("plot/HTO20kcell_6kcell_Classify2_singlet_",res,"_mindiff0.2_logfc2fold_1.2021.txt"),stringsAsFactors=FALSE)

plotlist=list()
print(table(Idents(dge)))
  print(i)
  print(table(markers$cluster))
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
size=sqrt(length(top2$gene))
p <- FeaturePlot(dge, top2$gene, min.cutoff = "q9", cols = c("lightblue", 
    "red"), pt.size = .5,ncol=ceiling(size), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist[[i]]=cowplot::plot_grid(plotlist = p)
}

pdf(paste0("plot/markerstop.pdf"),height=2*round(size),width=1.8*ceiling(size))
plotlist
dev.off()

# visualize known markers 1.12.2021
knownmarkers=c("PAX6","EOMES","TBR1","FOXG1","DLX2","NKX2-1","GSX2","BCL11B","SATB2","HOPX","PHETA1","VIM","FAM107A","PTPRZ1","NES","RELN","EMX2","EMX1","GLI3","SATB2")
# TBR2:    EOMES
# FAM109A: PHETA1
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
pdf("plot/clusters_ordered0_knownmarkers1_Feature.pdf",height=10,width=10)
print(plotlist)
dev.off()
png("plot/clusters_ordered0_knownmarkers1_Feature.png",res=300,height=3000,width=3000)
print(plotlist)
dev.off()

# visualize known markers 2.9.2021
knownmarkers=c("PCDH19","MKI67","FABP7","SFRP1","RSPO2","TTR","STMN2","EOMES","GAD2","CENPF","UBE2C")
# CENFP: CENPF
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 11

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

pdf("plot/clusters_ordered0_knownmarkers2_Violin.pdf",height=6,width=10)
VlnPlot(dge,gene,ncol=4,pt.size=-1)
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
pdf("plot/clusters_ordered0_knownmarkers2_Feature.pdf",height=10,width=10)
print(plotlist)
dev.off()
png("plot/clusters_ordered0_knownmarkers2_Feature.png",res=300,height=3000,width=3000)
print(plotlist)
dev.off()

# visualize known markers 2.9.2021
knownmarkers=c("RELN","BCL11B")
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 2

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

pdf("plot/clusters_ordered0_knownmarkers3_Violin.pdf",height=2,width=5)
VlnPlot(dge,gene,ncol=2,pt.size=-1)
dev.off()


# 6/30/2021
knownmarkers="ASH1L"
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 1

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

DefaultAssay(dge) <- "RNA"
pdf("plot/clusters_ordered0_ASH1L_Violin.pdf",height=2.5,width=4)
VlnPlot(dge,gene,pt.size=-1)
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

# 7/13/2021
knownmarkers="FEZF2"
knownmarkers[which(!(knownmarkers %in% rownames(dge)))]
length(knownmarkers) # 1

gene=knownmarkers[which(knownmarkers %in% rownames(dge))]

DefaultAssay(dge) <- "RNA"
pdf("plot/clusters_ordered0_FEZF2_Violin.pdf",height=2.5,width=4)
VlnPlot(dge,gene,pt.size=-1)
dev.off()
# remove legends and axes in FeaturePlot using combine=FALSE, then + NoLegend() + NoAxes()
p <- FeaturePlot(dge,gene,cols = c("lightblue", 
    "red"), combine = FALSE)
for(j in 1:length(p)) {
  p[[j]] <- p[[j]] + NoLegend() + NoAxes()
}
plotlist=cowplot::plot_grid(plotlist = p)
png("plot/clusters_ordered0_FEZF2_Feature.png",res=300,height=900,width=720)
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


## Rank cor of cluster centroids across 4 SOSRS

### using either union of HVG or top 50 markers for each cluster
res=c(paste0("RNA_snn_res.0.",c(2),"ordered"))
reps=rownames(dge@assays$HTO@data)

hvg.union=VariableFeatures(dgeall)
length(hvg.union) # 1158

top50markers=NULL
markers %>% group_by(cluster) %>% top_n(50, avg_logFC)  -> top50
top50markers=unique(top50$gene)
length(top50markers) #281

genelist=list(hvg.union,top50markers)
genelabels=c("HVG","Top50Markers")


### order cells by batch first, then by clusters of each batch
blockident=paste(dgeall$label,dgeall$RNA_snn_res.0.2ordered,sep="_")

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

print(cbind(table(Idents(dge))))

centroid=log(AverageExpression(dge)$RNA+1)
write.table(centroid,paste0("plot/HTO_singlet_rep_Centroid.txt"),row.names=T,col.names=T,quote=F,sep="\t")

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


heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))


}
dev.off()



###### Cell cycle index 
G1S=c("ACD","ACYP1","ADAMTS1","ANKRD10","APEX2","ARGLU1","ATAD2","BARD1","BRD7","C1orf63","C7orf41","C14orf142","CAPN7","CASP2","CASP8AP2","CCNE1","CCNE2","CDC6","CDC25A","CDCA7","CDCA7L","CEP57","CHAF1A","CHAF1B","CLSPN","CREBZF","CTSD","DIS3","DNAJC3","DONSON","DSCC1","DTL","E2F1","EIF2A","ESD","FAM105B","FAM122A","FLAD1","GINS2","GINS3","GMNN","HELLS","HOXB4","HRAS","HSF2","INSR","INTS8","IVNS1ABP","KIAA1147","KIAA1586","LNPEP","LUC7L3","MCM2","MCM4","MCM5","MCM6","MDM1","MED31","MRI1","MSH2","NASP","NEAT1","NKTR","NPAT","NUP43","ORC1","OSBPL6","PANK2","PCDH7","PCNA","PLCXD1","PMS1","PNN","POLD3","RAB23","RECQL4","RMI2","RNF113A","RNPC3","SEC62","SKP2","SLBP","SLC25A36","SNHG10","SRSF7","SSR3","TAF15","TIPIN","TOPBP1","TRA2A","TTC14","UBR7","UHRF1","UNG","USP53","VPS72","WDR76","ZMYND19","ZNF367","ZRANB2")                                                                                                          
S=c("ABCC5","ABHD10","ANKRD18A","ASF1B","ATAD2","BBS2","BIVM","BLM","BMI1","BRCA1","BRIP1","C5orf42","C11orf82","CALD1","CALM2","CASP2","CCDC14","CCDC84","CCDC150","CDC7","CDC45","CDCA5","CDKN2AIP","CENPM","CENPQ","CERS6","CHML","COQ9","CPNE8","CREBZF","CRLS1","DCAF16","DEPDC7","DHFR","DNA2","DNAJB4","DONSON","DSCC1","DYNC1LI2","E2F8","EIF4EBP2","ENOSF1","ESCO2","EXO1","EZH2","FAM178A","FANCA","FANCI","FEN1","GCLM","GOLGA8A","GOLGA8B","H1F0","HELLS","HIST1H2AC","HIST1H4C","INTS7","KAT2A","KAT2B","KDELC1","KIAA1598","LMO4","LYRM7","MAN1A2","MAP3K2","MASTL","MBD4","MCM8","MLF1IP","MYCBP2","NAB1","NEAT1","NFE2L2","NRD1","NSUN3","NT5DC1","NUP160","OGT","ORC3","OSGIN2","PHIP","PHTF1","PHTF2","PKMYT1","POLA1","PRIM1","PTAR1","RAD18","RAD51","RAD51AP1","RBBP8","REEP1","RFC2","RHOBTB3","RMI1","RPA2","RRM1","RRM2","RSRC2","SAP30BP","SLC38A2","SP1","SRSF5","SVIP","TOP2A","TTC31","TTLL7","TYMS","UBE2T","UBL3","USP1","ZBED5","ZWINT")                                                                             
G2M=c("ANLN","AP3D1","ARHGAP19","ARL4A","ARMC1","ASXL1","ATL2","AURKB","BCLAF1","BORA","BRD8","BUB3","C2orf69","C14orf80","CASP3","CBX5","CCDC107","CCNA2","CCNF","CDC16","CDC25C","CDCA2","CDCA3","CDCA8","CDK1","CDKN1B","CDKN2C","CDR2","CENPL","CEP350","CFD","CFLAR","CHEK2","CKAP2","CKAP2L","CYTH2","DCAF7","DHX8","DNAJB1","ENTPD5","ESPL1","FADD","FAM83D","FAN1","FANCD2","G2E3","GABPB1","GAS1","GAS2L3","H2AFX","HAUS8","HINT3","HIPK2","HJURP","HMGB2","HN1","HP1BP3","HRSP12","IFNAR1","IQGAP3","KATNA1","KCTD9","KDM4A","KIAA1524","KIF5B","KIF11","KIF20B","KIF22","KIF23","KIFC1","KLF6","KPNA2","LBR","LIX1L","LMNB1","MAD2L1","MALAT1","MELK","MGAT2","MID1","MIS18BP1","MND1","NCAPD3","NCAPH","NCOA5","NDC80","NEIL3","NFIC","NIPBL","NMB","NR3C1","NUCKS1","NUMA1","NUSAP1","PIF1","PKNOX1","POLQ","PPP1R2","PSMD11","PSRC1","RANGAP1","RCCD1","RDH11","RNF141","SAP30","SKA3","SMC4","STAT1","STIL","STK17B","SUCLG2","TFAP2A","TIMP1","TMEM99","TMPO","TNPO2","TOP2A","TRAIP","TRIM59","TRMT2A","TTF2","TUBA1A","TUBB","TUBB2A","TUBB4B","TUBD1","UACA","UBE2C","VPS25","VTA1","WSB1","ZNF587","ZNHIT2")                                        

M=c("AHI1","AKIRIN2","ANKRD40","ANLN","ANP32B","ANP32E","ARHGAP19","ARL6IP1","ASXL1","ATF7IP","AURKA","BIRC2","BIRC5","BUB1","CADM1","CCDC88A","CCDC90B","CCNA2","CCNB2","CDC20","CDC25B","CDC27","CDC42EP1","CDCA3","CENPA","CENPE","CENPF","CEP55","CFLAR","CIT","CKAP2","CKAP5","CKS1B","CKS2","CNOT10","CNTROB","CTCF","CTNNA1","CTNND1","DEPDC1","DEPDC1B","DIAPH3","DLGAP5","DNAJA1","DNAJB1","DR1","DZIP3","E2F5","ECT2","FAM64A","FOXM1","FYN","G2E3","GADD45A","GAS2L3","GOT1","GRK6","GTSE1","HCFC1","HMG20B","HMGB3","HMMR","HN1","HP1BP3","HPS4","HS2ST1","HSPA8","HSPA13","INADL","KIF2C","KIF5B","KIF14","KIF20B","KLF9","LBR","LMNA","MCM4","MDC1","MIS18BP1","MKI67","MLLT4","MZT1","NCAPD2","NCOA5","NEK2","NUF2","NUP35","NUP98","NUSAP1","ODF2","ORAOV1","PBK","PCF11","PLK1","POC1A","POM121","PPP1R10","PRPSAP1","PRR11","PSMG3","PTP4A1","PTPN9","PWP1","QRICH1","RAD51C","RANGAP1","RBM8A","RCAN1","RERE","RNF126","RNF141","RNPS1","RRP1","SEPHS1","SETD8","SFPQ","SGOL2","SHCBP1","SMARCB1","SMARCD1","SPAG5","SPTBN1","SRF","SRSF3","SS18","SUV420H1","TACC3","THRAP3","TLE3","TMEM138","TNPO1","TOMM34","TPX2","TRIP13","TSG101","TSN","TTK","TUBB4B","TXNDC9","TXNRD1","UBE2D3","USP13","USP16","VANGL1","WIBG","WSB1","YWHAH","ZC3HC1","ZFX","ZMYM1","ZNF207")   
MG1=c("AGFG1","AGPAT3","AKAP13","AMD1","ANP32E","ANTXR1","BAG3","BTBD3","CBX3","CDC42","CDK7","CDKN3","CEP70","CNIH4","CTR9","CWC15","DCP1A","DCTN6","DEXI","DKC1","DNAJB6","DSP","DYNLL1","EIF4E","ELP3","FAM60A","FAM189B","FOPNL","FOXK2","FXR1","G3BP1","GATA2","GNB1","GRPEL1","GSPT1","GTF3C4","HIF1A","HMG20B","HMGCR","HSD17B11","HSPA8","ILF2","JMJD1C","KDM5B","KIAA0586","KIF5B","KPNB1","KRAS","LARP1","LARP7","LRIF1","LYAR","MORF4L2","MRPL19","MRPS2","MRPS18B","MSL1","MTPN","NCOA3","NFIA","NFIC","NUCKS1","NUFIP2","NUP37","ODF2","OPN3","PAK1IP1","PBK","PCF11","PLIN3","PPP2CA","PPP2R2A","PPP6R3","PRC1","PSEN1","PTMS","PTTG1","RAD21","RAN","RHEB","RPL13A","SLC39A10","SNUPN","SRSF3","STAG1","SYNCRIP","TAF9","TCERG1","TLE3","TMEM138","TOB2","TOP1","TROAP","TSC22D1","TULP4","UBE2D3","VANGL1","VCL","WIPF2","WWC1","YY1","ZBTB7A","ZCCHC10","ZNF24","ZNF281","ZNF593")                                                                                             

print(c(length(G1S),length(S),length(G2M),length(M),length(MG1)))
#[1]  100 113 133 151 106
G1S=G1S[which(G1S %in% rownames(dge))]
S=S[which(S %in% rownames(dge))]
G2M=G2M[which(G2M %in% rownames(dge))]
M=M[which(M %in% rownames(dge))]
MG1=MG1[which(MG1 %in% rownames(dge))]
print(c(length(G1S),length(S),length(G2M),length(M),length(MG1)))
# [1]  99 112 133 150 106
labels=c("G1-S","S","G2-M","M","M-G1")
cellcycle=list(G1S,S,G2M,M,MG1)


#### Plot expression of cell cycle genes passing cor>0.3 with average expression pattern
cellcyclegene=cellcycle

### %expression of cell cycle genes in all cells 
exp0=colSums(as.matrix(expm1(GetAssayData(dge)[unlist(cellcyclegene), ])))/colSums(as.matrix(expm1(GetAssayData(dge))))
exp1=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[1]], ])))/colSums(as.matrix(expm1(GetAssayData(dge))))
exp2=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[2]], ])))/colSums(as.matrix(expm1(GetAssayData(dge))))
exp3=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[3]], ])))/colSums(as.matrix(expm1(GetAssayData(dge))))
exp4=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[4]], ])))/colSums(as.matrix(expm1(GetAssayData(dge))))
exp5=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[5]], ])))/colSums(as.matrix(expm1(GetAssayData(dge))))

ccinall=cbind(exp0,exp1,exp2,exp3,exp4,exp5)
write.table(ccinall,paste0("plot/cellcycle.txt"),col.names=T,row.names=T,sep="\t",quote=F)
dge$FracAllCellCycle <- exp0
which(exp0!=dge@meta.data$FracAllCellCycle) # named integer(0)
dge3=dge

###### Selecting active cycling cells by requiring % expression of all cycle genes > mean 
summary(exp0)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01738 0.08925 0.11371 0.12250 0.14125 0.59713  

### selected active cell cycling cells
length(which(exp0>mean(exp0)))/length(exp0)  # 0.4271589
length(which(exp0>mean(exp0))) # 6005
allcyclethresh=0.05
length(which(exp0>allcyclethresh))/length(exp0)  # [1] 0.6202162428
length(which(exp0>allcyclethresh))       # [1] 3926

dge=dge3
cyclingcell=colnames(dge)[which(exp0>allcyclethresh)]
dge=subset(dge3,subset = FracAllCellCycle > allcyclethresh)
dge # 27424 features across 3926 samples within 2 assays
dge2=dge
nCellperGene <- rowSums(as.matrix(GetAssayData(dge,slot="counts"))[,cyclingcell]>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use) # [1] 25977

### Number of cell cycle genes expressed in active cyling cells
# Minimum number of active cycling cells expressing each cell cycle gene
min(colSums(GetAssayData(dge)[G1S,]>0)) # 0
min(colSums(GetAssayData(dge)[S,]>0))   # 0
min(colSums(GetAssayData(dge)[G2M,]>0)) # 1
min(colSums(GetAssayData(dge)[M,]>0))   # 2
min(colSums(GetAssayData(dge)[MG1,]>0)) # 2

### %expression of all cell cycle genes in active cyling cells
exp0=colSums(as.matrix(expm1(GetAssayData(dge)[unlist(cellcyclegene), cyclingcell])))/colSums(as.matrix(expm1(GetAssayData(dge)[, cyclingcell])))
### %expression of cell cycle genes in active cycling cells
# total expression of cell cycle genes of each stage divided by total expression of all cell cycle genes
exp1=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[1]], cyclingcell])))/colSums(as.matrix(expm1(GetAssayData(dge)[unlist(cellcyclegene), cyclingcell])))
exp2=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[2]], cyclingcell])))/colSums(as.matrix(expm1(GetAssayData(dge)[unlist(cellcyclegene), cyclingcell])))
exp3=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[3]], cyclingcell])))/colSums(as.matrix(expm1(GetAssayData(dge)[unlist(cellcyclegene), cyclingcell])))
exp4=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[4]], cyclingcell])))/colSums(as.matrix(expm1(GetAssayData(dge)[unlist(cellcyclegene), cyclingcell])))
exp5=colSums(as.matrix(expm1(GetAssayData(dge)[cellcyclegene[[5]], cyclingcell])))/colSums(as.matrix(expm1(GetAssayData(dge)[unlist(cellcyclegene), cyclingcell])))

ccincyclingvssum=data.frame(GENE=names(exp0),exp0a=exp0,exp1a=exp1,exp2a=exp2,exp3a=exp3,exp4a=exp4,exp5a=exp5)

ccinallvsall=data.frame(GENE=rownames(ccinall),ccinall)
cc=merge(ccinallvsall,ccincyclingvssum,all=TRUE)
dim(cc) # 4027 13
rownames(cc)=cc[,1]
cc=as.matrix(cc[,-1])
cc=cc[rownames(ccinall),]

write.table(ccincyclingvssum,paste0("plot/cellcycle_cyclingcells.txt"),col.names=T,row.names=T,sep="\t",quote=F)
write.table(cc,paste0("plot/cellcycle5_allpluscyclingcells.txt"),col.names=T,row.names=T,sep="\t",quote=F)

dge=dge3
which(rownames(cc)!=rownames(dge@meta.data))
# integer(0)
dge$G1S <- exp1
dge$S <- exp2
dge$G2M <- exp3
dge$M <- exp4
dge$MG1 <- exp5
length(which(!(is.na(dge@meta.data$G1S)))) # [1] 3926
dge3=dge

  dge[["percent.mt"]] <- PercentageFeatureSet(dge, pattern = "^MT-")

dge.singlet=dge
dge3=dge

allcyclethresh=0.05
dge2=subset(dge, subset = FracAllCellCycle > allcyclethresh)
which(is.na(dge2$G1S))

###### per-cell attributes statistics
library(RColorBrewer)
myBrewerPalette <- brewer.pal(12,"Paired")
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
# used this color scheme for 13 clusters
tmp=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
myBrewerPalette=c(gg_color_hue(6),tmp[5:9],gg_color_hue(3),tmp[11:15])
# used this for 19 clusters once
gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette=gg_color_hue(7)

pdf(file=paste0(dgefile,"PerCellAttributes_ViolinPlot.pdf"),height=4,width=6)

dge=dge3
VlnPlot(dge, features="nFeature_RNA", ncol = 1,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="nCount_RNA", log=T, ncol = 1,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="percent.mt", ncol = 1,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="percent.x", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="percent.y", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="percent.autosome", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="GiniAll", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="GiniNon0", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="FracAllCellCycle", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')

dge=dge2
VlnPlot(dge, features="G1S", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="S", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="G2M", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="M", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="MG1", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
dev.off()



pdf(file=paste0(dgefile,"PerCellAttributes_ViolinPlot_4HTOs.pdf"),height=4,width=12)

dge=dge3
VlnPlot(dge, features="nFeature_RNA", ncol = 1,pt.size=-1,cols=myBrewerPalette)+facet_wrap(~dge$label,ncol=4)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="nCount_RNA", log=T, ncol = 1,pt.size=-1,cols=myBrewerPalette)+facet_wrap(~dge$label,ncol=4)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="percent.mt", ncol = 1,pt.size=-1,cols=myBrewerPalette)+facet_wrap(~dge$label,ncol=4)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
VlnPlot(dge, features="FracAllCellCycle", ncol = 1,pt.size =-1,cols=myBrewerPalette)+facet_wrap(~dge$label,ncol=4)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')

dge=dge2 
VlnPlot(dge, features="G1S", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+facet_wrap(~dge$label,ncol=4)+ theme(legend.position = 'none')
VlnPlot(dge, features="S", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+facet_wrap(~dge$label,ncol=4)+ theme(legend.position = 'none')
VlnPlot(dge, features="G2M", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+facet_wrap(~dge$label,ncol=4)+ theme(legend.position = 'none')
VlnPlot(dge, features="M", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+facet_wrap(~dge$label,ncol=4)+ theme(legend.position = 'none')
VlnPlot(dge, features="MG1", ncol = 1,pt.size =-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+facet_wrap(~dge$label,ncol=4)+ theme(legend.position = 'none')
dev.off()

redblue100.alpha<-rgb(read.table("/scratch/junzli_root/junzli/qzm/Dropseq_analysis/data_DGE/redblue100.txt",sep='\t',row.names=1,header=T),alpha=0.8)

dims="umap"

dge=dge2 #26694 genes across 1803 samples.
dgefile="plot/cyclingcells_"
labels=c("G1-S","S","G2-M","M","M-G1")
featurelist=list(dge@meta.data$G1S,dge@meta.data$S,dge@meta.data$G2M,dge@meta.data$M,dge@meta.data$MG1)

pdf(paste0(dgefile,"PerCellAttributes_heatmap_redblue0.8.pdf"),height=6,width=6)
par(mfrow=c(2,2),mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0),bg="grey")
for(i in 1:length(labels)){
for(dim in dims){
setname=labels[i]
feature=features.plot=featurelist[[i]]
names(feature)=rownames(dge@meta.data)

object=dge; dim.1 = 1; dim.2 = 2; cells.use = NULL; pt.size = 1;
                        pch.use = 16; reduction.use = dim;
                        use.imputed = FALSE; no.axes = TRUE; no.legend = FALSE


            cells.use <- set.ifnull(cells.use, colnames(object))
            dim.code <- translate.dim.code(reduction.use)
            dim.codes <- paste(dim.code, c(dim.1, dim.2), sep = "")
            data.plot <- FetchData(object, dim.codes)

            x1 <- paste(dim.code, dim.1, sep = "")
            x2 <- paste(dim.code, dim.2, sep = "")

            data.plot$x <- data.plot[, x1]
            data.plot$y <- data.plot[, x2]
            data.plot$pt.size <- pt.size
            #data.use <- data.frame(t(FetchData(object, features.plot, cells.use = cells.use,
            #                                   use.imputed = use.imputed)))


   data.plot$gene <- feature

st6<- data.plot
st6<-st6[order(st6[,6]),]
z<-st6[,6]

# redblue100 original
zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]

plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=paste0(setname,":",round(min(feature),2),"-",round(max(feature),2)),xlab="",ylab="")
}
}
dev.off()

