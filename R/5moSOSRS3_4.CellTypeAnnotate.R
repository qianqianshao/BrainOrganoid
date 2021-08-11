# cell type annotation for 13 ordered clusters of 3 SOSRSs at 5-month old by Qianyi

# 1. scHCL reference
# https://db.cngb.org/HCL/blast.html
# https://github.com/ggjlab/scHCL

# 2. Brain and organoid reference 190k
# https://cells.ucsc.edu/?ds=organoidatlas
# ref: Tanaka et al. 2020. Cell Rep.

# 3. Cortical Organoid reference
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132672
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7433012/



R

### color scheme
gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette=gg_color_hue(3)
library(RColorBrewer)
myBrewerPalette=brewer.pal(12,"Paired")
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])

home="/scratch/junzli_root/junzli/qzm/Dropseq_analysis/"
source(paste0(home,"data_DGE/Rcode_multiplot.R"))

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


###### Data loading
library(Seurat)
indiv=2
ss=sample[indiv]
  load(file=paste0(ss,"_singlet_CCA6allgenes.Robj"))
dge=dge.singlet
dgeall=dge.singlet
raw=GetAssayData(object=dgeall,slot="counts")
norm=GetAssayData(object=dgeall)
ident=Idents(dgeall)


# 1. cell type annotation by rank correlation with cluster centroids of scHCL reference
library(scHCL)
dim(ref.expr) # 6075 1841
ref=log(ref.expr+1)
genes=rownames(ref)[which(rownames(ref) %in% rownames(norm))]
length(genes) # 4353
ref1=ref[genes,]
centroid=log(AverageExpression(dgeall)$RNA+1)
scdata=centroid
centroid1=centroid[genes,]
cors=cor(ref1,centroid1,method="spearman")
dim(cors) #[1]  1841 10

num=3
gettissue <- function(x,Num=3){
  top_markers <-order(x,decreasing = T)[1:Num]
  return(top_markers)
}
numbers_plot=num
  cors_index <- apply(cors,2,gettissue,numbers_plot)
  cors_index <- sort(unique(as.integer(cors_index)))
  length(cors_index) # 3
  scblast.result <- apply(cors,2,function(x) rownames(cors)[which.max(x)])
  cors_in = cors[cors_index,]
  colnames(cors_in)<-colnames(scdata)
  cors_out = reshape2::melt(cors_in)[c(2,1,3)]
  colnames(cors_out)<- c("Cell","Cell type","Score")
library(dplyr)
  cors_out <- as.data.frame(cors_out %>% group_by(Cell) %>% top_n(n=numbers_plot,wt=Score))
  result <- list()
  cors[which(is.na(cors))]<- 0
  result[["cors_matrix"]] <- cors
  result[['top_cors']]<- numbers_plot
  result[['scHCL']]<- scblast.result
  result[['scHCL_probility']]<- cors_out
hcl_result <- result
jpeg(file="plot/scHCL_cluster_rankcor_plot.jpeg",res=300,height=2000,width=1600)
plotHCL(hcl_result, cluster_rows=TRUE,cluster_cols=TRUE,col_font_size = 10, row_font_size=8)
dev.off()
write.table(hcl_result$cors,"plot/scHCL_cluster_rankcor_all.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(hcl_result$scHCL,"plot/scHCL_cluster_rankcor_top1.txt",row.names=T,col.names=F,quote=F,sep="\t")
write.table(hcl_result$scHCL_probility,"plot/scHCL_cluster_rankcor_top3.txt",row.names=F,col.names=T,quote=F,sep="\t")
library(pheatmap) 
# modified visualization by removing clustering in rows and/or columns                          
plotHCL <- function(hcl_result,interactive_plot=F, cluster_rows=FALSE,cluster_cols=FALSE,numbers_plot=3, col_font_size = 1, row_font_size=8, show_col=T,show_bar=T, show_tree = T){
  data(ref.expr)
  cors <- hcl_result$cors_matrix
  cors_index <- apply(cors,2,gettissue,numbers_plot)
  cors_index <- sort(unique(as.integer(cors_index)))
  data = cors[cors_index,]
  cors_out = reshape2::melt(data)[c(2,1,3)]
  colnames(cors_out)<- c("Cell","Cell type","Score")
  cors_out <- as.data.frame(cors_out %>% group_by(Cell) %>% top_n(n=numbers_plot,wt=Score))
  hcl_result$scHCL_probility <- cors_out
  hcl_result$top_cors <- numbers_plot
  height=dim(data)[1]*10+230
  tree_height = 0
  if(isTRUE(show_tree)){tree_height=50}

    p<-pheatmap(
      data,
      cluster_rows=cluster_rows,
      cluster_cols=cluster_cols,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "ward.D",
     fontsize_col = col_font_size,
      fontsize_row = row_font_size,
      color = colorRampPalette(c("grey", "white", "red"))(100),
      cellheight = 10,
      show_colnames = show_col,
      border_color = NA,
      height = height,
      legend = show_bar,
      treeheight_col = tree_height,
      treeheight_row = tree_height
      )
    if(isTRUE(interactive_plot)){

      inter_data<-data[rev(p$tree_row$order),][,p$tree_col$order]
      height= length(p$tree_row$order)*10+230
      plot_ly(x=colnames(inter_data),y=rownames(inter_data),z = inter_data, colors = colorRamp(c("grey", "white","red")),height=height, type = "heatmap", showscale=show_bar) %>% layout(autosize=T,  margin=list(l=0,r=230,b=180,t=20,pad=4),font=list(size=row_font_size),xaxis=list(showticklabels=show_col),yaxis=list(side="right"))
    }
    else{
      p
    }

}
jpeg(file="plot/scHCL_cluster_rankcor_top10_noColcluster.jpeg",res=300,height=2000,width=1600)
plotHCL(hcl_result, cluster_rows=TRUE,numbers_plot=10,cluster_cols=FALSE,col_font_size = 10, row_font_size=8)
dev.off()
                          
           
                          
# 2. cell type annotation by rank correlation with cluster centroids of organoid atlas reference
centroidA=read.table(file="organoidatlas/brainorganoid_centroid_13celltypes_2k.txt",sep="\t")
centroidC=read.table(file="organoidatlas/brainorganoid_centroid_24clusters_2k.txt")

centroid1=log(AverageExpression(dgeall)$RNA+1)

ref=centroidA
dim(ref) # [1] 2000 13
genes=rownames(ref)[which(rownames(ref) %in% rownames(norm))]
length(genes) 
ref1=ref[genes,]
scdata=centroid1
centroid11=centroid1[genes,]
cors=cor(ref1,centroid11,method="spearman")
dim(cors) 

num=3
gettissue <- function(x,Num=3){
  top_markers <-order(x,decreasing = T)[1:Num]
  return(top_markers)
}
numbers_plot=num
  cors_index <- apply(cors,2,gettissue,numbers_plot)
  cors_index <- sort(unique(as.integer(cors_index)))
  length(cors_index) # 3
  scblast.result <- apply(cors,2,function(x) rownames(cors)[which.max(x)])
  cors_in = cors[cors_index,]
  colnames(cors_in)<-colnames(scdata)
  cors_out = reshape2::melt(cors_in)[c(2,1,3)]
  colnames(cors_out)<- c("Cell","Cell type","Score")
library(dplyr)
  cors_out <- as.data.frame(cors_out %>% group_by(Cell) %>% top_n(n=numbers_plot,wt=Score))
  result <- list()
  cors[which(is.na(cors))]<- 0
  result[["cors_matrix"]] <- cors
  result[['top_cors']]<- numbers_plot
  result[['scHCL']]<- scblast.result
  result[['scHCL_probility']]<- cors_out
hcl_result <- result

write.table(hcl_result$cors,"plot/brainorganoid_celltype_rankcor_all.txt",row.names=T,col.names=T,quote=F,sep="\t")
jpeg(file="plot/brainorganoid_celltype_rankcor_top13_noColcluster.jpeg",res=300,height=2000,width=1600)
plotHCL(hcl_result, cluster_rows=TRUE,numbers_plot=13,cluster_cols=FALSE,col_font_size = 10, row_font_size=8)
dev.off()


                          
# 3. cell type annotation by rank correlation with cluster centroids of cortical organoid reference
centroidA=read.table(file="organoidatlas/corticalorganoid_centroid_17celltypes.txt",sep="\t")
centroidC=read.table(file="organoidatlas/corticalorganoid_centroid_41clusters.txt")

ref=centroidA
# 1. using HVG from our data
genes=VariableFeatures(dgeall)
genes=genes[which(genes %in% rownames(ref))]
length(genes) 

# 2. using HVG from 41 cluster centroids of ref2
g.mean <- apply(centroidC,1,mean)
g.var <- apply(centroidC,1,var)
g.var <- g.var[g.mean!=0]
g.mean <- g.mean[g.mean!=0]
summary(g.var/g.mean)
summary(g.mean)
hvg2 <- names(which(g.mean > 150 & g.var/g.mean > 30))
length(hvg2) 
genes=hvg2
genes=genes[which(genes %in% rownames(dgeall))]
length(genes) 

centroid1=log(AverageExpression(dgeall)$RNA+1)
scdata=centroid1

centroid11=centroid1[genes,]
ref1=ref[genes,]
cors=cor(ref1,centroid11,method="spearman")
dim(cors) #[1]  41 7

num=3
gettissue <- function(x,Num=3){
  top_markers <-order(x,decreasing = T)[1:Num]
  return(top_markers)
}
numbers_plot=num
  cors_index <- apply(cors,2,gettissue,numbers_plot)
  cors_index <- sort(unique(as.integer(cors_index)))
  length(cors_index) # 3
  scblast.result <- apply(cors,2,function(x) rownames(cors)[which.max(x)])
  cors_in = cors[cors_index,]
  colnames(cors_in)<-colnames(scdata)
  cors_out = reshape2::melt(cors_in)[c(2,1,3)]
  colnames(cors_out)<- c("Cell","Cell type","Score")
library(dplyr)
  cors_out <- as.data.frame(cors_out %>% group_by(Cell) %>% top_n(n=numbers_plot,wt=Score))
  result <- list()
  cors[which(is.na(cors))]<- 0
  result[["cors_matrix"]] <- cors
  result[['top_cors']]<- numbers_plot
  result[['scHCL']]<- scblast.result
  result[['scHCL_probility']]<- cors_out
hcl_result <- result
write.table(hcl_result$cors,"plot/corticalorganoid_celltype_rankcor_all.txt",row.names=T,col.names=T,quote=F,sep="\t")

jpeg(file="plot/corticalorganoid_celltype_rankcor_top10_noColcluster.jpeg",res=300,height=2000,width=1600)
plotHCL(hcl_result, cluster_rows=TRUE,numbers_plot=10,cluster_cols=FALSE,col_font_size = 10, row_font_size=8)
dev.off()
