# cell type annotation for 7 ordered clusters of 4 SOSRSs at 1-month old on 2.4.2021 by Qianyi

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
names=c("3GEX","HTO")
run=c("NovaA-320","NovaA-320")
sample=c("2025-AT","2025-AT")
dataset=c("2025-AT-1-3GEX_CGTGACAT-TTTAGACC","2025-AT-1-HTO")
subject=c("Human1","Human1")
organ=c("Brain","Brain")
n=length(dataset)
datainfo=data.frame(run,dataset,names,subject,organ)

### Data loading
library(Seurat)
indiv=1
ss="1mo_6kcell_singlet"
load(file=paste0(ss,".Robj"))
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
# wget https://cells.ucsc.edu/organoidatlas/exprMatrix.tsv.gz
# wget https://cells.ucsc.edu/organoidatlas/meta.tsv
library(dplyr)
library(Seurat)
# use brain organoid data as reference
# metadata downloaded from https://cells.ucsc.edu/?ds=organoidatlas
meta <- read.table("organoidatlas/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
dim(meta) # [1] 190022      8
meta[1:2,]
                       nCount_RNA nFeature_RNA Cluster Cluster_name
Xiang_AAACATACCGGAGA-1       2480         1014       3          CN3
Xiang_AAACATACCTCGAA-1       6105         1994      15          AS2
                            Annotation       Dataset Protocol Age
Xiang_AAACATACCGGAGA-1 Cortical neuron Xiang 1m rep1    Xiang  1m
Xiang_AAACATACCTCGAA-1       Astrocyte Xiang 1m rep1    Xiang  1m

table(meta$Cluster)
    1     2     3     4     5     6     7     8     9    10    11    12    13 
 6197 16062 24989  7264  5248  4317  5341  5054  5153  8202  2414  9179  7068 
   14    15    16    17    18    19    20    21    22    23    24 
 8309 16256 10619 11367  2151  8550  1803  9852  5591  8753   283 
length(table(meta$Cluster_name)) # 24
length(table(meta$Annotation))   # 13
table(meta$Cluster_name,meta$Cluster) # 1-to-1 match
table(meta$Cluster_name,meta$Annotation) 
# combined 24 clusters to 13 major cell types
3 AS clusters as Astrocyte 
5 CN clusters as Cortical neuron
4 NEC clusters as Neuroepithelial cell
2 PGC clusters as Proteoglycan-expressing cell
2 UPRC clusters as unfolded protein responsible cell

# gene expression matrix #1. downloaded from https://cells.ucsc.edu/?ds=organoidatlas
library(data.table)
#mat <- fread("organoidatlas/exprMatrix.tsv.gz") gives error
mat <- read.table("organoidatlas/exprMatrix.tsv.gz")
dim(mat) # 2001 190023
mat[1:5,1:3]
genes = as.character(mat[,1])
genes = gsub(".+[|]", "", genes)
mat = data.frame(mat[,-1], row.names=genes)
dim(mat) # 2000 190022
tmp = t(mat)
tmp[1:5,1:3]
cells=tmp[,1]
colnames(mat)=cells
which(names(mat) != mat[1,]) # integer(0)
mat=mat[-1,]
dim(mat) # 2000 190022
mat[1:5,1:3]
      Xiang_AAACATACCGGAGA-1 Xiang_AAACATACCTCGAA-1 Xiang_AAACATACGGCGAA-1
UBE2C     0.0657753736989895     0.0826156240986122       1.60758241630298
TTR         1.42903881254195       2.87464653571497       2.24589722469453
GNRH1     0.0105743000665074     0.0780825795358962     0.0234568984944464
TNNC2   -0.00174491389132081    -0.0199264602041143    0.00832622304096964
TOP2A     0.0768601254774907    -0.0095935677037195       1.35816811537279

so <- CreateSeuratObject(counts = mat, project = "OrganoidAtlas", meta.data=meta)

save(so,file="organoidatlas/brainorganoid.Robj")

m1=apply(so@assays$RNA@data,1,mean)
sd1=apply(so@assays$RNA@data,1,sd)
m2=apply(so@assays$RNA@data,2,mean)
sd2=apply(so@assays$RNA@data,2,sd)
png(file="organoidatlas/brainorganoid_meansd.png",height=400,width=800)
par(mfrow=c(1,2),mar=c(4,4,2,0.5),mgp=c(2,0.5,0))
plot(m1,sd1,pch=16,cex=.5,col=rgb(0,0,0,0.3),main="Per Gene Mean and SD")
plot(m2,sd2,pch=16,cex=.5,col=rgb(0,0,0,0.3),main="Per Cell Mean and SD")
dev.off()


dge=so
pdf(file=paste0("organoidatlas/brainorganoid_PerCellAttributes_ViolinPlot.pdf"),height=8,width=6)
  plotlist=list()
Idents(dge) <- "Protocol"
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 1)
Idents(dge) <- "Annotation"
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 1)
Idents(dge) <- "Protocol"
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",log=T,pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 1)
Idents(dge) <- "Annotation"
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 1)
dev.off()

protocol2=as.character(dge$Protocol)
names(protocol2)=rownames(dge@meta.data)
protocol2[which(protocol2!="Fetal brain")] <- "Others"
dge$protocol2=as.factor(protocol2)
VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+facet_wrap(~dge$protocol2)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')

pdf(file=paste0("organoidatlas/brainorganoid_PerCellAttributes_ViolinPlot_facet.pdf"),height=12,width=18)
Idents(dge) <- "Annotation"
VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+facet_wrap(~dge$Protocol)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
dev.off()

plotlist=list()
dge$log10nGene=log10(dge$nFeature_RNA)
dge$log10nUMI=log10(dge$nCount_RNA)
png(file=paste0("organoidatlas/brainorganoid_PerCellAttributes_FeatureScatter.png"),res=300,height=2800,width=2000)
Idents(dge) <- "Annotation"
plotlist[[1]]=FeatureScatter(dge,"log10nGene","log10nUMI",pt.size=.5,cols=myBrewerPalette)
Idents(dge) <- "Protocol"
plotlist[[2]]=FeatureScatter(dge,"log10nGene","log10nUMI",pt.size=.5)
multiplot(plotlist,cols = 1)
dev.off()


table(dge$Protocol)
      Birey Fetal brain      Fiddes Giandomenic    Madhavan    Quadrato 
       4953        2338       13747       13009        3743       43706 
   Trujillo     Velasco       Xiang 
       3295       81524       23707 

table(dge$Annotation,dge$Protocol)
                                    Birey Fetal brain Fiddes Giandomenic
  Astrocyte                           588         203   2728        1090
  BMP responsible cell                 53          63    667          43
  Cilia-bearing cell                   13           4     16          15
  Cortical neuron                    2699        1253   1262        7737
  Glia progenitor cell                182          45   1226         494
  Intermediate                         45          68   1758          73
  Interneuron                          76         245    235         773
  Mesoderm                              4          13     30           0
  Neuroepithelial cell                889         237   4621        1174
  Neuron                              220          84    179         642
  Oligodendrocyte/OPC                  97          47    382         326
  Proteoglycan-expressing cell          9          13    340          33
  Unfolded protein responsible cell    78          63    303         609
                                   
                                    Madhavan Quadrato Trujillo Velasco Xiang
  Astrocyte                              517     6588      553   17465  5452
  BMP responsible cell                    48     3309      198    1114  3055
  Cilia-bearing cell                      11     1374        4      78   636
  Cortical neuron                       1527     5794     1421   32940  5127
  Glia progenitor cell                   421     2256      163    4641  1939
  Intermediate                            23     1221       18    1098   750
  Interneuron                            313      443       96    1449   687
  Mesoderm                                 5     1547        0     186    18
  Neuroepithelial cell                    99     4268      480    9978  3202
  Neuron                                 173      522      157    2778   586
  Oligodendrocyte/OPC                    139     1545       97    3693   742
  Proteoglycan-expressing cell             6    13782       24    1044   192
  Unfolded protein responsible cell      461     1057       84    5060  1321

dge=so
Idents(dge) <- "Annotation"
centroidA=log(AverageExpression(dge)$RNA+1)
write.table(centroidA,file="organoidatlas/brainorganoid_centroid_13celltypes_2k.txt",row.names=T,col.names=T,quote=F,sep="\t")
# use this as reference cell type centroids
Idents(dge) <- "Cluster_name"
centroidC=log(AverageExpression(dge)$RNA+1)
write.table(centroidC,file="organoidatlas/brainorganoid_centroid_24clusters_2k.txt",row.names=T,col.names=T,quote=F,sep="\t")
# use this as reference cell type centroids


centroidA=read.table(file="organoidatlas/brainorganoid_centroid_13celltypes_2k.txt",sep="\t")
centroidC=read.table(file="organoidatlas/brainorganoid_centroid_24clusters_2k.txt")

Idents(dgeall) <- dgeall$RNA_snn_res.0.2ordered
centroid1=log(AverageExpression(dgeall)$RNA+1)

ref=centroidA
dim(ref) # [1] 2000 13
genes=rownames(ref)[which(rownames(ref) %in% rownames(norm))]
length(genes) # 1919
ref1=ref[genes,]
scdata=centroid1
centroid11=centroid1[genes,]
cors=cor(ref1,centroid11,method="spearman")
dim(cors) # 24 7

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
# wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132672/suppl/GSE132672_allorganoids_withnew_matrix.txt.gz
library(dplyr)
library(Seurat)
# use cortical organoid data as reference
meta <- read.table("organoidatlas/Bhaduri_Organoid_metadata.txt", header=T, sep="\t", row.names=1)
dim(meta) # [1] 242349     14

meta[1:2,]
                          nGene      nUMI   Sample Line       Protocol Age
H1SWeek3_AAACCTGAGACAAAGG  5141 305542.57 H1SWeek3   H1 Least Directed   3
H1SWeek3_AAACCTGAGCACACAG  1748  53458.58 H1SWeek3   H1 Least Directed   3
                          Origin  Batch iPSC.or.hESC Cluster        Class
H1SWeek3_AAACCTGAGACAAAGG   hESC Batch1         hESC      29 Non-neuronal
H1SWeek3_AAACCTGAGCACACAG   hESC Batch1         hESC       5 Non-neuronal
                                 State        Type         Subtype
H1SWeek3_AAACCTGAGACAAAGG     Dividing Radial Glia Pan-radial glia
H1SWeek3_AAACCTGAGCACACAG Non-dividing Radial Glia    Hindbrain RG


table(meta$Cluster) # 41
table(meta$Class)
table(meta$State)
table(meta$Type)
table(meta$Subtype)
length(table(meta$Type))      # 10
length(table(meta$Subtype))   # 17
table(meta$Cluster,meta$Type) 
table(meta$Type,meta$Subtype) 
# combined 41 clusters to 17 cell subtypes, and futher combined to 10 major cell types
                          
# gene expression matrix #2. downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132672
# use fread since it is too large
library(data.table)
gse <- fread("organoidatlas/GSE132672_allorganoids_withnew_matrix.txt.gz",header=T)
mat=gse
mat=as.data.frame(mat)
dim(mat) #[1]  16774 242350
mat[1:5,1:3]
genes = as.character(mat[,1])
mat[,1] = gsub(".+[|]", "", mat[,1])
mat = data.frame(mat[,-1], row.names=genes)
dim(mat) #[1]  16774 242349

mat[1:5,1:2]
              H1SWeek3_AAACCTGAGACAAAGG H1SWeek3_AAACCTGAGCACACAG
FO538757.2                     1.768941                         0
AP006222.2                     0.000000                         0
RP11-206L10.9                 19.036036                         0
RP11-54O7.1                    0.000000                         0
RP11-54O7.3                    0.000000                         0

which(colnames(mat) != rownames(meta)) # integer(0)

library(Seurat)
so <- CreateSeuratObject(counts = mat, project = "CorticalOrganoid", meta.data=meta)

save(so,file="organoidatlas/corticalorganoid.Robj")

m1=apply(so@assays$RNA@data,1,mean)
sd1=apply(so@assays$RNA@data,1,sd)
m2=apply(so@assays$RNA@data,2,mean)
sd2=apply(so@assays$RNA@data,2,sd)
png(file="organoidatlas/corticalorganoid_meansd.png",res=300,height=800,width=1500)
par(mfrow=c(1,2),mar=c(4,4,2,0.5),mgp=c(2,0.5,0))
plot(m1,sd1,pch=16,cex=.5,col=rgb(0,0,0,0.3),main="Per Gene",xlab="Mean",ylab="SD")
plot(m2,sd2,pch=16,cex=.5,col=rgb(0,0,0,0.3),main="Per Cell",xlab="Mean",ylab="SD")
dev.off()
              
dge=so
pdf(file=paste0("organoidatlas/corticalorganoid_PerCellAttributes_ViolinPlot.pdf"),height=8,width=6)
  plotlist=list()
Idents(dge) <- "Subtype"
plotlist[[1]]=VlnPlot(dge, features = "nFeature_RNA",pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
plotlist[[2]]=VlnPlot(dge, features = "nCount_RNA",log=T,pt.size=-1,cols=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+ theme(legend.position = 'none')
multiplot(plotlist,cols = 1)
dev.off()

Idents(dge) <- "Subtype"
centroidA=log(AverageExpression(dge)$RNA+1)
write.table(centroidA,file="organoidatlas/corticalorganoid_centroid_17celltypes.txt",row.names=T,col.names=T,quote=F,sep="\t")
# use this as reference cell type centroids
Idents(dge) <- "Cluster"
centroidC=log(AverageExpression(dge)$RNA+1)
write.table(centroidC,file="organoidatlas/corticalorganoid_centroid_41clusters.txt",row.names=T,col.names=T,quote=F,sep="\t")
# use this as reference cell type centroids

centroidA=read.table(file="organoidatlas/corticalorganoid_centroid_17celltypes.txt",sep="\t")
centroidC=read.table(file="organoidatlas/corticalorganoid_centroid_41clusters.txt")

ref=centroidA
# 1. using HVG from our data
genes=VariableFeatures(dgeall)
genes=genes[which(genes %in% rownames(ref))]
length(genes) # 1093

# 2. using HVG from 41 cluster centroids of ref2
g.mean <- apply(centroidC,1,mean)
g.var <- apply(centroidC,1,var)
g.var <- g.var[g.mean!=0]
g.mean <- g.mean[g.mean!=0]
summary(g.var/g.mean)
summary(g.mean)
hvg2 <- names(which(g.mean > 150 & g.var/g.mean > 30))
length(hvg2) # 3,289
genes=hvg2
genes=genes[which(genes %in% rownames(dgeall))]
length(genes) # 3,133


Idents(dgeall) <- dgeall$RNA_snn_res.0.2ordered
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

                          
                          
                          
                          
