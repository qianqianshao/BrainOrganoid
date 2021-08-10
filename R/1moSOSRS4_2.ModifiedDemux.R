# scRNA-seq analysis for 4 SOSRSs at 1-month old with HTO on 1.2.2021 by Qianyi
# fit bimodal or trimodal Guassian distributions to define cutoff thresholds for HTO classification


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

load(file="HTO20kcell_6kcell_hashtag.Robj")

## HTO normalization modified by Jun 
# normalized each HTO’s count vector by using the "positive" cluster - the cells that have been provisionally assigned to be positive for this HTO.  The mean and variance for this subset of cells provide a better representation than the global mean and variance.  
hto2=read.table("plot/hto_norm.txt",header=T,row.names=1) 
a=t(hto2)
png(file="plot/hto_4counts_6kcell_2.png",res=300,height=1000,width=2000)
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


## modify HTO demultiplexing

# to fit HTO_A and HTO_B in two modes:
hto2=read.table("plot/hto_norm.txt",header=T,row.names=1) 
tmp1<-data.frame(cbind(density(hto2[,1])$x,density(hto2[,1])$y))
fit1 <- nls(X2~ (C1 * exp(-(X1-m1)**2/(2 * s1**2)) +C2 * exp(-(X1-m2)**2/(2 * s2**2))), data= tmp1, start = list(C1 = 0.7, m1 = -4, s1 = 1, C2 = 0.3, m2 =0, s2 = 1.5), algorithm="port") 
#the initial values were chosen by making a guess based on the density plot
plot(tmp1[,2],pch=16,cex=0.8,ylab="density",main="HTO_A")
lines(predict(fit1),col=2)
 

tmp2<-data.frame(cbind(density(hto2[,2])$x,density(hto2[,2])$y))
fit2 <- nls(X2~ (C1 * exp(-(X1-m1)**2/(2 * s1**2)) +C2 * exp(-(X1-m2)**2/(2 * s2**2))), data= tmp2, start = list(C1 = 0.7, m1 = -4, s1 = 1, C2 = 0.3, m2 =0, s2 = 1.5), algorithm="port") 
plot(tmp2[,2],pch=16,cex=0.8,ylab="density",main="HTO_B")
lines(predict(fit2),col=2)


# Or, to fit HTO_C and HTO_D to three modes:
tmp3<-data.frame(cbind(density(hto2[,3])$x,density(hto2[,3])$y))
fit3 <- nls(X2~ (C1 * exp(-(X1-m1)**2/(2 * s1**2)) +C2 * exp(-(X1-m2)**2/(2 * s2**2)) + C3 * exp(-(X1-m3)**2/(2 * s3**2))), data= tmp3, start = list(C1 = 0.7, m1 = -4, s1 = 1, C2 = 0.3, m2 =-2, s2 = 1.5, C3 = 0.3, m3 =0, s3 = 1.5), algorithm="port")
plot(tmp3[,2],pch=16,cex=0.8,ylab="density",main="HTO_C")
lines(predict(fit3),col=3)

tmp4 <- data.frame(cbind(density(hto2[,4])$x,density(hto2[,4])$y))
fit4 <- nls(X2~ (C1 * exp(-(X1-m1)**2/(2 * s1**2)) +C2 * exp(-(X1-m2)**2/(2 * s2**2))), data= tmp4, start = list(C1 = 0.7, m1 = -4, s1 = 1, C2 = 0.3, m2 =0, s2 = 1.5), algorithm="port") 
plot(tmp4[,2],pch=16,cex=0.8,ylab="density",main="HTO_D")
lines(predict(fit4),col=3)

png(file="plot/hto_6kcell_2_cutoff.png",res=300,height=1000,width=1200)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(1.5, .5, 0))
plot(tmp1[,2],pch=16,cex=0.8,ylab="density",main="HTO_A")
lines(predict(fit1),col=2)
plot(tmp2[,2],pch=16,cex=0.8,ylab="density",main="HTO_B")
lines(predict(fit2),col=2)
plot(tmp3[,2],pch=16,cex=0.8,ylab="density",main="HTO_C")
lines(predict(fit3),col=3)
plot(tmp4[,2],pch=16,cex=0.8,ylab="density",main="HTO_D")
lines(predict(fit4),col=2)
dev.off()

 
summary(fit3)
#The coefficients are:
     Estimate Std. Error  t value Pr(>|t|)   
C1  0.1377209  0.0007693  179.013   <2e-16 ***
m1 -5.1598902  0.0064095 -805.033   <2e-16 ***
s1  0.6550047  0.0052745  124.184   <2e-16 ***
C2  0.2430829  0.0010008  242.899   <2e-16 ***
m2 -2.8550851  0.0065013 -439.157   <2e-16 ***
s2 -0.9301494  0.0084282 -110.362   <2e-16 ***
C3  0.0690268  0.0007286   94.733   <2e-16 ***
m3 -0.3082631  0.0348374   -8.849   <2e-16 ***
s3  1.1499430  0.0239912   47.932   <2e-16 ***
# Here the cutoff can be selected based on 
# the position of the 2nd and the 3rd mode: -2.855, -0.308, 
# and their sd estimates: 0.93, 1.15.
 
f1<-(hto[,1]>45)&(hto[,2]<95)&(hto[,4]<128) &(hto[,3]<160)
f2<-(hto[,1]<45)&(hto[,2]>95) &(hto[,4]<128) &(hto[,3]<160)
f4<-(hto[,1]<45)&(hto[,2]<95)&(hto[,4]>128) &(hto[,3]<160)
f3<-(hto[,1]<45)&(hto[,2]<95)&(hto[,4]<128) &(hto[,3]>160)

sum(f1)
sum(f2)
sum(f3)
sum(f4)



# Visualize demultiplexing results
labels=read.table("plot/HTO_labels.txt",header=T,row.names=1) 
table(labels)
Double  Empty  HTO-A  HTO-B  HTO-C  HTO-D 
  1885     97    507    626    701   2193 

label=as.character(labels[,1])
names(label)=rownames(labels)
label=label[names(Idents(dge.hashtag))]
label[which(label=="Double")] <- "Doublet"
label[which(label=="Empty")] <- "Negative"

label2=label
label2[grepl("HTO",label)] <- "Singlet"

dge.hashtag$label <- label
dge.hashtag$label <- factor(dge.hashtag$label,levels=c("Doublet","HTO-A","HTO-B","HTO-C","HTO-D","Negative"))
dge.hashtag$label2 <- label2
dge.hashtag$label2 <- factor(dge.hashtag$label2,levels=c("Doublet","Negative","Singlet"))

Idents(dge.hashtag) <- "label"
p=list()
p[[1]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")
p[[2]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-A", feature2 = "hto_HTO-C")
p[[3]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-A", feature2 = "hto_HTO-D")
p[[4]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-B", feature2 = "hto_HTO-C")
p[[5]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-B", feature2 = "hto_HTO-D")
p[[6]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = "hto_HTO-C", feature2 = "hto_HTO-D")
png(file="plot/HTO_Classify2_pairs.png",res=300,height=2000,width=3600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(p,cols=3)
dev.off()

table(dge.hashtag$label,dge.hashtag$ID)
           Doublet HTO-A HTO-B HTO-C HTO-D Negative
  Doublet     1297   166   171   211    40        0
  HTO-A        143   364     0     0     0        0
  HTO-B         69     0   557     0     0        0
  HTO-C        148     0     0   553     0        0
  HTO-D        233     5    61    29   863     1002
  Negative      13     3    17    12     0       52


save(dge.hashtag,file="HTO20kcell_6kcell_hashtag.Robj")

p=list()
p[[1]]=FeatureScatter(dge, pt.size=.8,feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")
p[[2]]=FeatureScatter(dge, pt.size=.8,feature1 = "hto_HTO-A", feature2 = "hto_HTO-C")
p[[3]]=FeatureScatter(dge, pt.size=.8,feature1 = "hto_HTO-A", feature2 = "hto_HTO-D")
p[[4]]=FeatureScatter(dge, pt.size=.8,feature1 = "hto_HTO-B", feature2 = "hto_HTO-C")
p[[5]]=FeatureScatter(dge, pt.size=.8,feature1 = "hto_HTO-B", feature2 = "hto_HTO-D")
p[[6]]=FeatureScatter(dge, pt.size=.8,feature1 = "hto_HTO-C", feature2 = "hto_HTO-D")
png(file="plot/HTO_2Classify_pairs.png",res=300,height=2000,width=3600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(p,cols=3)
dev.off()



# Compare number of UMIs for singlets, doublets and negative cells
Idents(dge.hashtag) <- "label2"
pdf(file="plot/HTO_Classification2.pdf",height=4,width=5)
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
png(file="plot/HTO_Classify2_tSNE.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.hashtag.subset)
dev.off()
png(file="plot/HTO_Classify2_tSNE2.png",res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.hashtag.subset,group="label")
dev.off()

# You can also visualize the more detailed classification result by running Idents(object) <-
# 'HTO_classification' before plotting. Here, you can see that each of the small clouds on the
# tSNE plot corresponds to one of the 28 possible doublet combinations.
# Create an HTO heatmap, based on Figure 1C in the Cell Hashing paper.

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
pdf(file="plot/HTO_Classification2_heatmap.pdf",height=4,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
HTOHeatmap(dge.hashtag, classification="label",global.classification="label2",assay = "HTO", ncells = nrow(dge.hashtag@meta.data))
dev.off()

