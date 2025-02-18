# scRNA-seq analysis for 6 SOSRSs at 3-month old with HTOs on 2.4.2021 by Qianyi
# 6 replicates sequenced in 2 samples, each with 3 SOSRS (3 HTO tags)
# fit bimodal or trimodal Guassian distributions to define cutoff thresholds for HTO classification in each sample

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

# load all data before removing negatives and doublets
dgelist=list()
for(i in 1:2){
  load(file=paste0(gexset[i],"_hashtag.Robj"))
  dgelist[[i]]=dge.hashtag
}


## modify HTO demultiplexing

i=1;
cutmin=NULL;            # lowestest point on the data between the two peaks
cutoverlap=c(92,97,115); # intersect point of two gaussian curves
cut1data=cut2data=NULL;

i=2;
cutmin=NULL;            # lowestest point on the data between the two peaks
cutoverlap=c(95,108,107); # intersect point of two gaussian curves
cut1data=cut2data=NULL;

# do i=1 with the following; then repeat for i=2



counts=t(as.matrix(dge.hashtag@assays$HTO@counts))
hto=t(as.matrix(dge.hashtag@assays$HTO@data))
tmp1<-data.frame(cbind(density(hto[,1])$x,density(hto[,1])$y))
tmp2<-data.frame(cbind(density(hto[,2])$x,density(hto[,2])$y))
tmp3<-data.frame(cbind(density(hto[,3])$x,density(hto[,3])$y))
tmplist=list(tmp1,tmp2,tmp3)



for(j in 1:3){

gamma1 <- tmplist[[j]]

names(gamma1) <- c("x","y")
plot(y~x, data = gamma1, type = "l")
library(quantmod)
# Find extreme
gamma.max <- findPeaks(gamma1$y)[1:2] #rejects a number of extraneous peaks
abline(v = gamma1$x[gamma.max])
ep.min <- findValleys(gamma1$y)
plot(y~x, data = gamma1, type = "l", main = "Minima")
abline(v = gamma1$x[ep.min], col = "blue", lty = 2)
print(ep.min)
print(gamma1$x[ep.min])
cutmin=c(cutmin,ep.min[1])
# the peak looks more accurate
#Estimate the Ci
peak.heights <- gamma1$y[gamma.max]
#Estimate the mu_i
gamma.mu <- gamma1$x[gamma.max] #the same values as before
#Estimate the sigma_i from the FWHM
FWHM.finder <- function(ep.data, mu.index){
  peak.height <- ep.data$y[mu.index]
  fxn.for.roots <- ep.data$y - peak.height/2
  indices <- 1:nrow(ep.data)
  root.indices <- which(diff(sign(fxn.for.roots))!=0)
  tmp <- c(root.indices,mu.index) %>% sort
  tmp2 <- which(tmp == mu.index)
  first.root <- root.indices[tmp2 -1]
  second.root <- root.indices[tmp2]
  HWHM1 <- ep.data$x[mu.index] - ep.data$x[first.root]
  HWHM2 <- ep.data$x[second.root] - ep.data$x[mu.index]
  FWHM <- HWHM2 + HWHM1
  FWHM2 = 2*min(c(HWHM1,HWHM2))
  return(list(HWHM1 = HWHM1,HWHM2 = HWHM2,FWHM = FWHM,FWHM2 = FWHM2))
}
FWHM <- lapply(gamma.max, FWHM.finder, ep.data = gamma1)
gamma.sigma <- unlist(sapply(FWHM, '[', 'FWHM2'))/2.355
print(rbind(peak.heights,gamma.mu,gamma.sigma))

#Perform the fit
fit <- nls(y ~ (C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                  C2*exp(-(x-mean2)**2/(2 * sigma2**2))),
           data = gamma1,
           start = list(mean1 = gamma.mu[1],
                        mean2 = gamma.mu[2],
                        sigma1 = gamma.sigma[1],
                        sigma2 = gamma.sigma[2],
                        C1 = peak.heights[1],
                        C2 = peak.heights[2]),
           algorithm = "port")

#Plot the fit
dffit <- data.frame(x=seq(0, 1 , 0.001))
dffit$y <- predict(fit, newdata=dffit)
fit.sum <- summary(fit)
coef.fit <- fit.sum$coefficients[,1]
mu.fit <- coef.fit[1:2]
sigma.fit <- coef.fit[3:4]
C.fit <- coef.fit[5:6]



png(file=paste0("plot/HTO_",htoset[[i]][j],"_6kcell_cutoff.png"),res=300,height=1200,width=1200)
par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.5, .5, 0))

plot(y ~ x, data = gamma1, type = "l",main=paste0("HTO_",htoset[[i]][j]))
legend("topright", lty = c(1,1,1), col = c("black", "green", "blue","red"), c("Data", "Monoclonal", "Gamma", "Sum"))
lines(y ~ x, data = dffit, col ="red", cex = 0.2)
for(ii in 1:2){
  x <- dffit$x
  y <- C.fit[ii] *exp(-(x-mu.fit[ii])**2/(2 * sigma.fit[ii]**2))
  lines(x,y, col = ii + 2)
}
abline(v=gamma1$x[ep.min[1]],col=5)
# cutoverlap[j]=95
abline(v=gamma1$x[cutoverlap[j]],col=5)

plot(gamma1$y,pch=16,cex=0.8,ylab="density",main=paste0("HTO_",htoset[[i]][j]))
lines(predict(fit),col=2)
abline(v=ep.min[1],col=5)
abline(v=cutoverlap[j],col=5)

dev.off()

summary(fit)

# cutoff for normalized data
cut2data=c(cut1data,gamma1$x[cutoverlap[j]])
cut1data=c(cut1data,gamma1$x[ep.min[1]])

# calculate cutoffs for raw HTO counts
cut2=unique(hto[which(round(hto[,j],2)==round(gamma1$x[cutoverlap[j]],2)),j])
cut2cell=names(hto[which(round(hto[,j],2)==round(gamma1$x[cutoverlap[j]],2)),j])
if(length(cut2)==0){
 cut2=unique(hto[which(round(hto[,j],2)==round(ceiling(gamma1$x[cutoverlap[j]]),2)),j])
 cut2cell=names(hto[which(round(hto[,j],2)==round(ceiling(gamma1$x[cutoverlap[j]]),2)),j])
}
if(length(cut2)==0){
 cut2=unique(hto[which(round(hto[,j],2)==round(floor(gamma1$x[cutoverlap[j]]),2)),j])
 cut2cell=names(hto[which(round(hto[,j],2)==round(floor(gamma1$x[cutoverlap[j]]),2)),j])
}
cut2value=unique(counts[cut2cell,j])

cut1=unique(hto[which(round(hto[,j],2)==round(gamma1$x[ep.min[1]],2)),j])
cut1cell=names(hto[which(round(hto[,j],2)==round(gamma1$x[ep.min[1]],2)),j])
if(length(cut1)==0){
  cut1=unique(hto[which(round(hto[,j],2)==round(ceiling(gamma1$x[ep.min[1]]),2)),j])
  cut1cell=names(hto[which(round(hto[,j],2)==round(ceiling(gamma1$x[ep.min[1]]),2)),j])
}
if(length(cut1)==0){
  cut1=unique(hto[which(round(hto[,j],2)==round(floor(gamma1$x[ep.min[1]]),2)),j])
  cut1cell=names(hto[which(round(hto[,j],2)==round(floor(gamma1$x[ep.min[1]]),2)),j])
}
cut1value=unique(counts[cut1cell,j])

print(c(cutoverlap[j],gamma1$x[cutoverlap[j]],cut2value))
print(c(ep.min[1],gamma1$x[ep.min[1]],cut1value))
}


### classify cells based on the cutoffs
s1<-rownames(hto)[(hto[,1]>=cut2data[1])&(hto[,2]<cut2data[2])&(hto[,3]<cut2data[3])]
s2<-rownames(hto)[(hto[,1]<cut2data[1])&(hto[,2]>=cut2data[2])&(hto[,3]<cut2data[3])]
s3<-rownames(hto)[(hto[,1]<cut2data[1])&(hto[,2]<cut2data[2])&(hto[,3]>=cut2data[3])]

d12<-rownames(hto)[(hto[,1]>=cut2data[1])&(hto[,2]>=cut2data[2])&(hto[,3]<cut2data[3])]
d13<-rownames(hto)[(hto[,1]>=cut2data[1])&(hto[,2]<cut2data[2])&(hto[,3]>=cut2data[3])]
d23<-rownames(hto)[(hto[,1]<cut2data[1])&(hto[,2]>=cut2data[2])&(hto[,3]>=cut2data[3])]

d123<-rownames(hto)[(hto[,1]>=cut2data[1])&(hto[,2]>=cut2data[2])&(hto[,3]>=cut2data[3])]

neg<-rownames(hto)[(hto[,1]<cut2data[1])&(hto[,2]<cut2data[2])&(hto[,3]<cut2data[3])]

c(length(s1),length(s2),length(s3),length(d12),length(d13),length(d23),length(d123),length(neg))

label=c(rep(paste0("HTO-",htoset[[i]][1]),length(s1)),
rep(paste0("HTO-",htoset[[i]][2]),length(s2)),
rep(paste0("HTO-",htoset[[i]][3]),length(s3)),
rep(paste(paste0("HTO-",htoset[[i]][1:2]),collapse="_"),length(d12)),
rep(paste(paste0("HTO-",htoset[[i]][c(1,3)]),collapse="_"),length(d13)),
rep(paste(paste0("HTO-",htoset[[i]][2:3]),collapse="_"),length(d23)),
rep(paste(paste0("HTO-",htoset[[i]][1:3]),collapse="_"),length(d123)),
rep("Negative",length(neg))
  )
names(label)=c(s1,s2,s3,d12,d13,d23,d123,neg)
write.table(label,paste0("plot/hto_",i,"labels.txt"),row.names=T,col.names=F,quote=F,sep="\t")

## HTO classified by me by fitting two Gaussian curves
labels=read.table(paste0("plot/hto_",i,"labels.txt"),row.names=1) 
table(labels)

label=as.character(labels[,1])
names(label)=rownames(labels)
label=label[names(Idents(dge.hashtag))]
dge.hashtag$label1 <- label

label1=label
label[grepl("_",label)] <- "Doublet"
label2=label
label2[grepl("HTO",label)] <- "Singlet"

dge.hashtag$label <- label
dge.hashtag$label <- factor(dge.hashtag$label,levels=c("Doublet",paste0("HTO-",htoset[[i]]),"Negative"))
dge.hashtag$label2 <- label2
dge.hashtag$label2 <- factor(dge.hashtag$label2,levels=c("Doublet","Negative","Singlet"))

table(label)
table(label1)
table(label2)

Idents(dge.hashtag) <- "label"
save(dge.hashtag,file=paste0(gexset[i],"_hashtag.Robj"))


# Visualize demultiplexing results
Idents(dge.hashtag) <- "label"
p=list()
p[[1]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = paste0("hto_HTO-",htoset[[i]][1]), feature2 = paste0("hto_HTO-",htoset[[i]][2]))
p[[2]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = paste0("hto_HTO-",htoset[[i]][1]), feature2 = paste0("hto_HTO-",htoset[[i]][3]))
p[[3]]=FeatureScatter(dge.hashtag, pt.size=.8,feature1 = paste0("hto_HTO-",htoset[[i]][2]), feature2 = paste0("hto_HTO-",htoset[[i]][3]))
png(file=paste0("plot/hto_",i,"FitClassify_pairs.png"),res=300,height=1000,width=3600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
multiplot(p,cols=3)
dev.off()

Idents(dge.hashtag) <- "label"
feature1=paste0("hto_HTO-",htoset[[i]][1])
feature2=paste0("hto_HTO-",htoset[[i]][2])
feature3=paste0("hto_HTO-",htoset[[i]][3])
data = FetchData(object = dge.hashtag, vars = c(feature1, feature2, feature3))
ids=as.numeric(dge.hashtag$label)
ids1=ids2=ids3=ids
ids1[which(ids>3)]<-1
ids2[which(ids==3 | ids==5)]<-1
ids2[which(ids==4)]<-3
ids3[which(ids==2 | ids==5)]<-1
ids3[which(ids==3)]<-2
ids3[which(ids==4)]<-3
png(file=paste0("plot/hto_",i,"FitClassify_pairs2.png"),res=300,height=1000,width=3000)
par(mfrow=c(1,3),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(data[,1],data[,2],cex.lab=1.5,xlab=colnames(data)[1],ylab=colnames(data)[2],pch=16,col=ids1,cex=.5)
plot(data[,1],data[,3],cex.lab=1.5,xlab=colnames(data)[1],ylab=colnames(data)[3],pch=16,col=ids2,cex=.5)
plot(data[,2],data[,3],cex.lab=1.5,xlab=colnames(data)[2],ylab=colnames(data)[3],pch=16,col=ids3,cex=.5)
dev.off()

# zoomed-in view, order the singlets to be on top
Idents(dge.hashtag) <- "label"
ids=as.numeric(dge.hashtag$label)
names(ids)=names(dge.hashtag$label)
table(dge.hashtag$label2)
order=names(sort(dge.hashtag$label2))

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

table(dge.hashtag$label,dge.hashtag$ID)

# Compare number of UMIs for singlets, doublets and negative cells
Idents(dge.hashtag) <- "label2"
pdf(file=paste0("plot/hto_",i,"_GaussianFitClassification.pdf"),height=4,width=5)
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
# remove negative and visualize classification by fitting gaussians
png(file=paste0("plot/hto_",i,"FitClassify_PCA2.png"),res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.hashtag.subset,reduction="pca",group="label")
dev.off()
png(file=paste0("plot/hto_",i,"FitClassify_tSNE2.png"),res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.hashtag.subset,group="label")
dev.off()

pdf(file=paste0("plot/hto_",i,"FitClassification_heatmap.pdf"),height=4,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
HTOHeatmap(dge.hashtag, classification="label",global.classification="label2",assay = "HTO", ncells = nrow(dge.hashtag@meta.data))
dev.off()

# Cluster and visualize cells using the usual scRNA-seq workflow, and examine for the potential presence of batch effects.
dgelist2=list()
dge.hashtag=dgelist[[i]]

# Extract the singlets
Idents(dge.hashtag) <- "label2"
dge.singlet <- subset(dge.hashtag, idents = "Singlet")

# Select the top 1000 most variable features
dge.singlet <- FindVariableFeatures(dge.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
dge.singlet <- ScaleData(dge.singlet, features = VariableFeatures(dge.singlet))

# Run PCA
dge.singlet <- RunPCA(dge.singlet, features = VariableFeatures(dge.singlet))

dgelist2[[i]]=dge.singlet
# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
dge.singlet <- FindNeighbors(dge.singlet, reduction = "pca", dims = 1:10)
dge.singlet <- FindClusters(dge.singlet, resolution = seq(0.1,1,0.1), verbose = FALSE)
dge.singlet <- RunTSNE(dge.singlet, reduction = "pca", dims = 1:10)
dge.singlet <- RunUMAP(dge.singlet, reduction = "pca", dims = 1:10)

png(file=paste0("plot/hto_",i,"FitClassify_singlet_tSNE.png"),res=300,height=1000,width=1200)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge.singlet, reduction="tsne",group.by = "label")
dev.off()

dgelist2[[i]]=dge.singlet
save(dge.singlet,file=paste0(gexset[i],"_singlet.Robj"))



