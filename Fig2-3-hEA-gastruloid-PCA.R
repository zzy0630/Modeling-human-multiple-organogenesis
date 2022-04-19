# Merge hEA samples with mouse gastruloids
# PCA
# In R
library(DESeq2)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(RColorBrewer)
library(limma)

merge.ori <- read.table('./hEA_gastruloid.merge.csv', header = T, row.names = 1,sep=',')
# Load file 
# hEA_gastruloid.merge.csv is generated from Fig2-2-PrepDESeqMatrix.R

merge <- merge.ori[,c(1:9,13:27,31:50)] 
# Compare E8, Basal and hOM conditions with gastruloids

#---
org.samples <- c('Large-B','Large-B','Large-B','Large-BB','Large-BB','Large-BB','Large-C','Large-C','Large-C','Large-CC','Large-CC','Large-CC','Large-E8','Large-E8','Large-E8','Small-B','Small-B','Small-B','Small-BB','Small-BB','Small-BB','Small-C','Small-C','Small-C','Small-CC','Small-CC','Small-CC','Small-E8','Small-E8','Small-E8')
mouse.samples <- rep(c('H120','H144','H168','H24','H48','H72','H96'), each =2)
sample <- c(org.samples, mouse.samples)
condition <- factor(sample)
batch<-c(rep('A',30), rep('B',14))
batch <- factor(batch)
coldata <- data.frame(row.names=colnames(merge), condition, batch)
#---

dds <- DESeqDataSetFromMatrix(countData=merge, colData=coldata, design=~ condition)
# Normalize/transform
dds.vst <- vst(dds)
 assay(dds.vst) <- removeBatchEffect(assay(dds.vst), dds.vst$batch)
 plotPCA(dds.vst, intgroup=c('condition'),ntop=2000)

 pcaData <- plotPCA(dds.vst, intgroup=c('condition'),ntop=2000, returnData=T)
 write.csv(pcaData , file='org_gastruloid.PCAinfo.csv')

###############################
# Plot PCA # Fig. 2B
 #Set color
 ESC <- '#FFDDFF'
 D2C <-'#FF8B8B'
 D2B <- '#9DC3E6'
 D4B <- '#15629F'
 D4CC <- '#C00000'
 G24 <-'#E5F5E0'
 G48 <-'#C7E9C0'
 G72 <- '#A1D99B'
 G96 <-'#74C476'
 G120 <-'#41AB5D'
 G144 <-'#238B45'
 G168 <-'#006D2C'

 
 pca <- read.csv('org_gastruloid.PCAinfo.csv', head=T)
 pca$group <- factor(pca$group)
 
 
 ggplot(data=pca)+geom_vline(xintercept = 0, size=1.5, color='grey')+geom_hline(yintercept = 0,size=1.5, color='grey')+geom_point(aes(x=PC1, y=PC2, fill=group, shape=type), color='black', size=4, alpha=0.9)+theme_classic(base_size = 24)+scale_shape_manual(values=c(21,22,24))+scale_fill_manual(values=c(ESC,ESC,D2B,D2B,D2C,D2C,D4B,D4B,D4CC,D4CC,G120,G144,G168,G24,G48,G72,G96))+scale_x_continuous(limits = c(-80, 80),breaks = c(-60,-30,0,30,60))+scale_y_continuous(limits = c(-45,45),breaks=c(-40,-20,0,20,40))+guides(shape='none', fill='none')
 ggsave('org_gastruloid.PCA.tiff', dpi=300, width=4, height=3)
 
 ggplot(data=pca)+geom_vline(xintercept = 0, size=1.5, color='grey')+geom_hline(yintercept = 0,size=1.5, color='grey')+geom_point(aes(x=PC1, y=PC2, fill=group, shape=type), color='black', size=4, alpha=0.9)+theme_classic(base_size = 24)+scale_shape_manual(values=c(21,22,24))+scale_fill_manual(values=c(ESC,ESC,D2B,D2B,D2C,D2C,D4B,D4B,D4CC,D4CC,G120,G144,G168,G24,G48,G72,G96))+scale_x_continuous(limits = c(-80, 80),breaks = c(-60,-30,0,30,60))+scale_y_continuous(limits = c(-45,45),breaks=c(-40,-20,0,20,40))+guides(shape='none', fill='none')+geom_text(aes(x=PC1, y=PC2, label=group))
 ggsave('org_gastruloid.PCA.info.tiff', dpi=300, width=4, height=3)
 
 
 