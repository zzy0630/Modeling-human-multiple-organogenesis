# Heatmap shows marker gene expression
# In R
library(DESeq2)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(RColorBrewer)
library(reshape2)

org.ori <- read.table('./org.vstTransform.exp.csv', header=T, row.names=1, sep=',')
org <- org.ori[,c(1:9,13:27,31:36)] # keep E8, B, BB, C, CC samples

####################################
# Extract exp data for heatmap
# Pluripotency and developmental genes 
# fig. S5A
pluripotency <- c('POU5F1','NANOG','ZFP42','CDH1')
ecto<-c('SOX1','SOX3','OTX2','GBX2','PAX2','PAX6','CDH2')
mesendo <- c('GATA4','GATA6','EOMES')
meso <- c('TBXT','TBX6','MESP1','MESP2','HAND1','HAND2','FOXC1','FOXF1')
endo <- c('FOXA1','FOXA2','SOX17')
genes <- c(pluripotency, ecto, mesendo, meso, endo)
genes <- paste(genes, genes, sep ='|')

# Morphogens
# fig. S5B
wnts <- c('WNT1','WNT2','WNT3','WNT3A','WNT4','WNT5A','WNT5B','WNT6','WNT7A','WNT8A','WNT9A','WNT11','WNT16','DKK1')
bmps <- c('INHBA','BMP1','BMP2','BMP3','BMP4','BMP5','BMP6','BMP7','NODAL','NOG','CER1')
fgfs <- c('FGF2','FGF8','FGF12')
shh <- c('SHH')
morphogen <- c(wnts, bmps, fgfs, shh)
morphogen <- paste(morphogen, morphogen, sep ='|')

# Axial patterning
# fig. S5C
hox <- c('LHX2','SIX1','ONECUT2','LMX1B','PAX5','PAX7','PAX8','IRX3','MEOX1','HOXA1','HOXB1','HOXA2','HOXA3','HOXB4','PDX1','HOXC5','HOXB6','HOXC6','HOXB8','HOXC9','HOXB9','HOXD4','HOXB13','PAX3','MET','PAX1','NKX3-1')
hox <- paste(hox, hox, sep ='|')

# Extracellular matrix-related genes
# fig. S5D
ecm <- c('FN1','LAMA1','LAMB1','LAMB2','LAMC1','LAMA2','LAMC2','ITGA1','ITGA8','ITGA9','ITGA5','ITGA7','ITGA2','ITGA6','COL1A1','COL3A1','COL4A1','COL5A1','COL6A1','COL2A1','COL8A1','COL11A1','COL7A1','COL10A1','COL9A1')
ecm <- paste(ecm, ecm, sep ='|')

######################################
# Make heatmap

high<-'#FF333A'
low <- '#0512BD'
mid <- 'white'

samplename <- c('Large-B','Large-BB','Large-C','Large-CC','Large-E8','Small-B','Small-BB','Small-C','Small-CC','Small-E8')
suffix <- rep(c(1:3), length(samplename))
samplename <- rep(samplename, each=3)
samplename <- paste(samplename, suffix, sep = '.')
colnames(org)<-samplename

sampleorder <- c('Small-E8','Large-E8','Small-B','Small-BB','Large-B','Large-BB','Small-C','Small-CC','Large-C','Large-CC')
suffix <- rep(c(1:3), length(sampleorder))
sampleorder <- rep(sampleorder, each=3)
sampleorder <- paste(sampleorder, suffix, sep = '.')
sampleorder<-as.character(sampleorder)

exp.mat <- org[genes,]
exp.mat <- exp.mat[,sampleorder]
exp.mat <- t(scale(t(exp.mat)))
exp.df <- as.data.frame(exp.mat)

exp.df$gene <- rownames(exp.df)

exp.array <- melt(exp.df)

exp.array$variable<-factor(exp.array$variable, levels=sampleorder)
exp.array$gene <- factor(exp.array$gene, levels = rev(genes))


ggplot(data=exp.array)+geom_tile(aes(x=variable,y=gene, fill=value),color='white',size=1)+scale_fill_gradient2(low=low, high=high, mid=mid)+theme(axis.text.x = element_text(angle=90))+guides(fill='none')
ggsave('DevMarker.tiff', width=8, height=8,dpi=300) # fig. S5A

# HOX 
exp.mat <- org[hox,]
exp.mat <- exp.mat[,sampleorder]
exp.mat <- t(scale(t(exp.mat)))
exp.df <- as.data.frame(exp.mat)

exp.df$gene <- rownames(exp.df)

exp.array <- melt(exp.df)
exp.array$variable<-factor(exp.array$variable, levels=sampleorder)
exp.array$gene <- factor(exp.array$gene, levels = rev(hox))


ggplot(data=exp.array)+geom_tile(aes(x=variable,y=gene, fill=value),color='white',size=1)+scale_fill_gradient2(low=low, high=high, mid=mid)+theme(axis.text.x = element_text(angle=90))+guides(fill='none')
ggsave('BodyAxis.tiff', width=8, height=8,dpi=300) # fig. S5B

# Morphogen 
exp.mat <- org[morphogen,]
exp.mat <- exp.mat[,sampleorder]
exp.mat <- t(scale(t(exp.mat)))
exp.df <- as.data.frame(exp.mat)

exp.df$gene <- rownames(exp.df)

exp.array <- melt(exp.df)
exp.array$variable<-factor(exp.array$variable, levels=sampleorder)
exp.array$gene <- factor(exp.array$gene, levels = rev(morphogen))

ggplot(data=exp.array)+geom_tile(aes(x=variable,y=gene, fill=value),color='white',size=1)+scale_fill_gradient2(low=low, high=high, mid=mid)+theme(axis.text.x = element_text(angle=90))+guides(fill='none')
ggsave('Morphogen.tiff', width=8, height=8,dpi=300) # fig. S5C

# ECM 
exp.mat <- org[ecm,]
exp.mat <- exp.mat[,sampleorder]
exp.mat <- t(scale(t(exp.mat)))
exp.df <- as.data.frame(exp.mat)

exp.df$gene <- rownames(exp.df)

exp.array <- melt(exp.df)
exp.array$variable<-factor(exp.array$variable, levels=sampleorder)
exp.array$gene <- factor(exp.array$gene, levels = rev(ecm))


ggplot(data=exp.array)+geom_tile(aes(x=variable,y=gene, fill=value),color='white',size=1)+scale_fill_gradient2(low=low, high=high, mid=mid)+theme(axis.text.x = element_text(angle=90))+guides(fill='none')
ggsave('ECM.tiff', width=8, height=8,dpi=300) # fig. S5D


