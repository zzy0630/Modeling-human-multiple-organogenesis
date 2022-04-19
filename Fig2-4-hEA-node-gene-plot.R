# DESeq2 hEA sample RNA-seq
# In R
library(DESeq2)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(RColorBrewer)
library(reshape2)

ori <- read.table('./DESeq2_hEA/gene_count_matrix.csv', header=T, row.names=1, sep=',')

# Create variable table
samples <- c('Large-B','Large-B','Large-B','Large-BB','Large-BB','Large-BB','Large-C','Large-C','Large-C','Large-CB','Large-CB','Large-CB','Large-CC','Large-CC','Large-CC','Large-E8','Large-E8','Large-E8','Small-B','Small-B','Small-B','Small-BB','Small-BB','Small-BB','Small-C','Small-C','Small-C','Small-CB','Small-CB','Small-CB','Small-CC','Small-CC','Small-CC','Small-E8','Small-E8','Small-E8')

condition <- factor(samples)
coldata <- data.frame(row.names=colnames(org), condition)

# Create dds matrix
dds <- DESeqDataSetFromMatrix(countData=org, colData=coldata, design=~condition)

# Normalize/transform
# Using vst tranformation
dds.vst <- vst(dds)
save(dds.vst,file='dds.vst.Rdata')

# Get transformed data
org.vst.exp <- assay(dds.vst)
write.csv(org.vst.exp , file='org.vstTransform.exp.csv')

# PCA
pcaData <- plotPCA(dds.vst, intgroup=c('condition'), returnData=T)
write.csv(pcaData , file='org.PCAinfo.csv')

# DEseq
dds.deseq <- DESeq(dds)
save(dds.deseq, file='dds.deseq.Rdata')

###########################################

org <- read.table('./org.vstTransform.exp.csv', header=T, row.names=1, sep=',')


# Import node genes 
# Node genes are downloaded and selected from WikiPathway 
# See data S1 in this study
pluripotency <- read.csv('WikiPathway.pre-implantation.csv')
names(pluripotency) <- 'pluripotency'
ecto <- read.csv('WikiPathway.Ectoderm.csv')
names(ecto) <- 'ecto'
endo <- read.csv('WikiPathway.Endoderm.csv')
names(endo) <- 'endo'
meso <- read.csv('WikiPathway.Mesoderm.csv')
names(meso) <- 'meso'


ecto.only <- setdiff(setdiff(ecto$ecto, meso$meso), endo$endo)
meso.only <- setdiff(setdiff(meso$meso, ecto$ecto), endo$endo)
endo.only <- setdiff(setdiff(endo$endo, ecto$ecto), meso$meso)


stem.df <- data.frame(gene=pluripotency)
ecto.df <- data.frame(gene=ecto.only)
meso.df <- data.frame(gene=meso.only)
endo.df <- data.frame(gene=endo.only)

stem.df$germ <- 'stem'
ecto.df$germ <- 'ecto'
meso.df$germ <- 'meso'
endo.df$germ <- 'endo'

all.df <- rbind(stem.df,ecto.df, meso.df, endo.df)
rownames(all.df) <- paste(all.df$gene, all.df$gene, sep ='|')

samplename <- c('Large-B','Large-BB','Large-C','Large-CB','Large-CC','Large-E8','Small-B','Small-BB','Small-C','Small-CB','Small-CC','Small-E8')
suffix <- rep(c(1:3), length(samplename))
samplename <- rep(samplename, each=3)
samplename <- paste(samplename, suffix, sep = '.')
colnames(org)<-samplename

exp.mat <- org[rownames(all.df),]
sampleorder <- c('Small-E8','Small-B','Small-BB','Small-C','Small-CB','Small-CC','Large-E8','Large-B','Large-BB','Large-C','Large-CB','Large-CC')
suffix <- rep(c(1:3), length(sampleorder))
sampleorder <- rep(sampleorder, each=3)
sampleorder <- paste(sampleorder, suffix, sep = '.')
sampleorder<-as.character(sampleorder)


exp.mat <- exp.mat[,sampleorder]
exp.mat <- t(scale(t(exp.mat)))
exp.df <- na.omit(as.data.frame(exp.mat))

exp.df$gene <- rownames(exp.df)
exp.array <- melt(exp.df)
exp.array$germ <- all.df[exp.array$gene,2]
exp.zscore <- exp.array
exp.array$germ[which(exp.array$value < 0)] <- 'NA'
exp.array$variable<-gsub(pattern =".[[:digit:]]$", replacement = "", exp.array$variable)

x<-as.matrix(table(all.df$germ)*3)

# Make a integrated data frame
ESC <- '#FF99CC'
ecto <- '#0065B4'
meso <- '#9F0055'
endo <- '#FFCC00'
exp.zscore$variable<-gsub(pattern =".[[:digit:]]$", replacement = "", exp.zscore$variable)
exp.zscore$variable <- factor(exp.zscore$variable, levels=c('Small-E8','Small-B','Small-BB','Small-C','Small-CB','Small-CC','Large-E8','Large-B','Large-BB','Large-C','Large-CB','Large-CC'))
exp.zscore$germ <- factor(exp.zscore$germ, levels = c('stem','ecto','meso','endo'))

# Make boxplot # Fig. 2E
ggplot(data=exp.zscore)+geom_hline(yintercept = 0, lwd=0.8, color='black')+geom_boxplot(aes(x=variable,y=value, fill=germ), width=0.7, size=0.6)+theme_classic(base_size = 18)+scale_fill_manual(values = c(ESC,ecto,meso,endo))+scale_y_continuous(limits = c(-4,4), breaks=c(-4,-2,0,2,4))+guides(fill='none')
ggsave('devZscore.V2.tiff', dpi=300, width=8, height=3)




