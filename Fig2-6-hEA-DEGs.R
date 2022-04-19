# hEA DEGs
# In R
library(DESeq2)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(RColorBrewer)
library(fgsea)
library(dplyr)
library(tibble)


#load('dds.deseq.Rdata')
ori <- read.table('./DESeq2_hEA/gene_count_matrix.csv', header=T, row.names=1, sep=',')

# Create variable table
samplename <- c('Large-B','Large-BB','Large-C','Large-CB','Large-CC','Large-E8','Small-B','Small-BB','Small-C','Small-CB','Small-CC','Small-E8')
suffix <- rep(c(1:3), length(samplename))
samplename <- rep(samplename, each=3)
samplename <- paste(samplename, suffix, sep = '.')
colnames(org)<-samplename
rownames(org) <- gsub('\\|.*', '', rownames(org))

# Define function 
# compare hEA samples
# generate differentially expressed genes
makeDifftable <- function (a,b,c){
  samples <- c(a,b)
  suffix <- rep(c(1:3), length(samples))
  samples <- rep(samples, each=3)
  samples.suffix <- paste(samples, suffix, sep = '.')
  ori <- org[,samples.suffix]
  condition <- factor(samples)
  coldata <- data.frame(row.names=colnames(ori), condition)

  dds <- DESeqDataSetFromMatrix(countData=ori, colData=coldata, design=~condition)
  dds.deseq <- DESeq(dds)

  res <- results(dds.deseq)

  diff_gene_deseq2 <-subset(res,padj < 0.001 & (log2FoldChange >= 2 | log2FoldChange <= -2))
  diff_gene_deseq2_out <- cbind(row.names(diff_gene_deseq2),diff_gene_deseq2$baseMean,diff_gene_deseq2$log2FoldChange,
                              diff_gene_deseq2$lfcSE,diff_gene_deseq2$stat,diff_gene_deseq2$pvalue,diff_gene_deseq2$padj)
  colnames(diff_gene_deseq2_out)<-c('gene',colnames(diff_gene_deseq2))
  files <- paste(c,a,'-',b,'.csv',sep='')
  write.csv(diff_gene_deseq2_out, file=files)
}


# Compare Basal and hOM-treated samples with E8 samples 
# fig. S4A-C
lista <- c('Large-E8','Large-E8','Large-E8','Large-E8','Large-E8','Small-E8','Small-E8','Small-E8','Small-E8','Small-E8')
listb <- c('Large-B','Large-BB','Large-C','Large-CB','Large-CC','Small-B','Small-BB','Small-C','Small-CB','Small-CC')


for (i in c(1:length(lista))){
  makeDifftable(lista[i], listb[i],'./vsE8/')
}

# Compare hEA-CB with hEA-CC
# fig. S4D and E
lista <- c('Large-CB','Small-CB')
listb <- c('Large-CC','Small-CC')


for (i in c(1:length(lista))){
  makeDifftable(lista[i], listb[i],'./vsE8/')
}



# hEA size effect 
# fig. S6

lista <- c('Small-B','Small-BB','Small-C','Small-CB','Small-CC','Small-E8')
listb <- c('Large-B','Large-BB','Large-C','Large-CB','Large-CC','Large-E8')

for (i in c(1:length(lista))){
  makeDifftable(lista[i], listb[i],'./SizeEffect/')
}

