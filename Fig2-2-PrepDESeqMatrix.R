# Merge hEA DESeq matrix and gastruloid DESeq matrix
# in R
library(DESeq2)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(RColorBrewer)
library(limma)

org <- read.table('./DESeq2_hEA/gene_count_matrix.csv', header=T, row.names=1, sep=',')
mouse <- read.table('./DESeq2_gastruloid/gene_count_matrix.csv', header=T, row.names=1, sep=',')

### ### ###
# Convert gene name in mouse 
mouse.gene <- unlist(strsplit(rownames(mouse),'|', fixed=T))
mouse.all.gene <- mouse.gene[seq(from=1, to=length(mouse.gene), by=2)]
rownames(mouse) <- mouse.all.gene
mouse2human <- MouseHumanMapping(mouse.all.gene)
mouse <- mouse[names(mouse2human),]

org.gene <- unlist(strsplit(rownames(org),'|', fixed=T))
org.all.gene <- org.gene[seq(from=1, to=length(org.gene), by=2)]
rownames(org) <- org.all.gene

#------------------------------
rownames(mouse) <- toupper(rownames(mouse))
common <- intersect(rownames(org), rownames(mouse))
org <- org[common,]
mouse <- mouse[common,]
merge <- cbind (org, mouse)

#-------------------------------
write.csv(merge, 'hEA_gastruloid.merge.csv')

# 'Large-B','Large-B','Large-B','Large-BB','Large-BB','Large-BB','Large-C','Large-C','Large-C','Large-CB','Large-CB','Large-CB','Large-CC','Large-CC','Large-CC','Large-E8','Large-E8','Large-E8','Small-B','Small-B','Small-B','Small-BB','Small-BB','Small-BB','Small-C','Small-C','Small-C','Small-CB','Small-CB','Small-CB','Small-CC','Small-CC','Small-CC','Small-E8','Small-E8','Small-E8'
