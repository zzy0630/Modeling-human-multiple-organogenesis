# hEA gene correlation
# In R
library(DESeq2)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(RColorBrewer)
library(amap)
library(gplots)
library(dendextend)

# load('dds.vst.Rdata')
load('dds.deseq.Rdata')
load('dds.vst.Rdata')

normalized_count <- counts(dds.deseq, normalized=T)
normalized_count_mad <- apply(normalized_count, 1, mad)
normalized_count <- normalized_count[order(normalized_count_mad, decreasing = T),]

genes <- rownames(normalized_count)[0:3000]


org.rlog.exp <- assay(dds.vst)
org.rlog.exp <- org.rlog.exp[genes,]

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
spear_cor <- as.matrix(cor(org.rlog.exp, method = 'spearman'))
hc <- hcluster(t(org.rlog.exp), method = 'correlation')

dhc <- set(as.dendrogram(hc), 'branches_lwd',2)


### Correlation heatmap fig. S3
heatmap.2(spear_cor, Rowv = dhc, symm = T, trace = 'none', col=hmcol, margins = c(11,11))

tiff('./CorPlot/SpearmanCor.tiff', compression = 'lzw', units='in', width=8, height=8, pointsize = 8, res=600)
heatmap.2(spear_cor, Rowv = dhc, symm = T, trace = 'none', col=hmcol, margins = c(11,11))
dev.off()

tiff('./CorPlot/SpearmanCorNoKey.tiff', compression = 'lzw', units='in', width=8, height=8, pointsize = 8, res=600)
heatmap.2(spear_cor, Rowv = dhc, symm = T, trace = 'none', col=hmcol, margins = c(11,11), key='none')
dev.off()
