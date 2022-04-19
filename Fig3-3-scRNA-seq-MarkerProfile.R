# UMAP marker profile #
library(dplyr)
library(magrittr)
library(Seurat)
#
library(ggplot2)
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(FlexDotPlot)


#################################
# Main clusters # Fig. 3C
organoid <- readRDS('./organoid.mt20.umap.rds') # Read in

col.low <- '#d5dfff'
col.high <- '#cc0000'

clusters <- c(2,3,10,4,6,8,9,14,15,16,12,19) # Ectodermal lineage
cells <- names(Idents(subset(organoid, idents=clusters))) # Defining subset of cells

gene <- 'CCNO'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'DCX'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'RAX'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'OTX2'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'EMX2'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'FGF17'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()
gene <- 'OLIG3'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'HOXA2'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'GDF7'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'SOX10'
file <- paste('./Fig3/Seurat-mt20.ecto.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

################################
clusters <- c(0,1,5,7,11,13,17,18) # Mesendodermal lineage
cells <- names(Idents(subset(organoid, idents=clusters))) # Defining subset of cells

gene <- 'COL3A1'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'FOXA2'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'FOXF1'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'FOXC1'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'NKX2-5'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'SIX1'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'BARX1'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'PECAM1'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'MYH6'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'TBX1'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'ISL1'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'PRRX1'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

gene <- 'SHH'
file <- paste('./Fig3/Seurat-mt20.mesendo.',gene,'.tiff', sep='')
tiff(file, units='in', width=5, height=5, res=300, compression ='lzw')
FeaturePlot(object=organoid, features=gene, cols=c(col.low,col.high), order=T,cells=cells,pt.size=0.5)
dev.off()

################################
# Marker profile in sub-clusters

# Ectoderm # fig. S8C
org.Ecto.umap <- readRDS('./org.Ecto.rds')

gene <- c('PAX6','SOX2','DCX','RBFOX3','RAX','VSX2','OTX2','FEZF1','EMX2','FGF17','GBX2','OLIG3','BMP5','GDF7','FOXD3','SOX10','DCT','MKI67','FOXG1','PHOX2B','ISL1','CCNO','GDF15')

for (i in gene) {

file <- paste('./FigS9/Ecto.',i,'.tiff', sep='')
tiff(file, units='in', width=4, height=4, res=300, compression ='lzw')
p <- FeaturePlot(object=org.Ecto.umap, features=i, cols=c('pink','black'), order=T,pt.size=0.5)
print(p)
dev.off()
}

# Mesoderm # fig. S9C
org.Meso.umap <- readRDS('./org.Mesendo.rds')

gene <- c('TBX1','SIX1','LHX9','WT1','MYH6','MKI67','PRRX1','ISL1', 'GATA4','NKX2-5','HAND1','SOX2','CD70','IRX4','PAX2','PAX8','CDH5','BARX1','FOXD1','CA9','FOXC1','FOXC2','CALB2','SMYD1','PAX7','ZIC2','MEOX1','MEOX2','COX6A2','TBX3')

for (i in gene) {

file <- paste('./FigS9/Meso.',i,'.tiff', sep='')
tiff(file, units='in', width=4.5, height=4, res=300, compression ='lzw')
p <- FeaturePlot(object=org.Meso.umap, features=i, cols=c('pink','black'), order=T,pt.size=0.5)
print(p)
dev.off()
}

# Endoderm # fig. S9D
org.Endo.umap <- readRDS('./org.Endo.rds')

gene <- c('CTNNB1','EPCAM','OTX2','PAX1','SIX3','OSR1','TTR','AFP','HOXB4','HOXA5','HOXC5','FOXA2','SOX9','LGR5','MKI67','GDF15','PROX1','SHH','KLF5','TPPP3','CRLF1')

for (i in gene) {

file <- paste('./FigS9/Endo.',i,'.tiff', sep='')
tiff(file, units='in', width=3, height=2.5, res=300, compression ='lzw')
p <- FeaturePlot(object=org.Endo.umap, features=i, cols=c('pink','black'), order=T,pt.size=0.5)
print(p)
dev.off()
}








