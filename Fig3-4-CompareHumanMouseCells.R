# Compare hOM cells with human fetal tissue and mouse embryo
#
# In R
library(biomaRt)
library(dplyr)
library(Seurat)
library(Matrix)
library(cellMapper)
library(org.Hs.eg.db)

# @ 
# Human fetal tissue
# https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/

# local # 
count <- readRDS('./DATA/gene_count_sampled.RDS') # SAMPLED CELLS
celldata <- readRDS('./DATA/df_cell.RDS') # CELL ANNOTATIONS
genedata <- readRDS('./DATA/df_gene.RDS') # GENE ANNOTATIONS

row.names(count) <- sub("[.][0-9]*","",row.names(count))
fetal <- CreateSeuratObject(count)
meta.data <- data.frame(celldata, row.names=celldata$sample)
fetal <- AddMetaData(fetal, meta.data)

ids <- select(org.Hs.eg.db, keys=rownames(fetal), columns = c('ENSEMBL','SYMBOL'), keytype='ENSEMBL')
'select()' returned 1:many mapping between keys and columns
# head(ids)
ids<-na.omit(ids)
ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
pos=match(ids$ENSEMBL,rownames(fetal))
fetal=fetal[pos,]
# Load function: RenameGenesSeurat
fetal<-RenameGenesSeurat(obj=fetal, newnames=ids$SYMBOL) # Human fetal data successfully converted!
saveRDS(fetal, "./DATA/fetal.processed.rds") # temporal saving

fetal.count <- GetAssayData(fetal, slot='data')
celldata <- readRDS('./DATA/df_cell.RDS')
meta.data <- data.frame(celldata, row.names=celldata$sample)
fetal <- CreateSeuratObject(counts=fetal.count)
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
fetal <- AddMetaData(fetal, metadata=meta.data)
saveRDS(fetal, "./DATA/fetal.processed.OK.rds") # human fetal tissue atlas

# Mouse organogenesis atlas
# downloaded from MOCA 

# Prepare dataset
cds <- readRDS('./DATA/cds_cleaned_sampled_100k.RDS') 
# cds_cleaned_sampled_100k.RDS downloaded from MOCA
ref.metadata <- read.table('./DATA/cell_annotate.csv', head=T, sep=',', stringsAsFactors=F) 
# cell_annotate.csv downloaded from MOCA

# Extract necessary data
counts<-cds@assayData$exprs
metadata<-cds@phenoData@data
gene.data<-cds@featureData@data
rownames(counts) <- as.character(gene.data$gene_short_name)

# Convert mous gene name to human 
mouse2human <- MouseHumanMapping(rownames(counts))
counts <- counts[names(mouse2human),]
rownames(counts) <- mouse2human[rownames(counts)]

# Create Seurat Obj
moca <- CreateSeuratObject(counts, project='MOCA', meta.data=metadata)
# Warning: Non-unique features (rownames) present in the input matrix, making unique
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
save(moca, file = './DATA/moca_seurat_mouse2human.RData') # temporal saving
###
###

moca <- NormalizeData(moca)

ref.metadata <- subset(ref.metadata, sample %in% colnames(moca))
ref.clusters <- factor(ref.metadata$Main_cell_type)
names(ref.clusters) <- ref.metadata$sample

ref.sub.traj <- ref.metadata$Sub_trajectory_name
names(ref.sub.traj) <- ref.metadata$sample

moca$Main_Cell_Type <- ref.clusters[colnames(moca)]
moca$Sub_Trajectory <- ref.sub.traj[colnames(moca)]
moca <- SetIdent(moca, value = "Main_Cell_Type")
table(moca$Main_Cell_Type)

save(moca, file='./DATA/moca_seurat_mouse2human.celltype.RData')

#####################################
# load hOM data
org <- readRDS('./organoid.mt20.umap.subclustered.rds')
marker <- read.csv('./Fetal-MOCA-marker.csv', head=F) 
# Fetal-MOCA-marker.csv integrates human and mouse marker, see data S6 of this study
genes.use <- marker$gene

#####################################
# fig. S10
#####################################
# Compare hOM with human
fetal <- readRDS('./DATA/fetal.processed.OK.rds')


VariableFeatures(fetal) <- genes.use
fetal <- ScaleData(fetal, features = genes.use)
org <- ScaleData(org, features = genes.use)

org <- SetIdent(org, value ='sub_cluster')
fetal<-SetIdent(fetal, value='Main_cluster_name')

correlation.res <- MapClustersCor(GetAssayData(org, slot = "scale.data"), org$sub_cluster,
                                  GetAssayData(fetal, slot = "scale.data"), fetal$Main_cluster_name,
genes.use = genes.use, metric = "pearson")


cluster.uniq <- unique(org$sub_cluster)

correlation.res.mat <- correlation.res$cluster.correlations

correlation.res.mat <- correlation.res.mat[,cluster.uniq]

save(correlation.res.mat, file='./DATA/org-subcluster.fetal.77cell.mat.Rdata')
# correlation matrix

#####################################
# fig. S11
#####################################
# Compare hOM with mouse
load('./DATA/moca_seurat_mouse2human.celltype.RData')

VariableFeatures(moca) <- genes.use
moca <- ScaleData(moca, features = genes.use)
save(moca, file='./DATA/moca_seurat_mouse2human.processed.RData')

org <- ScaleData(org, features = genes.use)

correlation.res <- MapClustersCor(GetAssayData(org, slot = "scale.data"), org$sub_cluster,
                                  GetAssayData(moca, slot = "scale.data"), moca$Main_Cell_Type,
genes.use = genes.use, metric = "pearson")

cluster.uniq <- unique(org$sub_cluster)


org.ref.correlations <- correlation.res$cluster.correlations

org.ref.correlations <- org.ref.correlations[,cluster.uniq]

save(org.ref.correlations, file='./DATA/org-mouse.correlation.matrix')
# correlation matrix

#############################################
# Define MapClustersCor
MapClustersCor <- function (query.data, query.clusters, train.data, train.clusters, 
    genes.use = NULL, metric = "pearson") 
{
    stopifnot(metric %in% c("pearson", "cosine"))
    if (is.null(genes.use)) 
        genes.use <- intersect(rownames(train.data), rownames(query.data))
    else genes.use <- genes.use[genes.use %in% rownames(query.data) & 
        genes.use %in% rownames(train.data)]
    query.cluster.data <- t(apply(query.data[genes.use, ], 1, 
        function(x) tapply(x, query.clusters, mean)))
    ref.cluster.data <- t(apply(train.data[genes.use, ], 1, function(x) tapply(x, 
        train.clusters, mean)))
    cluster.correlations <- correlate_cols(query.cluster.data, 
        ref.cluster.data, metric = metric)
    return(list(query.cluster.data = query.cluster.data, ref.cluster.data = ref.cluster.data, 
        cluster.correlations = cluster.correlations))
}

correlate_cols <- function(data.1, data.2, metric = "pearson") {
  cor.mat <- sapply(1:ncol(data.1), function(i) {
    sapply(1:ncol(data.2), function(j) {
      if (metric == "pearson") {
        cor(data.1[,i],data.2[,j])
      } else if (metric == "cosine") {
        lsa::cosine(data.1[,i],data.2[,j])
      }
    })
  })
  colnames(cor.mat) <- colnames(data.1)
  rownames(cor.mat) <- colnames(data.2)

  return(cor.mat)
}




































