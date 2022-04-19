# Compare day-20 hOM with temporal ordered mouse embryos
# Fig. 4B
# MOCA: E9.5 E10.5 E11.5 E12.5 E13.5
# Early mouse embryo-gastrulation: E6.5-E8.5
# Early mouse embryo-gastrulation data were downloaded from: 
# curl https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz > atlas_data.tar.gz
# In R
library(Seurat)
library(cellMapper)
library(scran)

organoid <- readRDS('./organoid.mt20.umap.rds') # hOM cells
load('./DATA/moca_seurat_mouse2human.processed.RData') # E9.5-13.5 data

################################
# Prepare mouse early embryo data
# E6.5-8.5 data

mouse.data<-Read10X('./') # atlas data
metadata<-read.table('meta.tab', sep='\t', head=T)
rownames(metadata)<-metadata$cell

# Convert to human gene 
mouse2human <- MouseHumanMapping(rownames(mouse.data))
head(mouse2human)
#    Xkr4      Rp1    Sox17   Mrpl15   Lypla1    Tcea1 
#   "XKR4"    "RP1"  "SOX17" "MRPL15" "LYPLA1"  "TCEA1" 
mouse.data <- mouse.data[names(mouse2human),]
rownames(mouse.data) <- mouse2human[rownames(mouse.data)]
#

mouse <- CreateSeuratObject(counts = mouse.data, project = 'mouse', min.cells=3, min.features=200)
mouse <- AddMetaData(object=mouse, metadata=metadata)
# saveRDS(mouse, './gastrulation.raw.rds') # make backup
# mouse <- readRDS('./gastrulation.raw.rds') # read in

unique(mouse@meta.data$stage)
# [1] "E6.5"               "E7.5"               "E6.75"             
# [4] "E7.75"              "E7.0"               "E8.0"              
# [7] "E8.5"               "mixed_gastrulation" "E7.25"             
# [10] "E8.25"  
# no need to use them all, use only:
mouse <- subset(mouse, cells=rownames(subset(mouse@meta.data, stage %in% c('mixed_gastrulation', 'E6.5','E7.0','E7.5','E8.0','E8.5'))))

# 
mouse[['percent.mt']] <- PercentageFeatureSet(mouse, pattern='^MT-')
mouse <- subset(mouse, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 20)
# An object of class Seurat 
# 16831 features across 85506 samples within 1 assay 
# Active assay: RNA (16831 features, 0 variable features)
 mouse.list <- SplitObject(mouse, split.by='stage')
summary(mouse.list)
#      Length Class  Mode
# mixed_gastrulation 1      Seurat S4 
# E6.5 1      Seurat S4  
# E7.5 1      Seurat S4  
# E7.0 1      Seurat S4  
# E8.0 1      Seurat S4  
# E8.5 1      Seurat S4 

################################
# Make temporal mouse embryos

moca@meta.data$stage <- paste('E',moca@meta.data$day, sep='')
moca@meta.data$stage<-factor(moca@meta.data$stage, levels=c('E9.5','E10.5','E11.5','E12.5','E13.5'))
moca@meta.data$celltype <- moca@meta.data$Main_Cell_Type
moca.list <- SplitObject(moca, split.by='stage')
summary(moca.list)
#      Length Class  Mode
# E12.5 1      Seurat S4  
# E13.5 1      Seurat S4  
# E11.5 1      Seurat S4  
# E10.5 1      Seurat S4  
# E9.5  1      Seurat S4  


mouse.total.list <- SplitObject(mouse, split.by='stage')
moca.total.list <- SplitObject(moca, split.by='stage')

all.total.list<-list()
all.total.list <- mouse.total.list
for (i in names(moca.total.list)){
#  print(all.list[[i]])
	all.total.list[[i]] <- moca.total.list[[i]]
}

organoid@meta.data$celltype<-organoid@meta.data$seurat_clusters
organoid$stage <- 'organoid'
all.total.list[['organoid']] <- organoid

summary(all.total.list)


################################
# Ectoderm cosine similarity
# Ectoderm @ mouse
mouse.ecto <- subset(mouse, cells=rownames(subset(mouse@meta.data, celltype %in% c('Surface ectoderm','Spinal cord','Rostral neurectoderm','Caudal neurectoderm','NMP','Forebrain/Midbrain/Hindbrain','Neural crest','Primitive Streak'))))

# Ectoderm @ moca
moca.ecto <- subset(moca, cells=rownames(subset(moca@meta.data, celltype %in% c('Neural progenitor cells','Radial glia','Excitatory neurons','Inhibitory neurons','Postmitotic premature neurons','Inhibitory neuron progenitors','Neural Tube'))))

# Ectoderm @ organoid
organoid.ecto <- subset(organoid, cells=rownames(subset(organoid@meta.data, celltype %in% c(2,3,4,6,9,12,14))))

mouse.ecto.list <- SplitObject(mouse.ecto, split.by='stage')
moca.ecto.list <- SplitObject(moca.ecto, split.by='stage')

all.ecto.list <- mouse.ecto.list
for (i in names(moca.ecto.list)){
	all.ecto.list[[i]] <- moca.ecto.list[[i]]
}
all.ecto.list[['organoid']] <- organoid.ecto

summary(all.ecto.list)

load.list <- list()
emb.list <-list()
mouse.exp.list <- list()
organoid.compare.score.list <- list()
var.gene<-c()


for (i in names(all.ecto.list)){
	all.ecto.list[[i]] <- ScaleData(all.ecto.list[[i]], verbose=F)
	varf <- FindVariableFeatures(all.ecto.list[[i]], selection.method = "vst", nfeatures = 2000, verbose=F)
	gene.use <- VariableFeatures(varf)
	all.ecto.list[[i]] <- RunPCA(all.ecto.list[[i]], features=gene.use, npcs=20, verbose=F)
	load.list[[i]] <- Loadings(all.ecto.list[[i]], reduction='pca')
	emb.list[[i]] <- Embeddings(all.ecto.list[[i]], reduction='pca')
	stage.factor <- all.ecto.list[[i]]$stage
 	if (i %in% c('E9.5','E10.5','E11.5','E12.5','E13.5')){
 		stage.factor <- droplevels(stage.factor)
 	}

	mouse.exp.list[[i]] <- t(apply(t(emb.list[[i]]),1,function(x) tapply(x,stage.factor,mean)))
	
	print (i)
	
}

ecto.array <- data.frame()

for (i in names(mouse.exp.list)){
	if (i %in% c('E6.5','E7.5','E8.5','E9.5','E10.5','E11.5','E12.5','E13.5','organoid')){
	A <- mouse.exp.list[[i]]
	for (k in names(mouse.exp.list)){
#	print(k)
	if (k %in% c('E6.5','E7.5','E8.5','E9.5','E10.5','E11.5','E12.5','E13.5','organoid')){
		B <- mouse.exp.list[[k]]
		ecto.array[i,k] <- as.numeric(lsa::cosine(A[,],B[,]))
	}}}
}

write.csv(ecto.array, file='./DATA/ecto.total.CS.array.csv')

################################
# Mesoderm cosine similarity
# mesoderm @ mouse
mouse.meso <- subset(mouse, cells=rownames(subset(mouse@meta.data, celltype %in% c('Primitive Streak','Nascent mesoderm','Pharyngeal mesoderm','Paraxial mesoderm','Caudal Mesoderm','Endothelium','Mixed mesoderm','Intermediate mesoderm','Somitic mesoderm','Cardiomyocytes'))))

# mesoderm @ moca
moca.meso <- subset(moca, cells=rownames(subset(moca@meta.data, celltype %in% c('Endothelial cells','Limb mesenchyme','Chondrocytes & osteoblasts','Myocytes','Chondroctye progenitors','Stromal cells','Notochord cells','Cardiac muscle lineages','Connective tissue progenitors','Early mesenchyme'))))

# mesoderm @ organoid
organoid.meso <- subset(organoid, cells=rownames(subset(organoid@meta.data, celltype %in% c(0,1,5,7,11,13,18))))

mouse.meso.list <- SplitObject(mouse.meso, split.by='stage')
moca.meso.list <- SplitObject(moca.meso, split.by='stage')

all.meso.list<-list()
all.meso.list <- mouse.meso.list
for (i in names(moca.meso.list)){
#  print(all.list[[i]])
	all.meso.list[[i]] <- moca.meso.list[[i]]
}
all.meso.list[['organoid']] <- organoid.meso

summary(all.meso.list)

load.list <- list()
emb.list <-list()
mouse.exp.list <- list()
organoid.compare.score.list <- list()

for (i in names(all.meso.list)){
	print (i)
	all.meso.list[[i]] <- ScaleData(all.meso.list[[i]], verbose=F)
	varf <- FindVariableFeatures(all.meso.list[[i]], selection.method = "vst", nfeatures = 2000, verbose=F)
	gene.use <- VariableFeatures(varf)
	all.meso.list[[i]] <- RunPCA(all.meso.list[[i]], features=gene.use, npcs=20, verbose=F)
	load.list[[i]] <- Loadings(all.meso.list[[i]], reduction='pca')
	emb.list[[i]] <- Embeddings(all.meso.list[[i]], reduction='pca')
	stage.factor <- all.meso.list[[i]]$stage
 	if (i %in% c('E9.5','E10.5','E11.5','E12.5','E13.5')){
 		stage.factor <- droplevels(stage.factor)
 	}

	mouse.exp.list[[i]] <- t(apply(t(emb.list[[i]]),1,function(x) tapply(x,stage.factor,mean)))
	
}

meso.array <- data.frame()

for (i in names(mouse.exp.list)){
	if (i %in% c('E6.5','E7.5','E8.5','E9.5','E10.5','E11.5','E12.5','E13.5','organoid')){
	A <- mouse.exp.list[[i]]
	for (k in names(mouse.exp.list)){
#	print(k)
	if (k %in% c('E6.5','E7.5','E8.5','E9.5','E10.5','E11.5','E12.5','E13.5','organoid')){
		B <- mouse.exp.list[[k]]
		meso.array[i,k] <- as.numeric(lsa::cosine(A[,],B[,]))
	}}}
}

write.csv(meso.array, file='./DATA/meso.total.CS.array.adjusted.csv')

################################
# Endoderm cosine similarity
# endoderm @ mouse
mouse.endo <- subset(mouse, cells=rownames(subset(mouse@meta.data, celltype %in% c('Def. endoderm','Gut'))))

# endoderm @ moca
moca.endo <- subset(moca, cells=rownames(subset(moca@meta.data, celltype %in% c('Hepatocytes','Epithelial cells'))))

# endoderm @ organoid
organoid.endo <- subset(organoid, cells=rownames(subset(organoid@meta.data, celltype %in% c(17))))

mouse.endo.list <- SplitObject(mouse.endo, split.by='stage')
moca.endo.list <- SplitObject(moca.endo, split.by='stage')

all.endo.list<-list()
all.endo.list <- mouse.endo.list
for (i in names(moca.endo.list)){
	all.endo.list[[i]] <- moca.endo.list[[i]]
}
all.endo.list[['organoid']] <- organoid.endo

summary(all.endo.list)

load.list <- list()
emb.list <-list()
mouse.exp.list <- list()
organoid.compare.score.list <- list()

for (i in names(all.endo.list)){
	print (i)
	all.endo.list[[i]] <- ScaleData(all.endo.list[[i]], verbose=F)
	varf <- FindVariableFeatures(all.endo.list[[i]], selection.method = "vst", nfeatures = 2000, verbose=F)
	
	gene.use <- VariableFeatures(varf)
	all.endo.list[[i]] <- RunPCA(all.endo.list[[i]], features=gene.use, npcs=20, verbose=F)
	load.list[[i]] <- Loadings(all.endo.list[[i]], reduction='pca')
	emb.list[[i]] <- Embeddings(all.endo.list[[i]], reduction='pca')
	stage.factor <- all.endo.list[[i]]$stage
 	if (i %in% c('E9.5','E10.5','E11.5','E12.5','E13.5')){
 		stage.factor <- droplevels(stage.factor)
 	}

	mouse.exp.list[[i]] <- t(apply(t(emb.list[[i]]),1,function(x) tapply(x,stage.factor,mean)))

}

endo.array <- data.frame()

for (i in names(mouse.exp.list)){
	if (i %in% c('E6.5','E7.5','E8.5','E9.5','E10.5','E11.5','E12.5','E13.5','organoid')){
	A <- mouse.exp.list[[i]]
	for (k in names(mouse.exp.list)){
#	print(k)
	if (k %in% c('E6.5','E7.5','E8.5','E9.5','E10.5','E11.5','E12.5','E13.5','organoid')){
		B <- mouse.exp.list[[k]]
		endo.array[i,k] <- as.numeric(lsa::cosine(A[,],B[,]))
	}}}
}

write.csv(endo.array, file='./DATA/endo.total.CS.array.adjusted.csv')

################################




















