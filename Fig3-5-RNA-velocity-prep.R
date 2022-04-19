# Prepare RNA velocity embedding 
# In R
library(dplyr)
library(magrittr)
library(Seurat)


### Prepare main cluster files
organoid<-readRDS('./organoid.mt20.umap.rds')
organoid <- AddModuleScore(organoid, features='MKI67', name='KI67')
write.csv(Cells(organoid), file = "./velocyto/cellID_obs.csv", row.names = FALSE)
write.csv(Embeddings(organoid, reduction = "umap"), file = "./velocyto/cell_embeddings.csv")
write.csv(organoid@meta.data$seurat_clusters, file = "./velocyto/clusters.csv")


### Prepare ectoderm sub-cluster files
organoid<-readRDS('./org.Ecto.rds')
write.csv(Cells(organoid), file = "./velocyto/cellID_obs.Ecto.csv", row.names = FALSE)
write.csv(Embeddings(organoid, reduction = "umap"), file = "./velocyto/cell_embeddings.Ecto.csv")
write.csv(organoid@meta.data$seurat_clusters, file = "./velocyto/clusters.Ecto.csv")


### Prepare mesoderm sub-cluster files
organoid<-readRDS('./org.Meso.rds')
write.csv(Cells(organoid), file = "./velocyto/cellID_obs.Meso.csv", row.names = FALSE)
write.csv(Embeddings(organoid, reduction = "umap"), file = "./velocyto/cell_embeddings.Meso.csv")
write.csv(organoid@meta.data$seurat_clusters, file = "./velocyto/clusters.Meso.csv")


### Prepare endoderm sub-cluster files
organoid<-readRDS('./org.Endo.rds')
write.csv(Cells(organoid), file = "./velocyto/cellID_obs.Endo.csv", row.names = FALSE)
write.csv(Embeddings(organoid, reduction = "umap"), file = "./velocyto/cell_embeddings.Endo.csv")
write.csv(organoid@meta.data$seurat_clusters, file = "./velocyto/clusters.Endo.csv")
