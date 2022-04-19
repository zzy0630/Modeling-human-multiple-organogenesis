# Calculate RNA velocity of main clusters
# Fig. 3B
# In python3 (v.3.9.6)
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

sample= anndata.read_loom("./velocyto/MidOrganoid.loom") # velocyto generated loom
sample_obs = pd.read_csv("./velocyto/cellID_obs.csv") 
cell_clusters = pd.read_csv("./velocyto/clusters.csv") 
umap = pd.read_csv("./velocyto/cell_embeddings.csv") 

# Modify cell ID in different matrix
# 
sample.obs=sample.obs.rename(index = lambda x: x.replace('MidOrganoid:', ''))
sample.obs=sample.obs.rename(index = lambda x: x.replace('x', ''))
#
sample_obs.x=sample_obs.x.replace({"-1":""},regex=True)
# 
umap[["Unnamed: 0"]]=umap[["Unnamed: 0"]].replace({"-1":""},regex=True)
umap = umap.rename(columns = {"Unnamed: 0":'CellID'})
#
cell_clusters[["Unnamed: 0"]] = umap[["CellID"]]
#

#######
sample = sample[np.isin(sample.obs.index,sample_obs["x"])]
sample_index = pd.DataFrame(sample.obs.index)
sample_index = sample_index.rename(columns = {0:'CellID'})
umap_ordered = sample_index.merge(umap,on="CellID")
cell_clusters[["Unnamed: 0"]]=cell_clusters[["Unnamed: 0"]].replace({"_1":""},regex=True)

cell_clusters = cell_clusters.rename(columns = {"Unnamed: 0":'CellID'})
cell_clusters = sample_index.merge(cell_clusters,on="CellID")
cell_clusters = sample_index.merge(cell_clusters,on="CellID")

umap_ordered = umap_ordered.iloc[:,1:]
sample.obsm['X_umap'] = umap_ordered.values
# cell_clusters_ordered=cell_clusters.iloc[:,2]
cell_clusters_ordered=cell_clusters.iloc[:,1]
sample.obs['cell_clusters']=cell_clusters_ordered.values

######
# Preprocess the data 
scv.pp.filter_and_normalize(sample,min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(sample, n_pcs=30, n_neighbors=30)
# Run RNA Velocity
scv.tl.velocity(sample, mode = "stochastic")
scv.tl.velocity_graph(sample)
# scv.pl.velocity_embedding(sample, basis='X_umap',arrow_size=5)
# Modify colors
ident_colours = ["#F8766D","#EA8331" ,"#D89000" ,"#C09B00" ,"#A3A500", "#7CAE00" ,"#39B600",
"#00BB4E", "#00BF7D" ,"#00C1A3", "#00BFC4", "#00BAE0", "#00B0F6", "#35A2FF",
"#9590FF", "#C77CFF", "#E76BF3", "#FA62DB", "#FF62BC", "#FF6A98"]
col =['#d5dfff','#cc0000']
# Prject the velocities on UMAP cluster # Fig. 3B
scv.pl.velocity_embedding_stream(sample, basis='X_umap',color = "cell_clusters",palette = ident_colours)
# Project the velocities on UMAP KI67 # Fig. 3B
scv.pl.velocity_embedding_stream(sample, basis='X_umap',color = "MKI67", color_map=col)

###### 
# Pseudotime # Fig. 3B
# scv.pl.velocity_graph(sample, threshold=.1)
scv.tl.velocity_pseudotime(sample,n_dcs=1)
# computing terminal states
# WARNING: Uncertain or fuzzy root cell identification. Please verify.
#    identified 2 regions of root cells and 4 regions of end points .
#    finished (0:00:02) --> added
#    'root_cells', root cells of Markov diffusion process (adata.obs)
#    'end_points', end points of Markov diffusion process (adata.obs)
scv.pl.scatter(sample, color='velocity_pseudotime', cmap='gnuplot')


#######################################################
# RNA velocity of sub-clusters
ori= anndata.read_loom("./velocyto/MidOrganoid.loom")

####################
# Endo # fig. S9B
sample = ori
sample_obs = pd.read_csv("./velocyto/cellID_obs.Endo.csv") 
cell_clusters = pd.read_csv("./velocyto/clusters.Endo.csv") 
umap = pd.read_csv("./velocyto/cell_embeddings.Endo.csv") 

# Modify cell ID in different matrix
# 
sample.obs=sample.obs.rename(index = lambda x: x.replace('MidOrganoid:', ''))
sample.obs=sample.obs.rename(index = lambda x: x.replace('x', ''))
#
sample_obs.x=sample_obs.x.replace({"-1":""},regex=True)
# 
umap[["Unnamed: 0"]]=umap[["Unnamed: 0"]].replace({"-1":""},regex=True)
umap = umap.rename(columns = {"Unnamed: 0":'CellID'})
#
cell_clusters[["Unnamed: 0"]] = umap[["CellID"]]

#######
sample = sample[np.isin(sample.obs.index,sample_obs["x"])]
sample_index = pd.DataFrame(sample.obs.index)
sample_index = sample_index.rename(columns = {0:'CellID'})
umap_ordered = sample_index.merge(umap,on="CellID")
cell_clusters[["Unnamed: 0"]]=cell_clusters[["Unnamed: 0"]].replace({"_1":""},regex=True)

cell_clusters = cell_clusters.rename(columns = {"Unnamed: 0":'CellID'})
cell_clusters = sample_index.merge(cell_clusters,on="CellID")
cell_clusters = sample_index.merge(cell_clusters,on="CellID")

umap_ordered = umap_ordered.iloc[:,1:]
sample.obsm['X_umap'] = umap_ordered.values
# cell_clusters_ordered=cell_clusters.iloc[:,2]
cell_clusters_ordered=cell_clusters.iloc[:,1]
sample.obs['cell_clusters']=cell_clusters_ordered.values

######
# Run Velocity
scv.pp.filter_and_normalize(sample,min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(sample, n_pcs=30, n_neighbors=30)
scv.tl.velocity(sample, mode = "stochastic")
scv.tl.velocity_graph(sample)
# scv.pl.velocity_embedding(sample, basis='X_umap',arrow_size=5)
ident_colours = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

scv.pl.velocity_embedding_stream(sample, basis='X_umap',color = "cell_clusters",palette = ident_colours, legend_fontsize=0)

#######################
# Ectoderm # fig. S8B
sample = ori
sample_obs = pd.read_csv("./velocyto/cellID_obs.Ecto.csv") 
cell_clusters = pd.read_csv("./velocyto/clusters.Ecto.csv") 
umap = pd.read_csv("./velocyto/cell_embeddings.Ecto.csv") 

# Modify cell ID in different matrix
# 
sample.obs=sample.obs.rename(index = lambda x: x.replace('MidOrganoid:', ''))
sample.obs=sample.obs.rename(index = lambda x: x.replace('x', ''))
#
sample_obs.x=sample_obs.x.replace({"-1":""},regex=True)
# 
umap[["Unnamed: 0"]]=umap[["Unnamed: 0"]].replace({"-1":""},regex=True)
umap = umap.rename(columns = {"Unnamed: 0":'CellID'})
#
cell_clusters[["Unnamed: 0"]] = umap[["CellID"]]

#######
sample = sample[np.isin(sample.obs.index,sample_obs["x"])]
sample_index = pd.DataFrame(sample.obs.index)
sample_index = sample_index.rename(columns = {0:'CellID'})
umap_ordered = sample_index.merge(umap,on="CellID")
cell_clusters[["Unnamed: 0"]]=cell_clusters[["Unnamed: 0"]].replace({"_1":""},regex=True)

cell_clusters = cell_clusters.rename(columns = {"Unnamed: 0":'CellID'})
cell_clusters = sample_index.merge(cell_clusters,on="CellID")
cell_clusters = sample_index.merge(cell_clusters,on="CellID")

umap_ordered = umap_ordered.iloc[:,1:]
sample.obsm['X_umap'] = umap_ordered.values
# cell_clusters_ordered=cell_clusters.iloc[:,2]
cell_clusters_ordered=cell_clusters.iloc[:,1]
sample.obs['cell_clusters']=cell_clusters_ordered.values

######
# Run Velocity
scv.pp.filter_and_normalize(sample,min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(sample, n_pcs=30, n_neighbors=30)
scv.tl.velocity(sample, mode = "stochastic")
scv.tl.velocity_graph(sample)
# scv.pl.velocity_embedding(sample, basis='X_umap',arrow_size=5)
ident_colours = ["#F8766D", "#E58700", "#C99800", "#A3A500", "#6BB100", "#00BA38", "#00BF7D",
"#00C0AF", "#00BCD8", "#00B0F6" ,"#619CFF", "#B983FF" ,"#E76BF3", "#FD61D1",
"#FF67A4"]
scv.pl.velocity_embedding_stream(sample, basis='X_umap',color = "cell_clusters",palette = ident_colours, legend_fontsize=0)

# testing top likelihood genes # fig. S7B
# 
scv.tl.recover_dynamics(sample)
top_genes = sample.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(sample, var_names=top_genes, groupby='cell_clusters')

kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
scv.pl.scatter(sample, basis=top_genes[:15], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[15:30], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[30:45], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[45:50], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[45:60], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[60:75], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[75:90], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[90:105], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
# Recompute velocity
scv.tl.velocity(sample, diff_kinetics=True)
scv.tl.velocity_graph(sample)
scv.pl.velocity_embedding(sample, dpi=120, arrow_size=2, arrow_length=2)
scv.pl.velocity_embedding(sample, dpi=120, arrow_size=20, arrow_length=5)

####################
# Meso # fig. S9B
sample = ori
sample_obs = pd.read_csv("./velocyto/cellID_obs.Meso.csv") 
cell_clusters = pd.read_csv("./velocyto/clusters.Meso.csv") 
umap = pd.read_csv("./velocyto/cell_embeddings.Meso.csv") 

# Modify cell ID in different matrix
# 
sample.obs=sample.obs.rename(index = lambda x: x.replace('MidOrganoid:', ''))
sample.obs=sample.obs.rename(index = lambda x: x.replace('x', ''))
#
sample_obs.x=sample_obs.x.replace({"-1":""},regex=True)
# 
umap[["Unnamed: 0"]]=umap[["Unnamed: 0"]].replace({"-1":""},regex=True)
umap = umap.rename(columns = {"Unnamed: 0":'CellID'})
#
cell_clusters[["Unnamed: 0"]] = umap[["CellID"]]

#######
sample = sample[np.isin(sample.obs.index,sample_obs["x"])]
sample_index = pd.DataFrame(sample.obs.index)
sample_index = sample_index.rename(columns = {0:'CellID'})
umap_ordered = sample_index.merge(umap,on="CellID")
cell_clusters[["Unnamed: 0"]]=cell_clusters[["Unnamed: 0"]].replace({"_1":""},regex=True)

cell_clusters = cell_clusters.rename(columns = {"Unnamed: 0":'CellID'})
cell_clusters = sample_index.merge(cell_clusters,on="CellID")
cell_clusters = sample_index.merge(cell_clusters,on="CellID")

umap_ordered = umap_ordered.iloc[:,1:]
sample.obsm['X_umap'] = umap_ordered.values
# cell_clusters_ordered=cell_clusters.iloc[:,2]
cell_clusters_ordered=cell_clusters.iloc[:,1]
sample.obs['cell_clusters']=cell_clusters_ordered.values

######
# Run Velocity
scv.pp.filter_and_normalize(sample,min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(sample, n_pcs=30, n_neighbors=30)
scv.tl.velocity(sample, mode = "stochastic")
scv.tl.velocity_graph(sample)
# scv.pl.velocity_embedding(sample, basis='X_umap',arrow_size=5)
ident_colours = ["#F8766D", "#E58700", "#C99800", "#A3A500", "#6BB100", "#00BA38", "#00BF7D",
"#00C0AF", "#00BCD8", "#00B0F6" ,"#619CFF", "#B983FF" ,"#E76BF3", "#FD61D1",
"#FF67A4"]
scv.pl.velocity_embedding_stream(sample, basis='X_umap',color = "cell_clusters",palette = ident_colours, legend_fontsize=0)

# testing top likelihood genes # fig. S7C
# 
scv.tl.recover_dynamics(sample)
top_genes = sample.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(sample, var_names=top_genes, groupby='cell_clusters')

kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
scv.pl.scatter(sample, basis=top_genes[:15], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[15:30], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[30:45], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[45:50], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[45:60], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[60:75], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[75:90], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
scv.pl.scatter(sample, basis=top_genes[90:105], ncols=5, add_outline='fit_diff_kinetics', **kwargs)
# Recompute velocity
scv.tl.velocity(sample, diff_kinetics=True)
scv.tl.velocity_graph(sample)
scv.pl.velocity_embedding(sample, dpi=120, arrow_size=2, arrow_length=2)
scv.pl.velocity_embedding(sample, dpi=120, arrow_size=20, arrow_length=5)
















