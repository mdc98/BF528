#Individual Project 
library(dplyr)
library(patchwork)
library(Seurat)
library(tidyverse)
library(ggpubr)

cells <- readRDS("/projectnb/bf528/users/frazzled/project_4/programmer/pbmc_seurat2.rds")

cells_markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#p adjusted marker genes 
sig_marker_genes <- cells_markers[cells_markers$p_val_adj<0.05,]

#saving csv files 
marker_genes <- write_csv(cells_markers,"Marker_Genes")
sig_genes <- write_csv(sig_marker_genes,"PADJ_Marker_Genes")


#naming clusters
cluster_names <- c("Alpha","Beta","Delta","Gamma","Epsilon","Acinar","Ductal","Quiescent steliate","Activated steilate","Endothelial",
                   "Macrophage")
names(cluster_names) <- levels(cells)
cells <- RenameIdents(cells,cluster_names)


#get top 10 
top_5 <- cells_markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_logFC)

#heatmap 
map <- DoHeatmap(cells, features=top_5$gene, size=5, angle=75)





#UMAP 
UMAP_markers <-RunUMAP(cells, dims = 1:10) #label by cell type 
UMAP <- DimPlot(UMAP_markers, reduction= "umap", label = TRUE)
UMAP


