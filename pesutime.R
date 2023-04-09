library(monocle3)
library(Seurat)
library(SummarizedExperiment)
setwd('~/Dropbox/Single_cell/figs/local_figure/')
immune.combined = readRDS('placenta_exp_combine_singleR_pre_ct.rds')
DefaultAssay(immune.combined) <- "RNA"
PAS_removed = immune.combined[,immune.combined@meta.data$orig.ident=='PAS']
Norm_removed = immune.combined[,immune.combined@meta.data$orig.ident=='Norm']

seurat_obj = Norm_removed
seurat_obj = NormalizeData(seurat_obj)

expression_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "data",assay = "RNA"))
cell_metadata <- seurat_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(seurat_obj@assays$RNA))
rownames(gene_annotation) <- gene_annotation$gene_short_name


cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch")
cds <- reduce_dimension(cds,preprocess_method ='PCA',reduction_method ='UMAP')
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "SingleR_celltype")
ggsave('./pesudotime_Norm_ct.png',width = 4,height = 4)

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "SingleR_celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
library(ggplot2)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
ggsave('./pesudotime_Norm_pesudo.png',width = 4,height = 4)
