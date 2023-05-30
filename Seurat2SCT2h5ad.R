library(SeuratDisk)
library(Seurat)
sc=readRDS('../scRNA_anno_roughly.rds')
seurat <- SCTransform(sc, verbose = FALSE)
SaveH5Seurat(seurat, filename = "snrna_fetal_brain.h5Seurat")
Convert("snrna_fetal_brain.h5Seurat", dest = "h5ad")
colnames(sc@meta.data)

##
Adding data for integrated
Adding data for integrated
Adding scale.data for integrated
Adding variable features for integrated
No feature-level metadata found for integrated
Adding counts for SCT
Adding data for SCT
Adding scale.data for SCT
Adding variable features for SCT
No feature-level metadata found for SCT
Writing out SCTModel.list for SCT
Adding cell embeddings for pca
Adding loadings for pca
No projected loadings for pca
Adding standard deviations for pca
No JackStraw data for pca
Adding cell embeddings for umap
No loadings for umap
No projected loadings for umap
No standard deviations for umap
No JackStraw data for umap
Adding cell embeddings for tsne
No loadings for tsne
No projected loadings for tsne
No standard deviations for tsne
No JackStraw data for tsne
> Convert("snrna_fetal_brain.h5Seurat", dest = "h5ad")
Validating h5Seurat file
Adding scale.data from SCT as X
Adding data from SCT as raw
Transfering meta.data to obs
Adding dimensional reduction information for tsne (global)
Adding dimensional reduction information for umap (global)
> colnames(sc@meta.data)
[1] "orig.ident"           "nCount_RNA"           "nFeature_RNA"
[4] "organ"                "batch"                "percent.MT"
[7] "integrated_snn_res.1" "seurat_clusters"
##
