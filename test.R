#!/usr/bin/Rscript
#annotate the combined sample of PAS and normal control by SingleR with ref data from HCA
library(Seurat)
print('loading ref data')
placenta_seu = readRDS('/home/liyaqi/DATA/placenta/seurat/placenta_seu.rds')
ct = read.csv('/home/liyaqi/DATA/placenta/FetalMaternal-Placenta-10x_cell_type_2020-03-12.csv')
rownames(ct)=ct$cell_name
rownames(placenta_seu@meta.data)=placenta_seu@meta.data$cell_names
placenta_seu = AddMetaData(object=placenta_seu,metadata=ct)
print(unique(placenta_seu@meta.data$annotated_cell_identity.ontology_label))

placenta_seu=placenta_seu[,placenta_seu@meta.data[,'annotated_cell_identity.ontology_label']!='NA']
placenta_seu=NormalizeData(placenta_seu)

print('loading combined data')
combine=readRDS('/home/liyaqi/DATA/placenta/seurat/Fcombined_remove.rds')
test=as.matrix(GetAssayData(combine, slot = "data",assay = "RNA"))
ref=as.matrix(GetAssayData(placenta_seu, slot = "data",assay ="RNA"))
labels=placenta_seu@meta.data$annotated_cell_identity.ontology_label
library(SingleR)

print('SingleR')
result=SingleR(test = test,ref = ref, labels = labels)
print(unique(result$labels))
saveRDS(placenta_seu,'/home/liyaqi/DATA/placenta/placenta_seu_normliz_ct_part.rds')
write.csv(result,'/home/liyaqi/DATA/placenta/result.csv')
combine@meta.data$SingleR_celltype = result$labels
saveRDS(combine,'/home/liyaqi/DATA/placenta/placenta_exp_combine_singleR_pre_ct.rds')
library(ggplot2)
dev.new()
p=Seurat::DimPlot(combined, group.by = "SingleR_celltype", label = T, repel = T)
plot(p)
ggsave('/home/liyaqi/DATA/placenta/SingleR_combine.png',width = 15,height = 10)
