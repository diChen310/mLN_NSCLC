####normal lymph node scRNA,2024-09-09####
library(Seurat)
ref.geo = read.delim(file='./variables/GSE131907_Lung_Cancer_cell_annotation.txt')

geo.data = readRDS(file='./variables/GSE131907_matrix.rds')
geo.data = CreateSeuratObject(counts = geo.data)
rownames(ref.geo)=ref.geo$Index


geo.data@meta.data$orig.ident = ref.geo[rownames(geo.data@meta.data),'Sample']
geo.data@meta.data = cbind(geo.data@meta.data,ref.geo[rownames(geo.data@meta.data),c(-1,-3)])



geo.data.ln = geo.data[,geo.data$Sample_Origin %in% c('nLN')]
saveRDS(geo.data.ln,file = './variables_v2/geo.data.ln_normal.Rds')



geo.data = readRDS(file='./variables/GSE131907_matrix.rds')
geo.data = CreateSeuratObject(counts = geo.data)
rownames(ref.geo)=ref.geo$Index


geo.data@meta.data$orig.ident = ref.geo[rownames(geo.data@meta.data),'Sample']
geo.data@meta.data = cbind(geo.data@meta.data,ref.geo[rownames(geo.data@meta.data),c(-1,-3)])



geo.data.ln = geo.data[,geo.data$Sample_Origin %in% c('nLN')]
saveRDS(geo.data.ln,file = './variables_v3/geo.data.ln_normal.Rds')

