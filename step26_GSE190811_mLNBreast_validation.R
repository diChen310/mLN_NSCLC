
####GSE190811####

library(Seurat)
library(SeuratData)
####Integrate all samples####
st.objs = c()
cell_ids = c()
st.input.folder = '/share/dichen/data/GSE190811/'
st.folders = list.dirs(path = st.input.folder,recursive = F)

for(st.f in st.folders){
  st.obj = Load10X_Spatial(st.f,filename = 'feature_bc_matrix.h5')
  names = strsplit(st.f,'//')[[1]][2]
  
  st.obj$region = names
  
  st.objs = append(st.objs,st.obj)
  print(paste('Read',st.f,'Finished!'))
  cell_ids = append(cell_ids,names)
}
st.obj = Read10X('./data/V1_nLN/filtered_feature_bc_matrix/')####Manully change
st.obj = CreateSeuratObject(st.obj, assay = 'Spatial')
st.obj$slice = 1
st.obj$region = 'nLN'####Manully change
img = Seurat::Read10X_Image(image.dir = './data/V1_nLN/spatial')####Manully change
Seurat::DefaultAssay(object = img) = 'Spatial'
img = img[colnames(x = st.obj)]
st.obj[['image']] = img
st.objs = append(st.objs,st.obj)
cell_ids = append(cell_ids,'nLN')
st.m = merge(st.objs[1][[1]],st.objs[c(-1,-5)],add.cell.ids = cell_ids[-5])
names(st.m@images) = cell_ids[-5]

st.m = PercentageFeatureSet(st.m, "^MT-", col.name = "percent_mito")
st.m = PercentageFeatureSet(st.m, "^HB[^(P)]", col.name = "percent_hb")


st.m$orig.ident = st.m$region
Idents(st.m)=st.m$orig.ident
P1=VlnPlot(st.m, features = "nCount_Spatial", pt.size = 0) + NoLegend()
P2=VlnPlot(st.m, features = "nFeature_Spatial", pt.size = 0) + NoLegend()
P3=VlnPlot(st.m, features = "percent_mito", pt.size = 0)+ NoLegend()
P4=VlnPlot(st.m, features = "percent_hb", pt.size = 0)+ NoLegend()

P1 / P2 / P3 / P4
meta.data = st.m@meta.data
meta.data$Source = ifelse(meta.data$orig.ident %in% c('A','B','C','D'),'mLN_Breast','nLN')
st.m$Source = 'mLN_Breast'
VlnPlot(st.m,features = c('CST1','IGHG4','CXCL13','FTH1'),group.by = 'Source',pt.size = 0,adjust = 19)



image.info = data.frame()
for(i in 1:length(unique(st.m$region))){
  
  image.i = st.m@images[[i]]
  image.i = image.i@coordinates
  
  image.info = rbind(image.info,image.i)
}

st.m = AddMetaData(st.m,image.info)
meta.data = st.m@meta.data

st.m = st.m[,st.m$nFeature_Spatial>100 & !is.na(st.m$tissue)]
st.m = SCTransform(st.m,variable.features.n = 2000,assay = "Spatial")
st.m = RunPCA(st.m, verbose = T, npcs = 30)
st.m = harmony::RunHarmony(st.m,group.by.vars = 'region',assay.use = 'SCT')
st.m = RunUMAP(st.m,reduction = "harmony", dims = 1:30)

colData <- st.m@meta.data

colData$spot <- rownames(colData)  
colData$in_tissue = colData$tissue
#colData <- colData[colData$in_tissue > 0, ]

counts <- Matrix::as.matrix(st.m@assays$Spatial@counts)

rowData <- data.frame(gene_id = rownames(counts),gene_name = rownames(counts))

sce <- SingleCellExperiment(assays = list(counts = counts), 
                            rowData = rowData, colData = colData)
metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "Visium"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

clusterPlot(sce, label="region", color = NA)  #make sure no overlap between samples
sce$row[sce$region == "B"] = 
  100 + sce$row[sce$region == "B"]
sce$row[sce$region == "C"] = 
  200 + sce$row[sce$region == "C"]
sce$row[sce$region == "D"] = 
  300 + sce$row[sce$region == "D"]

# sce$row[sce$region == "nLN"] = 
#   400 + sce$row[sce$region == "nLN"]

####3. preprocess####
sce = spatialPreprocess(sce, n.PCs = 50) #lognormalize, PCA


reducedDim(sce,'HARMONY') = st.m@reductions$harmony@cell.embeddings
reducedDim(sce,'UMAP.HARMONY') = st.m@reductions$umap@cell.embeddings

colnames(reducedDim(sce, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce$region))) +
  geom_point(size=0.1) +
  labs(color = "Sample") +
  theme_bw()
####4. clustering####


sce = spatialCluster(sce, use.dimred = "HARMONY", q = 20, 
                     platform="Visium", d=20,
                     init.method="mclust", model="t", gamma=3,
                     nrep=1000, burn.in=100) #use HARMONY
print(Sys.time())#16:00 to ?"2024-12-19 16:45:53 +08"
clusterPlot(sce, color = NA,
            size=0.1) + #plot clusters
  labs(title = "BayesSpace joint clustering")
#saveRDS(sce,'./variables_v2/sce_bayesSpace_validation.Rds')


st.m$bs.cluster = paste0('SC',sce[['spatial.cluster']])

Idents(st.m)=st.m$bs.cluster

de_markers = FindAllMarkers(st.m, assay = 'SCT',logfc.threshold = 0.1, 
                            test.use = "wilcox", 
                            min.pct = 0.01, 
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)
de_markers$cluster = paste0(de_markers$cluster,'m')
top_markers = filter(de_markers,p_val_adj<0.01) %>% group_by(cluster) %>% top_n(6,avg_log2FC) %>% filter(pct.1>0.3)

DoHeatmap(st.m[,sample(colnames(st.m),2500)],features = unique(top_markers$gene))
DotPlot(st.m,assay='SCT',features =  unique(top_markers$gene))+coord_flip()+
  theme_cowplot(font_size = 6)+
  scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))+
  scale_size_continuous(range = c(0,2.5))+
  theme(axis.text.x = element_text(angle = 90))

# 
# write.csv(de_markers,file = './tables_v2/top markers of spatial clusters_bayesspace SC21.csv')
# 


st.m = AddModuleScore(st.m,features = list(CXCL13.markers),name = 'CXCL13T_Score')

ggplot(st.m@meta.data,aes(x=region,fill=bs.cluster))+geom_bar(position = 'fill',width = 0.75)+
  theme_cowplot()+
  scale_fill_manual(values = cols.cl)+
  theme(axis.text.x = element_text(angle = 90))




####Merge nLN####

st.nln = st.obj
st.nln$bs.cluster = 'nLN'
st.nln = SCTransform(st.nln,variable.features.n = 2000,assay = 'Spatial')
st.m.2 = merge(st.m,st.nln)

st.nln = SCTransform(st.nln,variable.features.n = 2000,assay = 'Spatial')
st.m.2 = AddModuleScore(st.m.2,features = list(CXCL13.markers),name = 'CXCL13T_Score')
st.m = AddModuleScore(st.m,features = list(IFI30.top),name = 'IFI30_Score')

VlnPlot(st.m,features = c('CXCL13T_Score1'),
        group.by = 'bs.cluster',pt.size = 0,adjust = 5)
SpatialFeaturePlot(st.m.2,features = 'CXCL13T_Score1',pt.size.factor = 1.25,images = c('A','B','C','D'))


SpatialFeaturePlot(st.m.2,features = 'IFI30_Score1',pt.size.factor = 1.25,images = c('A','B','C','D'))

st.m$LS.like = ifelse(st.m$bs.cluster %in% c('SC6','SC20'),'LS-like','Others')
Idents(st.m) = st.m$LS.like
markers_ls = FindMarkers(st.m,ident.1 = 'LS-like',logfc.threshold = 0.3,only.pos = T)


SpatialDimPlot(st.m,group.by = 'LS.like',pt.size.factor = 1.25,images = c('A','B','C','D'))


VlnPlot(st.m,features = c('CXCL13T_Score1'),group.by = 'LS.like',pt.size = 0)+ggpubr::stat_compare_means()
SpatialFeaturePlot(st.m.2,features = c('CCL21','CXCL13','MS4A1','CXCL13T_Score1'),pt.size.factor = 1.25,
                   images = c('A','B','C','D'),
                   ncol = 4)
SpatialFeaturePlot(st.m.2,features = c('CCL21','CXCL13','MS4A1','CXCL13T_Score1'),pt.size.factor = 1.25,
                   images = c('image'),
                   ncol = 1)
write.csv(st.m@meta.data,file = './tables_v3/Fig S_mLN_Breast_part1.csv')
meta.data = read.csv(file = './tables_v3/Fig S_mLN_Breast_part1.csv',row.names = 1)
st.m = AddMetaData(st.m,meta.data)



SpatialFeaturePlot(st.m.2,features = 'IFI30',pt.size.factor = 1.25,images = c('A','B','C','D','image'))
SpatialFeaturePlot(st.m.2,features = 'FTH1',pt.size.factor = 1.25,images = c('A','B','C','D','image'))
SpatialFeaturePlot(st.m.2,features = 'SPP1',pt.size.factor = 1.25,images = c('A','B','C','D','image'))

st.m.2$Source = ifelse(st.m.2$orig.ident %in% c('A','B','C','D'),'mLN_Breast','nLN')
st.m.2$LS.like = ifelse(st.m.2$bs.cluster %in% c('SC6','SC20','nLN'),'LS-like','Others')
VlnPlot(st.m.2[,st.m.2$LS.like == 'LS-like'],features = c('CXCL13T_Score1'),group.by = 'Source',pt.size = 0)+ggpubr::stat_compare_means()

st.m.2 = AddModuleScore(st.m.2,features = list(c(IFI30.top,'SPP1')),name = 'IFI30MP_Score')


VlnPlot(st.m.2,features = c('IFI30','SPP1','IFI30MP_Score1'),group.by = 'Source',pt.size = 0,adjust = 3,ncol = 4)


st.m.2$IFI30 = st.m.2@assays$SCT@data['IFI30',]
st.m.2$SPP1 = st.m.2@assays$SCT@data['SPP1',]
st.m.2$FTH1 = st.m.2@assays$SCT@data['FTH1',]
write.csv(st.m.2@meta.data[,c('orig.ident','Source','IFI30','SPP1','FTH1','CXCL13T_Score1','IFI30MP_Score1')],file = './tables_v3/Fig S_mLN_Breast_part2.csv')


# saveRDS(st.m.2,file = './variables_v2/st.m_mLNBreast&nLN.Rds')
# 
# st.m.2 = readRDS(file = './variables_v2/st.m_mLNBreast&nLN.Rds')
