.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')

prepare_coord = function(st.f.input,names=names(st.f.input@images)){
  
  N=length(st.f.input@images)
  
  st.f.coord = data.frame()
  
  for(i in 1:N){
    print(i)
    i.coord = st.f.input@images[[i]]@coordinates
    
    i.coord$sample_name= names[i]
    
    st.f.coord=rbind(st.f.coord,i.coord)
  }
  return(st.f.coord)
}
library(CARD)
library(monocle3)
library(Seurat)
library(jcolors)
library(RColorBrewer)
library(reshape2)
library(dplyr)

setwd('/share/dichen/lungNodeM')

####1. predict cell by CARD####

## set the colors. 
cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559","#32709f",
                 '#b8a3c4','#60428d','#9f5633','grey','#c92e2a','#c79795','#d9c17c',"#d0896b",'#5d94b4',"black","blue",'cyan')

cols.celltype =cols.default[1:9]
names(cols.celltype)=c('T','B','Plasma','Epi','Mast','Myeloid','Fibro','Endo','NK')

cds = readRDS(file = './variables_v2/inter_cds_refine.Rds')


#cds.mLNT = cds[,cds[['PatientID']] == 'P6']
cds.mLNT = cds[,cds[['Source']] == 'PLM']

marker.pi = top_markers(cds.mLNT,genes_to_test_per_group = 100,cores=12,group_cells_by = 'subCellType')

marker.pi = filter(marker.pi,marker_test_p_value<0.001)
saveRDS(marker.pi,file = './variables_v2/cds_PLM_subCellType_markers.Rds')
#cds.pi = cds[,cds[['From']] == 'nLN'] for NLN
#saveRDS(marker.pi,file = './variables_v2/cds_mLNT_LUSC_markers.Rds')
sc.data = assay(cds.mLNT)
# meta.data = as.data.frame(colData(cds.mLNT))
# meta.data$majorCellType = as.character(meta.data$majorCellType)
# table(meta.data$majorCellType)
# st.obj = Read10X('./result-1/P7-T-WK_outs//filtered_feature_bc_matrix/')####Manully change
# 
# st.obj = Read10X('./ST-2/1.Count/P5T/filtered_feature_bc_matrix/')
# st.obj = CreateSeuratObject(st.obj, assay = 'Spatial')
# st.obj$slice = 1
# st.obj$region = 'P5T'####Manully change
# img = Seurat::Read10X_Image(image.dir = paste0('./ST-2/1.Count/P5T/','/spatial'))####Manully change
# 
# img = Seurat::Read10X_Image(image.dir = paste0('./result-1/P7-T-WK_outs/','/spatial'))####Manully change
# 
# Seurat::DefaultAssay(object = img) = 'Spatial'
# img = img[colnames(x = st.obj)]
# st.obj[['image']] = img
# 
# SpatialFeaturePlot(st.obj,features = c('CXCL13','CST1'))
# 
# SpatialFeaturePlot(st.obj,features = c('CCL5','GZMA','PDCD1','CTLA4'))
# 
# table(cds.mLNT[['majorCellType']])


meta.all = as.data.frame(colData(cds.mLNT))
meta.all$subCellTypeV3 = meta.all$subCellType
sel.subtypes = c('ALDOA+ T',"ARHGAP15+ T",'CD4 Naive T','CD4 Resident Effector Memory',
                 'CD4 Effector Memory', "CXCL13+ T", "Exhausted T",'Treg',
                 "CD8 Cytotoxic T","CD8 Exhausted T",'NK',"¦Ã¦Ä T",
                 'B','Naive B','GC B','Activated B',
                 'Plasma','IgG low plasma','ZBTB20+ B')
meta.all$subCellTypeV3[meta.all$subCellTypeV3 %in% sel.subtypes == F] = 'Others'
table(meta.all$subCellTypeV3)
cds.mLNT[['subCellTypeV3']] = meta.all$subCellTypeV3
st.obj = readRDS(file='./variables_v2/st.m_bayesSpace.Rds')


####

st.obj = NormalizeData(st.obj)
#st.obj = FindVariableFeatures(st.obj,nfeatures = 3000)

st.obj.count = st.obj@assays$Spatial@counts[intersect(rownames(st.obj),unique(marker.pi$gene_id)),]
dim(st.obj.count) #2676 48535,2024-10-31

image.info = data.frame()
for(i in 1:length(unique(st.m$region))){
  
  image.i = st.m@images[[i]]
  image.i = image.i@coordinates
  
  image.info = rbind(image.info,image.i)
}

st.obj.coord = image.info
colnames(st.obj.coord)[c(2,3)]=c('x','y')

CARD_obj = createCARDObject(
  sc_count = as.matrix(sc.data),
  sc_meta = meta.all,
  spatial_count = st.obj.count,
  spatial_location = st.obj.coord,
  ct.varname = "subCellTypeV3",
  ct.select = unique(meta.all$subCellTypeV3),
  sample.varname = "orig.ident",
  minCountGene = 50,
  minCountSpot = 3) 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
print(CARD_obj@Proportion_CARD[1:2,])
saveRDS(CARD_obj,file = './variables_v2/CARD_obj.Rds')
CARD_obj = readRDS(file = './variables_v2/CARD_obj.Rds')
# CARD_obj2 = createCARDObject(
#   sc_count = as.matrix(sc.data),
#   sc_meta = meta.data,
#   spatial_count = st.obj.count,
#   spatial_location = st.obj.coord,
#   ct.varname = "subCellType",
#   ct.select = unique(meta.data$majorCellType),
#   sample.varname = "orig.ident",
#   minCountGene = 100,
#   minCountSpot = 5) 
# CARD_obj2 = CARD_deconvolution(CARD_object = CARD_obj2)
# print(CARD_obj2@Proportion_CARD[1:2,])


# p1 <- CARD.visualize.pie(
#   proportion = CARD_obj@Proportion_CARD,
#   spatial_location = CARD_obj@spatial_location,
#   colors = cols.celltype) ### You can choose radius = NULL or your own radius number
# print(p1)
# 
# 
# ## visualize the spatial distribution of the cell type proportion
# p2 <- CARD.visualize.prop(
#   proportion = CARD_obj@Proportion_CARD,        
#   spatial_location = CARD_obj@spatial_location, 
#   ct.visualize = c('T','NK','B','Plasma','Epi','Mast','Myeloid','Fibro','Endo'),                 ### selected cell types to visualize
#   colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
#   NumCols = 4)                             ### point size in ggplot2 scatterplot  
# print(p2)

prediction.scores = as.data.frame(CARD_obj@Proportion_CARD)
# prediction.scores$celltype_prediction <- NA
# dim(prediction.scores)
# for(i in 1:nrow(prediction.scores)){
#   prediction.scores$celltype_prediction[i] <- colnames(prediction.scores)[prediction.scores[i,1:(ncol(prediction.scores)-1)] == max(as.double(prediction.scores[i,1:(ncol(prediction.scores)-1)]))]
# }
st.obj = AddMetaData(st.obj,metadata = prediction.scores)

SpatialFeaturePlot(st.obj,images = c('image'), features = c('T','NK','B','Plasma','Epi','Mast','Myeloid','Fibro','Endo'),
                   pt.size.factor = 1.5,stroke = 0,
                   ncol = 3, alpha = c(0.1,1),crop = TRUE,image.alpha = 0)

SpatialFeaturePlot(st.obj,images = c('image'), features = c('LYZ','EPCAM','CD3D','CD3E','CD79A','CD79B','IGHG4','MZB1','COL1A1'),
                   pt.size.factor = 2,stroke = 0,
                   ncol = 3, alpha = c(0.6),crop = TRUE,image.alpha = 0)
SpatialFeaturePlot(st.obj,images = c('image'), features = c('CCL19','CCL21','CD3D','CD3E','CD79A','CD79B','IGHG4','MZB1','COL1A1'),
                   pt.size.factor = 2,stroke = 0,
                   ncol = 3, alpha = c(0.6),crop = TRUE,image.alpha = 0)

SpatialDimPlot(st.obj,images = c('image'),cols = cols.celltype,group.by = 'celltype_prediction',pt.size.factor = 1.2,alpha = 1,image.alpha = 0)


Idents(st.obj) = st.obj$celltype_prediction
st.diff = FindMarkers(st.obj,ident.1 = 'Endo')
saveRDS(prediction.scores,file = './variables_v2//SubCellTypeAll_predScores.Rds')
prediction.scores = readRDS(file = './variables_v2//SubCellTypeAll_predScores.Rds')
VlnPlot(st.obj,features = c('CXCL13..T'),group.by = 'cluster.name',pt.size = 0)

meta.data = st.obj@meta.data
meta.data$cluster.name = factor(meta.data$cluster.name,levels = c('LS-like','Epithelial','Fibroblast','Plasma', 'MT-high',
                                                                  'Plasma+Myeloid','Myeloid'))
ggplot(meta.data,aes(x=cluster.name,fill=cluster.name,y=CXCL13..T))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+
  ggpubr::stat_compare_means()+
  scale_fill_manual(values = cols.cl)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

ggplot(meta.data,aes(x=cluster.name,fill=cluster.name,y=Plasma))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+
  ggpubr::stat_compare_means()+
  scale_fill_manual(values = cols.cl)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

ggplot(meta.data,aes(x=cluster.name,fill=cluster.name,y=GC.B))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+
  facet_wrap(~From)+
  scale_fill_manual(values = cols.cl)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

meta.data.LS = meta.data[meta.data$cluster.name == 'LS-like',]
ggplot(meta.data.LS,aes(x=From,fill=From,y=GC.B))+geom_boxplot(outlier.size = 0.2)+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(2,3),c(1,3)))+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())

meta.data.LS.2 = melt(meta.data.LS,measure.vars = c("CXCL13..T","CD4.Resident.Effector.Memory","CD4.Effector.Memory"  ,
                                                    "CD4.Naive.T","Treg", "CD8.Cytotoxic.T","ALDOA..T",
                                                    "CD8.Exhausted.T","Exhausted.T","ARHGAP15..T"  ,                "¦Ã¦Ä.T"  ))
ggplot(meta.data.LS.2,aes(x=variable,fill=variable,y=value))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+
  facet_wrap(~From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

meta.data.LS.3 = melt(meta.data.LS,measure.vars = c( "Plasma",   "Activated.B" , "ZBTB20..B",
                                                     "GC.B", "Naive.B" ))
ggplot(meta.data.LS.3,aes(x=variable,fill=variable,y=value))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())
SpatialFeaturePlot(st.obj,images = c('P6L','P6T'),
                   features = c('CXCL13..T','GC.B'),
                   pt.size.factor = 1.3,
                   stroke = 0,
                   alpha = 0.8,image.alpha = 0.6,
                   ncol = 2,
                   combine = T)



####For B/Plasma -- Fibro co-localization####
merge.res$pair = paste(merge.res$Var1,merge.res$Var2)
merge.res.inv = merge.res[merge.res$pair %in% c('Plasma Fibro','Fibro Plasma','B Fibro','Fibro B'),]
merge.res.inv$pair = dplyr::recode(merge.res.inv$pair,
                                   'Fibro B'='B Fibro',
                                   'Fibro Plasma' = 'Plasma Fibro')
ggplot(merge.res.inv,aes(x=pair,y=log2(value+1),fill=pair))+geom_boxplot()+
  stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank())
####3. Merge object####
st.objs = c()
cell_ids = c()
st.input.folder = './result-1/'
st.folders = list.dirs(path = st.input.folder,recursive = F)
st.folders = st.folders[grepl('outs',st.folders)]
st.folders = st.folders[c(1,2,7,8)]

for(st.f in st.folders){
  st.obj = Read10X(data.dir = paste0(st.f,'/filtered_feature_bc_matrix'))
  names = strsplit(st.f,'//')[[1]][2]
  names = strsplit(names,'-')[[1]][-3]
  st.obj = CreateSeuratObject(st.obj, project = paste(names,collapse = '') ,assay = 'Spatial')
  st.obj$slice = 1
  st.obj$region = substr(names[2],1,1)
  img = Seurat::Read10X_Image(image.dir = paste0(st.f,'/spatial'))
  Seurat::DefaultAssay(object = img) = 'Spatial'
  img = img[colnames(x = st.obj)]
  st.obj[['image']] = img
  st.objs = append(st.objs,st.obj)
  print(paste('Read',st.f,'Finished!'))
  cell_ids = append(cell_ids,paste(names,collapse = ''))
}

st.folders.2 = list.dirs(path = './ST-2/1.Count/',recursive = F)

for(st.f in st.folders.2){
  st.obj = Read10X(data.dir = paste0(st.f,'/filtered_feature_bc_matrix'))
  names = strsplit(st.f,'//')[[1]][2]
  st.obj = CreateSeuratObject(st.obj, project = names ,assay = 'Spatial')
  st.obj$slice = 1
  st.obj$region = substr(names,3,3)
  img = Seurat::Read10X_Image(image.dir = paste0(st.f,'/spatial'))
  Seurat::DefaultAssay(object = img) = 'Spatial'
  img = img[colnames(x = st.obj)]
  st.obj[['image']] = img
  st.objs = append(st.objs,st.obj)
  print(paste('Read',st.f,'Finished!'))
  cell_ids = append(cell_ids,names)
}


st.m = merge(st.objs[1][[1]],st.objs[-1],add.cell.ids = cell_ids)
names(st.m@images) = cell_ids

st.m = PercentageFeatureSet(st.m, "^MT-", col.name = "percent_mito")
st.m = PercentageFeatureSet(st.m, "^HB[^(P)]", col.name = "percent_hb")

P1=VlnPlot(st.m, features = "nCount_Spatial", pt.size = 0) + NoLegend()
P2=VlnPlot(st.m, features = "nFeature_Spatial", pt.size = 0) + NoLegend()
P3=VlnPlot(st.m, features = "percent_mito", pt.size = 0)+ NoLegend()
P4=VlnPlot(st.m, features = "percent_hb", pt.size = 0)+ NoLegend()

P1 / P2 / P3 / P4

st.f =  subset(st.m, subset = nFeature_Spatial > 250  & nCount_Spatial > 200  & percent_mito < 25 )




pred.all = data.frame()

tags = c('p3t'='P3T','p3l'='P3L12','p5t'='P5T','p5l'='P5L','p6t'='P6T','p6l'='P6L','p7t'='P7T','p7l'='P7L5')
for(pi in c('p3','p5','p6','p7')){
  for(from in c('t','l')){
    pred.i = readRDS(file = paste0('/share/dichen/lungNodeM/',pi,from,'_predScores.Rds'))
    tag = as.character( tags[paste0(pi,from)])
    rownames(pred.i) = paste0(tag,'_',rownames(pred.i))
    pred.all = rbind(pred.all,pred.i)
  }
}

st.f = AddMetaData(st.f,metadata=pred.all)

saveRDS(st.f,file = './variables_v2/stIntegration-withCARDPred.Rds')


cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559","#32709f",
                 '#b8a3c4','#60428d','#9f5633','grey','#c92e2a','#c79795','#d9c17c',"#d0896b",'#5d94b4',"black","blue",'cyan')

cols.celltype =cols.default[1:9]
names(cols.celltype)=c('T','B','Plasma','Epi','Mast','Myeloid','Fibro','Endo','NK')

SpatialDimPlot(st.f,images = c('P3L12','P5L','P6L','P7L5','P3T','P5T','P6T','P7T'),
               cols = cols.celltype,
               group.by = 'celltype_prediction',
               pt.size.factor = 1.3,
               stroke = 0,
               alpha = 1,image.alpha = 0,
               label.box = F,
               ncol = 4,
               combine = T)

SpatialFeaturePlot(st.f,images = c('P3L12','P3T'),
               features = c('EPCAM','CD3D','ACTA2','MS4A1','MZB1','LYZ','PECAM1','TPSAB1','KLRB1'),
               pt.size.factor = 1.3,
               stroke = 0,
               alpha = 1,image.alpha = 0,
               ncol = 6,
               combine = T)

SpatialFeaturePlot(st.f,images = c('P5L7','P5T'),
                   features = c('EPCAM','CD3D','ACTA2','MS4A1','MZB1','LYZ','PECAM1','TPSAB1','KLRB1'),
                   pt.size.factor = 1.3,
                   stroke = 0,
                   alpha = 1,image.alpha = 0,
                   ncol = 6,
                   combine = T)

SpatialFeaturePlot(st.f,images = c('P5L7','P5T'),
                   features = c('EPCAM','CD3D','ACTA2','MS4A1','MZB1','LYZ','PECAM1','TPSAB1','KLRB1'),
                   pt.size.factor = 1.3,
                   stroke = 0,
                   alpha = 1,image.alpha = 0,
                   ncol = 6,
                   combine = T)

SpatialFeaturePlot(st.f,images = c('P5L7','P5T'),
                   features = c('EPCAM','CD3D','ACTA2','MS4A1','MZB1','LYZ','PECAM1','TPSAB1','KLRB1'),
                   pt.size.factor = 1.3,
                   stroke = 0,
                   alpha = 1,image.alpha = 0,
                   ncol = 6,
                   combine = T)



#####neighborhood analysis, update 2024-08-28#### 



st.i = st.f[,st.f$orig.ident == 'P3T']
spots.info = st.f@images[['P3T']]@coordinates
colnames(spots.info)=c('in_tissue','array_row','array_col','pxl_col','pxl_row')
spots.info$pairxy=paste(spots.info$array_row,spots.info$array_col)

n.spot = ncol(st.i)
#spot.matrix = matrix(0,nrow = n.spot,ncol = n.spot,dimnames = list(colnames(st.i),colnames(st.i)))


spot.matrix = Reduce(cbind,lapply(colnames(st.i),function(a){
  
  coord.i = spots.info[a,]
  row.i = as.integer(coord.i['array_row'])
  col.i = as.integer(coord.i['array_col'])
  candi.coord.i=data.frame(rowc=c(row.i,row.i,row.i-1,row.i-1,row.i+1,row.i+1),
                           colc=c(col.i-2,col.i+2,col.i-1,col.i+1,col.i-1,col.i+1))
  candi.coord.i$pairxy = paste(candi.coord.i$rowc,candi.coord.i$colc)
  
  neis.a = rownames(spots.info[spots.info$pairxy %in% as.character(candi.coord.i$pairxy),])
  
  one.row = data.frame(row.names = colnames(st.i),V1=rep(0,ncol(st.i)))
  one.row[neis.a,'V1'] = 1
  
  return(one.row)
  
}))

spot.matrix = as.matrix(spot.matrix)

clusters = unique(st.f$celltype_prediction)
clusters = clusters[!is.na(clusters)]
meta.data = st.i@meta.data
cluster.matrix = Reduce(cbind,lapply(clusters,function(a){
  
  cluster.spots = rownames(meta.data[meta.data$celltype_prediction == a,])
  cluster.spots = cluster.spots[cluster.spots != 'NA']
  one.row = data.frame(row.names = colnames(st.i),V1=rep(0,ncol(st.i)))
  one.row[cluster.spots,'V1'] = 1
  colnames(one.row) = a
  
  return(one.row)
  
}))



cluster.matrix = as.matrix(cluster.matrix)

CC.nei.matrix = t(cluster.matrix) %*% spot.matrix %*% cluster.matrix

diag(CC.nei.matrix) = 0.5*diag(CC.nei.matrix)

