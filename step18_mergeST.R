.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')
library(Seurat)
library(SeuratData)
####Integrate all samples####
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



st.folder = '/share/dichen/lungST/data/E-MTAB-13530/h5 files/'
img.folder = '/share/dichen/lungST/data/E-MTAB-13530/'
st.files = list.files(path = st.folder)
st.files = st.files[c(-2,-3,-4,-6,-7,-8,-9,-12,-13,-16,-18,-19,-21,-22,-23,-24,-25,-26,-27,-28)]

for(st.file.i in st.files){
  st.file.i = st.files[3]
  st.i = Read10X_h5(filename = paste0(st.folder,st.file.i))
  name = strsplit(st.file.i,'-')[[1]][1]
  st.i = CreateSeuratObject(st.i, assay = 'Spatial')
  st.i$slice = 1
  st.i$region = name
  st.i$From = 'Ana'
  img = Seurat::Read10X_Image(image.dir = paste0(img.folder,name,'-spatial'))####Manully change
  Seurat::DefaultAssay(object = img) = 'Spatial'
  img = img[colnames(x = st.i)]
  st.i[['image']] = img
  st.objs = append(st.objs,st.i)
  print(paste('Read',name,'Finished!'))
  cell_ids = append(cell_ids,name)
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

meta.data = st.m@meta.data
meta.data[meta.data$orig.ident == 'P3L12','region'] = 'P3L'
meta.data[meta.data$orig.ident == 'P3T','region'] = 'P3T'
meta.data[meta.data$orig.ident == 'P7L5','region'] = 'P7L'
meta.data[meta.data$orig.ident == 'P7T','region'] = 'P7T'
meta.data[meta.data$orig.ident == 'P5L','region'] = 'P5L'
meta.data[meta.data$orig.ident == 'P5T','region'] = 'P5T'
meta.data[meta.data$orig.ident == 'P6L','region'] = 'P6L'
meta.data[meta.data$orig.ident == 'P6T','region'] = 'P6T'
meta.data$orig.ident = meta.data$region
meta.data$Source = ifelse(substr(meta.data$orig.ident,1,2) %in% c('P3','P5','P6','P7'),'PLM','Ana')

meta.data$From = ''

meta.data[substr(meta.data$orig.ident,1,2) %in% c('P3','P5','P6','P7') == F,'From'] = 'mLN- PT'
meta.data[meta.data$orig.ident %in% c('P3L','P5L','P6L','P7L') ,'From'] = 'mLN'
meta.data[meta.data$orig.ident %in% c('P3T','P5T','P6T','P7T') ,'From'] = 'mLN+ PT'

st.m$region = meta.data$region
st.m$orig.ident = meta.data$orig.ident
st.m$Source = meta.data$Source
st.m$From = meta.data$From
SpatialFeaturePlot(st.m,features = 'CST1' ,ncol = 4)

VlnPlot(st.m,features = c('CST1','IGHG4','CXCL13'),group.by = 'From',pt.size = 0,adjust = 19)
saveRDS(st.m,file = '/share/dichen/lungNodeM/st_E13530&Our.Rds')