####Read all scRNA objects from all five folders####
.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')
library(vctrs)
library(Seurat)
library(dplyr)
#library(scDblFinder)

sc.objs = c()
cell_ids = c()
sc.input.folder = './result/'
sc.folders = list.dirs(path = sc.input.folder,recursive = F)
for(sc.f in sc.folders[1:3]){
  
  names = strsplit(sc.f,'//')[[1]][2]
  names = strsplit(names,'_|-')[[1]]
  names = names[names %in% c('RNA','result') == F]
  sc.obj = read.delim(file=gzfile(paste0(sc.f,'/expression_matrix/',paste(names,collapse = '_'),'_matrix.tsv.gz')))
  sc.obj = CreateSeuratObject(sc.obj, names.field = 2,min.cells = 3, project = paste(names,collapse = ''))
  
  print(paste('Read',sc.f,'Finished!'))
  cell_ids = append(cell_ids,paste(names,collapse = ''))
  sc.objs = append(sc.objs,sc.obj)
}
for(sc.f in sc.folders[4:11]){
  
  names = strsplit(sc.f,'//')[[1]][2]
  names = strsplit(names,'_|-')[[1]]
  names = names[names %in% c('RNA','result') == F]
  #sc.obj = read.delim(file=gzfile(paste0(sc.f,'/expression_matrix/',paste(names,collapse = '_'),'_matrix.tsv.gz')))
  
  sc.obj = read.delim(file=gzfile(paste0(sc.f,'/expression_matrix/',paste(names,collapse = '-'),'_matrix.tsv.gz')))
  sc.obj = CreateSeuratObject(sc.obj, min.cells = 3, project = paste(names,collapse = ''))
  
  print(paste('Read',sc.f,'Finished!'))
  cell_ids = append(cell_ids,paste(names,collapse = ''))
  sc.objs = append(sc.objs,sc.obj)
}
for(sc.f in sc.folders[12:14]){
  
  names = strsplit(sc.f,'//')[[1]][2]
  names = strsplit(names,'_|-')[[1]]
  names = names[names %in% c('L','S') == F]
  sc.obj = Read10X(data.dir = paste0(sc.f,'/filtered_feature_bc_matrix'))
  sc.obj = CreateSeuratObject(sc.obj, min.cells = 3, project = paste(names,collapse = ''))
  
  print(paste('Read',sc.f,'Finished!'))
  cell_ids = append(cell_ids,paste(names,collapse = ''))
  sc.objs = append(sc.objs,sc.obj)
}
sc.folders.2 = list.dirs(path = './new/C1656909322225315840/',recursive = F)
for(sc.f in sc.folders.2[1:4]){
  names = strsplit(sc.f,'//')[[1]][2]
  names = strsplit(names,'_|-')[[1]]
  names = names[names %in% c('L','S') == F]
  
  sc.obj = Read10X(data.dir = paste0(sc.f,'/filtered_feature_bc_matrix'))
  sc.obj = CreateSeuratObject(sc.obj, min.cells = 3, project = paste(names,collapse = ''))
  
  print(paste('Read',sc.f,'Finished!'))
  cell_ids = append(cell_ids,paste(names,collapse = ''))
  sc.objs = append(sc.objs,sc.obj)
  
}
for(sc.f in sc.folders.2[5:6]){
  names = strsplit(sc.f,'//')[[1]][2]
  names = strsplit(names,'_|-')[[1]]
  names = names[names %in% c('L','S') == F]
  
  sc.obj = Read10X(data.dir = paste0(sc.f,'/05.count/',paste(names,collapse = '-'),'_filtered_feature_bc_matrix'))
  sc.obj = CreateSeuratObject(sc.obj, min.cells = 3, project = paste(names,collapse = '-'))
  
  print(paste('Read',sc.f,'Finished!'))
  cell_ids = append(cell_ids,paste(names,collapse = ''))
  sc.objs = append(sc.objs,sc.obj)
}
sc.folders.3 = list.dirs('./new/mengdewei0628082/X101SC22064646-Z01-F001_release_20230628/Result-X101SC22064646-Z01-J001-B1-1/2.Celescope/',recursive = F)
for(sc.f in sc.folders.3){
  names = strsplit(sc.f,'//')[[1]][2]
  #names = strsplit(names,'_|-')[[1]]
  names = names[names %in% c('L','S') == F]
  
  sc.obj = Read10X(data.dir = paste0(sc.f,'/',names,'_filtered_feature_bc_matrix'))
  sc.obj = CreateSeuratObject(sc.obj, min.cells = 3, project = paste(names,collapse = '-'))
  
  print(paste('Read',sc.f,'Finished!'))
  cell_ids = append(cell_ids,paste(names,collapse = ''))
  sc.objs = append(sc.objs,sc.obj)
}

sc.folders.4=c('./new2/朴海龙，刘宇-6个样本/表达矩阵-6例样本/P23090703//MLX622562T/20231013/Matrix',
               './new2/朴海龙，刘宇-6个样本/表达矩阵-6例样本/P23090703//MLX622562LN5/20231013/Matrix',
               "./new2/朴海龙，刘宇-6个样本/表达矩阵-6例样本/P23090703//TDM-m/20231013/Matrix" ,
               "./new2/朴海龙，刘宇-6个样本/表达矩阵-6例样本/P23090703//TDM-ln/20231013/Matrix",
               "./new2/朴海龙，刘宇-6个样本/表达矩阵-6例样本/P23090703//ZCD-Ca/20231013/Matrix",
               "./new2/朴海龙，刘宇-6个样本/表达矩阵-6例样本/P23090703//ZCD-Ln/20231013/Matrix")

for(sc.f in sc.folders.4){
  
  names = strsplit(sc.f,'//')[[1]][2]
  names = strsplit(names,'/')[[1]][1]
  
  sc.obj = Read10X(data.dir = sc.f)
  sc.obj = CreateSeuratObject(sc.obj, min.cells = 3, project = paste(names,collapse = ''))
  
  print(paste('Read',sc.f,'Finished!'))
  cell_ids = append(cell_ids,paste(names,collapse = ''))
  sc.objs = append(sc.objs,sc.obj)
}


sc.folders.5 = c('/share/dichen/lungNodeM/new2/刘宇-4个转录组/countMatrix/GLR-T_filtered_feature_bc_matrix/',
                 '/share/dichen/lungNodeM/new2/刘宇-4个转录组/countMatrix/GLR-LN_filtered_feature_bc_matrix/',
                 '/share/dichen/lungNodeM/new2/刘宇-4个转录组/countMatrix/BJ23002985-WW-T_filtered_feature_bc_matrix/',
                 '/share/dichen/lungNodeM/new2/刘宇-4个转录组/countMatrix/BJ23002986-WW-N_filtered_feature_bc_matrix/')
names.5=c('GLR_T','GLR_LN','WW_T','WW_LN')
for(i in 1:4){
  sc.f = sc.folders.5[i]
  names = names.5[i]
  
  sc.obj = Read10X(data.dir = sc.f)
  sc.obj = CreateSeuratObject(sc.obj, min.cells = 3, project = names)
  
  print(paste('Read',sc.f,'Finished!'))
  cell_ids = append(cell_ids,names)
  sc.objs = append(sc.objs,sc.obj)
}
print(Sys.time())#2024-09-05 9:53 to 10:06

#2024-09-13 13:31 to 13:43

###L10, L3A poor
