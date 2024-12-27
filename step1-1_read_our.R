####QC for all scRNA-seq objects, 2024-09-05
.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')
library(ggplot2)
library(Seurat)
library(dplyr)
library(scDblFinder)
###sc.objs from step0
####Find doubles####
for(i in 1:length(sc.objs)){
  
  sc.obj = sc.objs[[i]]
  db.r = scDblFinder(as.SingleCellExperiment(sc.obj))
  sc.objs[[i]]$scDblFinder.class = db.r$scDblFinder.class
  sc.objs[[i]]$scS = db.r$scDblFinder.score
  
  print(paste('scDdlFinder',cell_ids[i],'Finished!'))
  
}
print(Sys.time())

#2024-09-13 13:43 to 14:16:20
####Change names####
sc.names.info = read.csv(file = './documents/scRNA names.csv',row.names = 1)
sc.names.info$New_ident[sc.names.info$New_ident == 'Pl5.T'] = 'P15.T'
sc.names.info$New_ident[sc.names.info$New_ident == 'Pl6.L'] = 'P16.L'

for(i in 1:length(sc.objs)){
  
  sc.objs[[i]]$orig.ident = sc.names.info[as.character(sc.objs[[i]]$orig.ident),'New_ident']
  print(unique(sc.objs[[i]]$orig.ident))
}

for(i in 1:length(sc.objs)){
  
  print(sc.objs[[i]]$orig.ident)
}
saveRDS(sc.objs,file='./variables_v2/sc.objs.withDbScores.Rds')#"2024-09-13 14:34:31 +08"

####Merge into one####
cell_ids = as.character(unlist(lapply(1:length(sc.objs),function(a){
  unique(sc.objs[a][[1]]$orig.ident)
})))
sc.m = merge(sc.objs[1][[1]],sc.objs[-1],add.cell.ids = cell_ids)
print(unique(sc.m$orig.ident))

####Calculate mito and hb percent####

sc.m = PercentageFeatureSet(sc.m, "^HB[^(P)]", col.name = "percent_hb")
sc.m = PercentageFeatureSet(sc.m, "^MT-", col.name = "percent_mito")

Idents(sc.m)=sc.m$orig.ident

#### Plot QC for all samples####

P1=VlnPlot(sc.m, features = "nCount_RNA", pt.size = 0,log = T) + NoLegend()
P2=VlnPlot(sc.m, features = "nFeature_RNA", pt.size = 0,log = T) + NoLegend()
P3=VlnPlot(sc.m, features = "percent_mito", pt.size = 0,sort = 'increasing')+ NoLegend()
P4=VlnPlot(sc.m, features = "percent_hb", pt.size = 0)+ NoLegend()

P1 / P2 
P3 / P4

P3+ggplot2::geom_hline(yintercept = c(15,20,25,30),lty=3,color='red')


sc.m = sc.m[,sc.m$orig.ident %in% c('P0.L10','P0.3A') == F] ###remove poor or un-clear samples


####Filter the scRNA object one by one####
####
# Quality control of cells was performed using scater (v1.14.6) whereby cells with total read counts or number of genes 
# detected greater than 3 median absolute deviations (MADs) below the median (both calculated on the log10-scale) 
# or cells where the percentage of mitochondrial counts was 3 MADs above the median were removed from the filtered matrices downloaded from the 10x website.
#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02552-3



# sc.objs.filter = c()
# cell_ids.filter = c()
# 
# for(i in 1:length(sc.objs)){
#   cell_id = cell_ids[i]
#   if(cell_id %in% c('P0.L10','P0.3A') ==F){
#     sc.obj = sc.objs[[i]]
#     
#     print(paste(c(cell_id,'Raw size',dim(sc.obj)),collapse = ','))
#     
#     mt.th = 25
#     nFeature.th = 300
#     sc.obj = PercentageFeatureSet(sc.obj, "^MT-", col.name = "percent_mito")
#     
#     sc.obj = sc.obj[,sc.obj$percent_mito < mt.th]
#     
#     sc.obj = sc.obj[,sc.obj$nFeature_RNA > nFeature.th]
#     
#     print(paste(c(cell_id,'Filter by MT and nFeature(',mt.th,nFeature.th,')',dim(sc.obj)),collapse = ','))
#     
#     
#     #sc.obj = sc.obj[,sc.obj$scDblFinder.class == 'singlet']
#     
#     print(paste(c(cell_id,'Filter doublets',dim(sc.obj)),collapse = ','))
#     
#     cell_ids.filter = append(cell_ids.filter,cell_id)
#     sc.objs.filter = append(sc.objs.filter,sc.obj)
#     
#   }
#   
# }
# 
# rm(sc.m)
# gc()


saveRDS(sc.objs.filter,file='./variables_v2/sc.objs.filter.Rds')




