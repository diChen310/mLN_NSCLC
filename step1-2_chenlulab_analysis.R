#####Read chenlulab data####2024-0905
.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')


library(Seurat)
library(dplyr)
library(scDblFinder)
library(circlize)
library(RColorBrewer)
library(jcolors)
library(celda)
library(harmony)
library(corrplot)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(ggsci)
library(clusterProfiler)
library(ggplot2)
library(cowplot)
library(data.table)

meta.data = read.csv(file = './data/meta.csv')
meta.data[meta.data$SampleID %in% c('BT1247'),'Disease']='LUSC_Normal' 
meta.data[meta.data$Stage == 'I','N']='N0'
meta.data[meta.data$Stage == 'I','M']='M0'

sc.RNA = fread(file = './data/raw_data.csv',sep = ',')
sc.RNA = as.data.frame(sc.RNA)
rownames(sc.RNA)=sc.RNA$V1
sc.RNA[1:5,1:5]
sc.RNA = sc.RNA[,-1]

sc.obj = CreateSeuratObject(counts = as.sparse(sc.RNA),project = "chenlu")


meta.data.u = meta.data[!duplicated(meta.data$SampleID),]
table(meta.data.u$Tissue,meta.data.u$Disease)
meta.data.u.sel = meta.data.u[meta.data.u$Disease %in% c('LUAD','LUSC','LUAD_Normal','LUSC_Normal'),]
meta.data.u.sel = meta.data.u.sel[meta.data.u.sel$N != '',]
meta.data.u.sel = meta.data.u.sel[meta.data.u.sel$M == 'M0',] ##88
  

sel.cells = meta.data[meta.data$SampleID %in% unique(meta.data.u.sel$SampleID),'Cells']
rownames(meta.data) = meta.data$Cells
sel.cells.meta.data = meta.data[sel.cells,]

sc.obj = sc.obj[,sel.cells]
sc.obj=AddMetaData(sc.obj,sel.cells.meta.data)

rm(sc.RNA)
gc()

sc.obj = PercentageFeatureSet(sc.obj, "^MT-", col.name = "percent_mito")



gc()
Idents(sc.obj)=sc.obj$SampleID
VlnPlot(sc.obj, features = "nCount_RNA", pt.size = 0,log = T) + NoLegend()
VlnPlot(sc.obj, features = "nFeature_RNA", pt.size = 0,log = T) + NoLegend()
VlnPlot(sc.obj, features = "percent_mito", pt.size = 0,sort = 'increasing')+ NoLegend()

meta.data.sel = sc.obj@meta.data

saveRDS(sc.obj,file = './variables_v2/chenluData.Rds')

