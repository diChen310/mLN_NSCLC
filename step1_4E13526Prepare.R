####Step1, E-MTAB-13526 scRNA-seq data####
.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')
library(monocle3)
library(Seurat)
library(jcolors)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggsci)
library(pheatmap)
library(reshape2)
library(ggpubr)

sc.data = readRDS(file = '/share/dichen/lungNodeM/variables/E13526Data.Rds')

sc.data[['majorCellType']] = dplyr::recode(sc.data[['Cell.types.v25']],
                                           'Cycling anti-inflammatory macrophages'='Myeloid',
                                           'Cycling plasma B cells'= 'B',                
                                           'Cycling T cells'='T',
                                           'mo-DC2'='Myeloid',                                
                                           'Cycling exhausted cytotoxic T cells'='T',
                                           'Cycling myeloid cells,2,0 (to remove)'='Unassigned', 
                                           'TNF+ B cells'='B',
                                           'LYZ+ B cells'='B',                          
                                           'Downregulated B cells'='B',
                                           'Plasma B cells'='Plasma',                        
                                           'B cells'='B',
                                           'Mast cells'='Mast',
                                           'Mast cells,3 (to remove)'='Unassigned',
                                           'Cycling mast cells'='Mast',                    
                                           'Mast cells,9 (to remove)'='Unassigned',
                                           'Pro-inflammatory macrophages'='Myeloid',          
                                           'Immature myeloid cells'='Myeloid',
                                           'STAB1+ anti-inflammatory macrophages'='Myeloid',  
                                           'cDC2'='Myeloid',                                    
                                           'Anti-inflammatory macrophages'='Myeloid',          
                                           'Anti-inflammatory alveolar macrophages'='Myeloid',  
                                           'CAMLs'='Myeloid',                                  
                                           'NK cells (higher cytotoxicity)'='T',
                                           'Cytotoxic T cells' = 'T',                    
                                           'NK cells (lower cytotoxicity)' = 'T',         
                                           'Exhausted cytotoxic T cells' ='T',           
                                           'AT2 cells' ='Epithelial',                             
                                           'Transitioning epithelial cells'  ='Epithelial',          
                                           'Atypical epithelial cells'  ='Epithelial',                
                                           'Ciliated epithelial cells'  ='Epithelial',               
                                           'Activated adventitial fibroblasts' ='Fibro',     
                                           'Lymphatic endothelial cells'    ='Endo',       
                                           'Fibroblasts'='Fibro',                           
                                           'Cycling AT2 cells'   ='Epithelial',                      
                                           'Cycling epithelial cells'    ='Epithelial',               
                                           'Naive T cells' ='T',                        
                                           'Downregulated T cells'  ='T',             
                                           'Tregs'='T',                                 
                                           'Exhausted T cells' ='T' ,                    
                                           'pDCs'    = 'Myeloid'                              )


meta.data = as.data.frame(colData(sc.data))

ggplot(meta.data,aes(x=patient_sample,y=n_genes))+geom_boxplot()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))

ggplot(meta.data,aes(x=patient_sample,y=doublet_scores))+geom_boxplot()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(meta.data,aes(x=patient_sample,fill=Cell.types.v25))+geom_bar(position = 'fill')+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))

ggplot(meta.data,aes(x=Cell.types.v25,y=doublet_scores))+geom_boxplot()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))

table(meta.data$patient_sample)


ggplot(meta.data,aes(x=patient_sample,fill=majorCellType))+geom_bar(position = 'fill')+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))
ggplot(meta.data,aes(x='',fill=majorCellType))+geom_bar(position = 'fill')+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))


N0.patients = c('Patient 4','Patient 6','Patient 7','Patient 10','Patient 11',
                'Patient 12','Patient 13','Patient 14','Patient 15','Patient 16',
                'Patient 17','Patient 19','Patient 21','Patient 22','Patient 23',
                'Patient 24','Patient 25')
remove.patient = 'Patient 9' ###N0/N1

sc.data[['N0']] = ifelse(sc.data[['patient']] %in% N0.patients,'N0','Not N0')


sc.data.seurat = CreateSeuratObject(counts = assay(sc.data))
sc.data.seurat = AddMetaData(sc.data.seurat,as.data.frame(colData(sc.data)))

VlnPlot(sc.data.seurat,features = c('FTH1','FTL','CST1'),group.by = 'N0',
        pt.size=0)
VlnPlot(sc.data.seurat[,sc.data.seurat$majorCellType == 'Epithelial'],features = c('FTH1','FTL','CLDN4'),group.by = 'N0',
        pt.size=0)


sc.data.seurat = sc.data.seurat[,sc.data.seurat$patient != remove.patient]
