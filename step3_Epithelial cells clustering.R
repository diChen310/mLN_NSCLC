
.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')
####Epithelial cells ####
library(Seurat)
library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggpubr)
library(cowplot)

####1. Read previous results####
cols.default = c('#e69f00','#54b6e9','#009e73','#f0e442','#0072b2','#d55e00',
                 '#cc79a7','#666666','#ad7700','#1c91d4','#007756','#d5c711',
                 '#005685', '#a04700','#b14380','#4d4d4d','#ffbe2d','#80c7ef',
                 '#00f6b3','#f4eb71')

cds = readRDS(file = './variables_v2/inter_cds.Rds')

epi.data = cds[,cds[['majorCellType']]=='Epi']

table(epi.data[['From']])
# mLN+ N mLN+ PT     mLN  mLN- N mLN- PT     nLN 
# 7693   72332   27457   16484   57478       1 
epi.data = epi.data[,epi.data[['From']]!='nLN'] #remove only one nLN epithelial cell 19094 169124

rm(cds)
####2. Re-clustering####
epi.data <- preprocess_cds(epi.data, num_dim = 50,method="PCA",norm_method="log") 
plot_pc_variance_explained(epi.data)

epi.data <- align_cds(epi.data, num_dim = 20, alignment_group = "PatientID")
epi.data <- reduce_dimension(epi.data,cores = 20)
epi.data = cluster_cells(epi.data, resolution=1e-5)
cols.cl = colorRampPalette(cols.default)(length(unique(monocle3::clusters(epi.data))))

plot_cells(epi.data,group_label_size = 4,cell_size = 0.1,cell_stroke = 0)+ggplot2::scale_color_manual(values = cols.cl)

epi.data[['cluster.epi']]= paste0('Epi_C',clusters(epi.data))

marker_test_res = top_markers(epi.data, group_cells_by="cluster.epi", genes_to_test_per_group=100,
                             reference_cells=1000, cores=32)
write.csv(marker_test_res,file = './tables_v2/Epi assigned cell type marker.csv')
marker_test_res$cell_group = paste(marker_test_res$cell_group,'M')

marker_test_res.epi = read.csv(file = './tables_v2/Epi assigned cell type marker.csv',row.names = 1)
#####Check QC####

meta.data.show = epi.data@colData
meta.data.show = as.data.frame(meta.data.show)

# 
ggplot2::ggplot(meta.data.show,ggplot2::aes(x=cluster.epi,y=nFeature_RNA,fill=cluster.epi))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')

ggplot2::ggplot(meta.data.show,ggplot2::aes(x=cluster.epi,y=percent_mito,fill=cluster.epi))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')


ggplot2::ggplot(meta.data.show,ggplot2::aes(x=cluster.epi,y=scS,fill=cluster.epi))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')

ggplot2::ggplot(meta.data.show,ggplot2::aes(x=orig.ident,fill=cluster.epi))+
  ggplot2::geom_bar(position = 'fill')+facet_wrap(~Source,ncol = 1,scales = 'free')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8))+
  scale_fill_manual(values = cols.cl)



colData(epi.data)$assigned_cell_type <- dplyr::recode(colData(epi.data)$cluster.epi,
                                                      "Epi_C1"="Basal_C1",#KRT6A
                                                      "Epi_C2"="Basal_C2",
                                                      "Epi_C3"="AT2_C1",#WFDC2,SFTA2,SFTPB
                                                      "Epi_C4"="Basal_C3",
                                                      'Epi_C5'='AT2_C2',
                                                      'Epi_C6'='AT2_C3',
                                                      "Epi_C7"="Ciliated",##CAPS,TPPP3
                                                      "Epi_C8"="AT2_C4",
                                                      "Epi_C9"="Basal_C4",#KRT6A, KRT5
                                                      "Epi_C10"="AT2_C5",
                                                      "Epi_C11"="AT2_C6",
                                                      "Epi_C12"="PIP+ Epi",
                                                      'Epi_C13'='Club',#SCGB3A1,SCGB1A1
                                                      "Epi_C14"="AT2_C7",#??ALB
                                                      "Epi_C15"="AT2_C8",
                                                      "Epi_C16"="SUCNR1+ Epi",
                                                      "Epi_C17"="AT1",#AGER, EMP2
                                                      "Epi_C18"="Remove_lowRNA")



epi.data[['subCellType']] = unlist(lapply(as.character(colData(epi.data)$assigned_cell_type),function(a){
  if(grepl('_',a)){
    return(strsplit(a,'_')[[1]][1])
  }else{
    return(a)
  }
}))

cols.cell = cols.default[1:length(unique(SummarizedExperiment::colData(epi.data)$subCellType))]
names(cols.cell)=unique(SummarizedExperiment::colData(epi.data)$subCellType)[order(unique(SummarizedExperiment::colData(epi.data)$subCellType))]
plot_cells(epi.data,group_cells_by = 'subCellType',color_cells_by = 'subCellType',
           cell_stroke = 0,cell_size = 0.1,label_cell_groups = F,
           rasterize = T,
           group_label_size = 3)+scale_color_manual(values = cols.cell )+theme_void()



cols.diseases = jcolors('pal5')[1:4]
names(cols.diseases)=unique(epi.data[['Disease']])
plot_cells(epi.data, color_cells_by="Disease", label_cell_groups=FALSE,cell_size = 0.1,cell_stroke = 0)+
  facet_wrap(~Disease,ncol = 4)+
  scale_color_manual(values = cols.diseases)+
  theme(legend.position = 'none')+
  theme_void()


cols.From = as.character(jcolors('pal6')[1:6])
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')
plot_cells(epi.data,  color_cells_by="From",
           cell_size=0.2,
           cell_stroke=0,
           rasterize = T,
           label_cell_groups = F)+
  ggplot2::scale_color_manual(values = cols.From)+
  theme_void()
saveRDS(epi.data,file = './variables_v2/epi.data.Rds')


table(epi.data[['sorting']],epi.data[['cluster.epi']])

table(epi.data[['orig.ident']],epi.data[['cluster.epi']])

table(epi.data[['From']],epi.data[['cluster.epi']])

epi.data[['From.2']] = ifelse(epi.data[['From']] %in% c('mLN','mLN+ PT'), 'mLN+','mLN-')
marker_test_res.From.2 = top_markers(epi.data, group_cells_by="From.2", genes_to_test_per_group=200,
                              reference_cells=1000, cores=32)

