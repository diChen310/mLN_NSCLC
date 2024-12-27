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
cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559",'#5d94b4',
                 '#b8a3c4','#60428d','#c79795','#c92e2a','#9f5633','#d9c17c',"#d0896b","black",'grey',"blue",'cyan',"#32709f")

cds.all = readRDS(file = './variables_v2/inter_cds.Rds')
stromal.data = cds.all[,cds.all[['majorCellType']] %in% c('Fibro','Endo')]

rm(cds.all)

gc()

####Clustering####
stromal.data <- preprocess_cds(stromal.data, num_dim = 50,method="PCA",norm_method="log") 

stromal.data <- align_cds(stromal.data, num_dim = 20, alignment_group = "PatientID")
stromal.data <- reduce_dimension(stromal.data,cores = 20)
stromal.data = cluster_cells(stromal.data, resolution=1e-4)

cols.cl = colorRampPalette(cols.default)(length(unique(monocle3::clusters(stromal.data))))
stromal.data[['cluster.stromal']]= paste0('Stromal_C',clusters(stromal.data))
plot_cells(stromal.data,group_label_size = 4,cell_size = 0.3,cell_stroke = 0)+
  ggplot2::scale_color_manual(values = cols.cl)

marker_test_res = top_markers(stromal.data,group_cells_by="cluster.stromal", genes_to_test_per_group=100,
                              reference_cells=1000, cores=32)
write.csv(marker_test_res,file = './tables_v2/Stromal cluster marker.csv')


stromal.data[['assigned_cell_type']] = recode(stromal.data[['cluster.stromal']],
                                              "Stromal_C1"="mCAF_C1",
                                              "Stromal_C2"="Veins",#
                                              "Stromal_C3"="iCAF",
                                              "Stromal_C4"="vCAF",#NOTCH3, COL18A1
                                              "Stromal_C5"="Arteries",#HEY1,IGFBP3
                                              "Stromal_C6"="Remove_twoCellTypes",#Endo & Fibro
                                              "Stromal_C7"="LEC",
                                              'Stromal_C8'='mCAF_C2',#CST1
                                              "Stromal_C9"="Remove_twoCellTypes",#
                                              "Stromal_C10"="Capillaries",#CA4
                                              'Stromal_C11'="Remove_twoCellTypes",#
                                              'Stromal_C12'='tCAF',
                                              'Stromal_C13'='Remove_twoCellTypes',
                                              'Stromal_C14'='Pericyte',#RGS5,NOTCH3,COX4I2
                                              "Stromal_C15"="Remove_twoCellTypes",
                                              'Stromal_C16'="Remove_twoCellTypes",#
                                              'Stromal_C17'='Proliferation stromal',
                                              'Stromal_C18'='Remove_patientSpecific',
                                              'Stromal_C19'='Remove_lowNumber')


stromal.data[['subCellType']] = unlist(lapply(as.character(colData(stromal.data)$assigned_cell_type), function(a){
  if(grepl('_',a)){return( strsplit(a,'_')[[1]][1])}else{return(a)}
}))

saveRDS(stromal.data,file = './variables_v2/stromal.data.Rds')


####Re-clustering after removing Remove####

stromal.data = stromal.data[,stromal.data[['subCellType']] != 'Remove'] #14787 15595

stromal.data <- preprocess_cds(stromal.data, num_dim = 50,method="PCA",norm_method="log") 

stromal.data <- align_cds(stromal.data, num_dim = 20, alignment_group = "PatientID")
stromal.data <- reduce_dimension(stromal.data,cores = 20)
stromal.data = cluster_cells(stromal.data, resolution=1e-3)

cols.cl = colorRampPalette(cols.default)(length(unique(monocle3::clusters(stromal.data))))
stromal.data[['cluster.stromal.v2']]= paste0('Stromal_C',clusters(stromal.data))
plot_cells(stromal.data,group_label_size = 4,cell_size = 1,cell_stroke = 0)+
  ggplot2::scale_color_manual(values = cols.cl)

marker_test_res.r2 = top_markers(stromal.data,group_cells_by="cluster.stromal.v2", genes_to_test_per_group=100,
                              reference_cells=1000, cores=32)
write.csv(marker_test_res.r2,file = './tables_v2/Stromal cluster marker.r2.csv')

meta.data.show = stromal.data@colData
meta.data.show = as.data.frame(meta.data.show)
meta.data.show$Disease2 = substr(meta.data.show$Disease,1,4)


ggplot2::ggplot(meta.data.show,ggplot2::aes(x=From,fill=cluster.stromal.v2))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+facet_wrap(~Disease2,ncol = 1)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values =  cols.cl)+
  theme(axis.text  = element_text(size = 8),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))


marker_test_res.r2$cell_group = paste(marker_test_res.r2$cell_group,'M')

stromal.data[['assigned_cell_type']] = recode(stromal.data[['cluster.stromal.v2']],
                                              "Stromal_C1"="Veins_C1",#ACKR1
                                              "Stromal_C2"="iCAF_C1",#C3,CFD
                                              "Stromal_C3"="mCAF_C1",#COL1A1
                                              "Stromal_C4"="vCAF",#NOTCH3, COL18A1
                                              "Stromal_C5"= 'Capillaries',#CA4
                                              "Stromal_C6"="Veins_C2",#ACKR1
                                              "Stromal_C7"="mCAF_C2",
                                              'Stromal_C8'='Arteries_C1',#HEY1, IGFBP3
                                              "Stromal_C9"="mCAF_C3",#
                                              "Stromal_C10"="LEC_C1",#CCL21,TFF3,PROX1
                                              'Stromal_C11'="Remove_twoCellTypes",#T cell markers
                                              'Stromal_C12'='Veins_C3',
                                              'Stromal_C13'='SMC',#ACTA2,MYH11
                                              'Stromal_C14'='mCAF_C4',#
                                              "Stromal_C15"="mCAF_C5",
                                              'Stromal_C16'="iCAF_C2",#
                                              'Stromal_C17'='mCAF_C6',
                                              'Stromal_C18'='Veins_C4',
                                              'Stromal_C19'='Endo progenitor',#CD34, KDR
                                              'Stromal_C20'="Remove_twoCellTypes",#T cell markers
                                              'Stromal_C21'="Aerocyte",#CA4
                                              'Stromal_C22'='Remove_twoCellTypes',#B
                                              'Stromal_C23'='tCAF',
                                              'Stromal_C24'='Remove_twoCellTypes',#Epi, myeloid
                                              'Stromal_C25'='Proliferation stromal',#Epi
                                              'Stromal_C26'='Pericyte',#RGS5,NOTCH3,COX4I2
                                              'Stromal_C27'='CST1+ Fibro',
                                              'Stromal_C28'='Remove_twoCellTypes',#Epi
                                              'Stromal_C29'='Remove_twoCellTypes',#T
                                              'Stromal_C30'='Remove_twoCellTypes',#Epi
                                              'Stromal_C31'='Remove_twoCellTypes'#Epi
                                              )


stromal.data[['subCellType']] = unlist(lapply(as.character(colData(stromal.data)$assigned_cell_type), function(a){
  if(grepl('_',a)){return( strsplit(a,'_')[[1]][1])}else{return(a)}
}))

stromal.data = stromal.data[,stromal.data[['subCellType']] != 'Remove'] #14787 15595

saveRDS(stromal.data,file = './variables_v2/stromal.data.refine.Rds')


