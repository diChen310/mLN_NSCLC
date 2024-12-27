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


####1. Read previous results####
cols.default = c("#36977b","#a0ca8a",'#fce8a9',"#ebad71","#bb8f7c",'#32709f',"#639f95","#fafa64","#e58b8e","#d65559",'#5d94b4',
                 '#c79795','#60428d',"#479544",'#c92e2a','#9f5633','#d9c17c',"#d0896b","black",'grey',"blue",'cyan')

cds.all = readRDS(file = './variables_v2/inter_cds.Rds')

Myeloid.data = cds.all[,cds.all[['majorCellType']] =='Myeloid']

rm(cds.all)

gc()


####2. Re-clustering####
#selected_genes = rownames(Myeloid.data)[ Matrix::rowSums(Myeloid.data@assays@data$counts > 0)> 50]

Myeloid.data <- monocle3::preprocess_cds(Myeloid.data, num_dim = 50,method="PCA",norm_method="log") 
monocle3::plot_pc_variance_explained(Myeloid.data)

Myeloid.data <- monocle3::align_cds(Myeloid.data, num_dim = 20, alignment_group = "PatientID")
Myeloid.data <- monocle3::reduce_dimension(Myeloid.data,cores = 20)
Myeloid.data = monocle3::cluster_cells(Myeloid.data, resolution=1e-5)
# 1e-5, cluster number 24
cols.cl = colorRampPalette(cols.default)(length(unique(monocle3::clusters(Myeloid.data)))) 

monocle3::plot_cells(Myeloid.data,group_label_size = 4,cell_size = 0.1,cell_stroke = 0)+ggplot2::scale_color_manual(values = cols.cl)

Myeloid.data[['cluster.Myeloid']]= paste0('Myeloid_C',monocle3::clusters(Myeloid.data))

monocle3::plot_cells(Myeloid.data,color_cells_by = 'From',group_label_size = 4,cell_size = 0.1,cell_stroke = 0)

#plot_cells(Myeloid.data,color_cells_by = 'partition',group_label_size = 4,cell_size = 0.1,cell_stroke = 0)
plot_cells(Myeloid.data,genes = c('FABP4','SPP1'),cell_size = 0.1,cell_stroke = 0,scale_to_range = F)+
  scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))

plot_cells(Myeloid.data,genes = c('CD1C'),cell_size = 0.1,cell_stroke = 0,scale_to_range = F)+
  scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))


#####Find markers####

marker_test_res = monocle3::top_markers(Myeloid.data,
                              group_cells_by="cluster.Myeloid", genes_to_test_per_group=100,
                              reference_cells=1000, cores=32)
write.csv(marker_test_res,file = './tables_v2/Myeloid assigned cell type marker.csv')
marker_test_res$cell_group = paste(marker_test_res$cell_group,'M')
table(Myeloid.data[['From']],Myeloid.data[['cluster.Myeloid']])


#####Annotations####
######################################CCL21!!!!!!!!!!!!!!! check#########################################
colData(Myeloid.data)$assigned_cell_type = Myeloid.data[['cluster.Myeloid']]
colData(Myeloid.data)$assigned_cell_type <- dplyr::recode(colData(Myeloid.data)$assigned_cell_type,
                                                          
                                                          "Myeloid_C1"="MP_C1",
                                                          'Myeloid_C2'='FABP4+ MP_C1',
                                                          'Myeloid_C3'='SPP1+ MP_C1',
                                                          'Myeloid_C4'='SPP1+ MP_C1',
                                                          'Myeloid_C5'='CD1C+ DC_C1',#CCL17, CD1C, CD1A, CD1B
                                                          "Myeloid_C6"="GPNMB+ MP",
                                                          'Myeloid_C7'='CD1C+ DC_C2',#
                                                          'Myeloid_C8'='CD1C+ DC_C3',
                                                          'Myeloid_C9'='CD1C+ DC_C4',
                                                          'Myeloid_C10'='FCGBP+ MP_C1',#
                                                          'Myeloid_C11'='SPP1+ MP_C2',
                                                          'Myeloid_C12'='CD1C+ DC_C5',
                                                          'Myeloid_C13'='Proliferation MP',
                                                          'Myeloid_C14'='Monocyte_C1',#AIF, CD14
                                                          'Myeloid_C15'='CD1C+ DC_C6',
                                                          'Myeloid_C16'='pDC',
                                                          'Myeloid_C17'='SPP1+ MP_C3',#
                                                          'Myeloid_C18'='CD1C- DC',#
                                                          "Myeloid_C19"="LAMP3+ DC",#
                                                          "Myeloid_C20"="cDC",#
                                                          'Myeloid_C21'='Neu_C1',#
                                                          'Myeloid_C22'='Remove_lowRNA',#
                                                          "Myeloid_C23"="Remove_lowRNA",#
                                                          "Myeloid_C24"="FABP4+ MP_C2"
)
Myeloid.data[['subCellType']] = unlist(lapply(as.character(colData(Myeloid.data)$assigned_cell_type), function(a){
  if(grepl('_',a)){return( strsplit(a,'_')[[1]][1])}else{return(a)}
}))

cols.cell = cols.default[1:length(unique(colData(Myeloid.data)$subCellType))]
#cols.cell[6]='#fbc543'
names(cols.cell)=unique(colData(Myeloid.data)$subCellType)[order(unique(colData(Myeloid.data)$subCellType))]


plot_cells(Myeloid.data,group_cells_by = 'subCellType',color_cells_by = 'subCellType',
              cell_stroke = 0,cell_size = 0.5,label_cell_groups = F,
              group_label_size = 4)+scale_color_manual(values = cols.cell )+
  theme_void()

saveRDS(Myeloid.data,file = './variables_v2/myeloid.data.Rds')



####3. re-analyze the marker genes after remove####

Myeloid.data = readRDS(file = './variables_v2/myeloid.data.Rds') # 14787 278503


Myeloid.data = Myeloid.data[,Myeloid.data[['subCellType']] != 'Remove'] # 14787 276315
plot_cells(Myeloid.data,group_cells_by = 'cluster.Myeloid',color_cells_by = 'cluster.Myeloid',
           cell_stroke = 0,cell_size = 0.5,label_cell_groups = T,
           group_label_size = 4)+
  theme_void()
cc=table(Myeloid.data[['From']],Myeloid.data[['cluster.Myeloid']])/apply(table(Myeloid.data[['From']],Myeloid.data[['cluster.Myeloid']]),1,sum)
heatmap(cc)

marker_test_res.r2 = top_markers(Myeloid.data,group_cells_by="cluster.Myeloid", genes_to_test_per_group=100,
                                 reference_cells=1000, cores=32)
marker_test_res.r2$cell_group = paste(marker_test_res.r2$cell_group,'M')

myeloid.sc = CreateSeuratObject(counts = assay(Myeloid.data))
myeloid.sc = AddMetaData(myeloid.sc,as.data.frame(colData(Myeloid.data)))
Idents(myeloid.sc)=myeloid.sc$cluster.Myeloid
mc6.marker = FindMarkers(myeloid.sc,ident.1 = 'Myeloid_C6',only.pos = T,logfc.threshold = 1)

colData(Myeloid.data)$assigned_cell_type = Myeloid.data[['cluster.Myeloid']]
colData(Myeloid.data)$assigned_cell_type <- dplyr::recode(colData(Myeloid.data)$assigned_cell_type,
                                                          "Myeloid_C1"="cDC2_C1", #HLA-DPA1,CST3 based on cds.all
                                                          'Myeloid_C2'='FABP4+ MP_C1',
                                                          'Myeloid_C3'='SPP1+ MP',
                                                          'Myeloid_C4'='SPP1+ MP',
                                                          'Myeloid_C5'='cDC2_C2',#CCL17, CD1C, CD1A, CD1B,FCER1A
                                                          "Myeloid_C6"="APOE+ MP",#
                                                          'Myeloid_C7'='cDC2_C3',#
                                                          'Myeloid_C8'='cDC2_C4',#CD1C based on cds.all
                                                          'Myeloid_C9'='cDC2_C5',
                                                          'Myeloid_C10'='cDC2_C6',#
                                                          'Myeloid_C11'='Gran_C2',
                                                          'Myeloid_C12'='cDC2_C7',
                                                          'Myeloid_C13'='Proliferation Myeloid',
                                                          'Myeloid_C14'='Gran_C1',#S100A12
                                                          'Myeloid_C15'='cDC2_C8',
                                                          'Myeloid_C16'='pDC',
                                                          'Myeloid_C17'='SPP1+ MP',#SPP1 based on cds.all
                                                          'Myeloid_C18'='CD1C- DC',#
                                                          "Myeloid_C19"="LAMP3+ DC",#
                                                          "Myeloid_C20"="cDC1",#
                                                          'Myeloid_C21'='Neu',#
                                                          "Myeloid_C24"="FABP4+ MP_C2"
)
Myeloid.data[['subCellType']] = unlist(lapply(as.character(colData(Myeloid.data)$assigned_cell_type), function(a){
  if(grepl('_',a)){return( strsplit(a,'_')[[1]][1])}else{return(a)}
}))



saveRDS(Myeloid.data,file = './variables_v2/myeloid.data.refine.Rds')
