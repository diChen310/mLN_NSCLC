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
cols.default = c('#e69f00','#54b6e9','#009e73','#f0e442','#0072b2','#d55e00',
                 '#cc79a7','#666666','#ad7700','#1c91d4','#007756','#d5c711',
                 '#005685', '#a04700','#b14380','#4d4d4d','#ffbe2d','#80c7ef',
                 '#00f6b3','#f4eb71',"black",'grey',"blue",'cyan')
cds.all = readRDS(file = './variables_v2/inter_cds.Rds')

#lym.data = cds.all[,cds.all[['majorCellType']] %in% c('T')]

lym.data = cds.all[,cds.all[['majorCellType']] %in% c('T','B','NK','Plasma')]


rm(cds.all)

gc()


####2. Re-clustering####
#selected_genes = rownames(lym.data)[ Matrix::rowSums(lym.data@assays@data$counts > 0)> 100]
lym.data <- monocle3::preprocess_cds(lym.data, num_dim = 50,method="PCA",norm_method="log") 
monocle3::plot_pc_variance_explained(lym.data)

lym.data <- monocle3::align_cds(lym.data, num_dim = 20, alignment_group = "PatientID")
lym.data <- monocle3::reduce_dimension(lym.data,cores = 20)
lym.data = monocle3::cluster_cells(lym.data, resolution=5e-5)

#1e-5, cluster 33
#5e-5, cluster 73
cols.cl = colorRampPalette(cols.default)(length(unique(monocle3::clusters(lym.data)))) 

monocle3::plot_cells(lym.data,group_label_size = 4,cell_size = 0.1,cell_stroke = 0)+ggplot2::scale_color_manual(values = cols.cl)

lym.data[['cluster.lymphocyte']]= paste0('lymphocyte_C',monocle3::clusters(lym.data))

monocle3::plot_cells(lym.data,color_cells_by = 'From',group_label_size = 4,cell_size = 0.1,cell_stroke = 0)


#####Find markers####

marker_test_res = monocle3::top_markers(lym.data,
                                        group_cells_by="cluster.lymphocyte", genes_to_test_per_group=100,
                                        reference_cells=1000, cores=32)
write.csv(marker_test_res,file = './tables_v2/lymphocyte assigned cell type marker_v2.csv')
marker_test_res.t = read.csv(file = './tables_v2/lymphocyte assigned cell type marker_v2.csv',row.names = 1)
marker_test_res$cell_group = paste(marker_test_res$cell_group,'M')

marker_test_res.2 = top_markers(lym.data[substr(rownames(lym.data),1,2) %in% c('MT','RP')==F & rownames(lym.data) != 'MALAT1',],
                              group_cells_by="cluster.lymphocyte", genes_to_test_per_group=100,
                              reference_cells=1000, cores=32)



table(lym.data[['From']])

table(lym.data[['From']],lym.data[['cluster.lymphocyte']])

plot_cells(lym.data,genes =c('CD8A','CD8B','CD4','IL7R'),scale_to_range = F)
marker.genes = c('CD3E','CD3D','CD4','CD8A','CD8B',
                 'FCGR3A','XCL1','KLRF1',
                 'S100A4','GPR183',
                 'IL7R','SELL','TCF7','CCR7','LEF1',
                 'LAG3','TIGIT','PDCD1','HAVCR2','CTLA4',
                 'IL2','GZMA','GNLY','PRF1','GZMB','GZMK','IFNG','NKG7',
                 'IL2RA','FOXP3','TGFB1',
                 'CXCL13','CXCR5','PDCD1',
                 'CREM','NR4A2',
                 'CXCR4','KLRB1',
                 'MKI67',
                 'NCAM1')
marker.genes = intersect(marker.genes,rownames(lym.data))

plot_genes_by_group(lym.data,marker.genes,
                    group_cells_by="cluster.lymphocyte",
                    ordering_type="cluster_row_col",
                    max.size=4)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))

#####Annotations####

colData(lym.data)$assigned_cell_type = lym.data[['cluster.lymphocyte']] # res 1e-5

colData(lym.data)$assigned_cell_type <- dplyr::recode(colData(lym.data)$assigned_cell_type,
                                                      
                                                      "lymphocyte_C1"="CD8 Exhausted T_C1",#LAG3, TIGIT, PDCD1, HAVCR2
                                                      'lymphocyte_C2'='CD8 Cytotoxic T_C1',
                                                      'lymphocyte_C3'='Proliferation lymphocyte_C1',
                                                       'lymphocyte_C4'='Naive B_C1',#CCR7,  SELL, 
                                                      'lymphocyte_C5'='CD4 ALDOA+ T_C1',
                                                      "lymphocyte_C6"="CD4 ALDOA+ T_C2",
                                                      'lymphocyte_C7'='CD4 Naive T_C1',#
                                                      'lymphocyte_C8'='Treg_C1',
                                                      'lymphocyte_C9'='CD4 Naive T_C2',
                                                      'lymphocyte_C10'='Plasma_C1',#
                                                      'lymphocyte_C11'='CD4 resident effector memory_C1',#£¿£¿LTB IL7R
                                                      'lymphocyte_C12'='NK_C1',
                                                      'lymphocyte_C13'='NK_C2',
                                                      'lymphocyte_C14'='CD4 resident effector memory_C2',#??????IL7R CD40LG IL32 CD69
                                                        'lymphocyte_C15'='CD4 CXCL13 T_C1',
                                                      'lymphocyte_C16'='Treg_C2',
                                                      'lymphocyte_C17'='B_C1',#
                                                      'lymphocyte_C18'='CD8 Exhausted T_C2',#LAG3
                                                      "lymphocyte_C19"="CD4 effector memory T_C1",#CCL5,S100A11, VIM, ANXA1, CRIP1
                                                      "lymphocyte_C20"="NK_C3",
                                                      'lymphocyte_C21'='Memory B_C1',#IgG,FCER2
                                                      'lymphocyte_C22'='B_C2',#
                                                      "lymphocyte_C23"="CD8 Cytotoxic T_C3",#
                                                      "lymphocyte_C24"="CD4 effector memory T_C2",#S100A11, VIM, ANXA1, CRIP1
                                                      "lymphocyte_C25"="CD4 effector memory T_C3",#S100A11, VIM, ANXA1, CRIP1
                                                      'lymphocyte_C26'='NCAM1+ NK_C1',#
                                                      'lymphocyte_C27'='B_C3',#
                                                      'lymphocyte_C28'='Treg_C3',
                                                      'lymphocyte_C29'='Remove_lowRNA',
                                                      'lymphocyte_C30'='Plasma_C2',
                                                      'lymphocyte_C31'='CD4 CXCL13 T_C2',
                                                      'lymphocyte_C32'='Remove_twoCellTypes',#Epi
                                                      'lymphocyte_C33'='CD8 Cytotoxic T_C4',
                                                      "lymphocyte_C34"="CD4 resident effector memory_C3",#CD69,LTB
                                                      'lymphocyte_C35'='B_C4',
                                                      'lymphocyte_C36'='Memory B_C2',#IgG,FCER2
                                                      'lymphocyte_C37'='Plasma_C3',#SSR4,FKBP11,DERL3
                                                      'lymphocyte_C38'='CD8 Exhausted T_C2',#LAG3
                                                      "lymphocyte_C39"="NK_C5",#?
                                                        'lymphocyte_C40'='Memory B_C3',#IgG,FCER2
                                                      'lymphocyte_C41'='NCAM1+ NK_C2',#?CD8?
                                                      'lymphocyte_C42'='CD8 Cytotoxic T_C6',
                                                      'lymphocyte_C43'='GC B_C1',#RGS13, AICDA
                                                      'lymphocyte_C44'='Remove_twoCellTypes',#CD68?
                                                      'lymphocyte_C45'='CD4 Naive T_C3',
                                                      'lymphocyte_C46'='B_C5',
                                                      'lymphocyte_C47'='B_C6',
                                                      'lymphocyte_C48'='Plasma_C4',
                                                      'lymphocyte_C49'='CD4 CXCL13 T_C3',
                                                      'lymphocyte_C50'='Plasma_C5',
                                                      'lymphocyte_C51'='B_C7',
                                                      'lymphocyte_C52'='Remove_twoCellTypes',
                                                      'lymphocyte_C53'='Anergic T',#CBLB
                                                      'lymphocyte_C54'='Remove_highRNA',
                                                      'lymphocyte_C55'='CD4 resident effector memory_C4',
                                                      'lymphocyte_C56'='Remove_twoCellTypes',
                                                      'lymphocyte_C57'='Remove_twoCellTypes',
                                                      'lymphocyte_C58'='Remove_lowRNA',
                                                      'lymphocyte_C59'='Remove_twoCellTypes',
                                                      'lymphocyte_C60'='Remove_twoCellTypes',
                                                      'lymphocyte_C61'='NK_C6',
                                                      'lymphocyte_C62'='Early pro B_C1',#PAX5, EBF1
                                                      'lymphocyte_C63'='Other T_C1',
                                                      'lymphocyte_C64'='Proliferation lymphocyte_C2',
                                                      'lymphocyte_C65'='Remove_twoCellTypes',
                                                      'lymphocyte_C66'='Remove_twoCellTypes',
                                                      'lymphocyte_C67'='Remove_twoCellTypes',
                                                      'lymphocyte_C68'='Remove_twoCellTypes',
                                                      'lymphocyte_C69'='Remove_twoCellTypes',
                                                      'lymphocyte_C70'='CD4 ALDOA+ T_C1',
                                                      'lymphocyte_C71'='IgG low plasma_C1',
                                                      'lymphocyte_C72'='Other T_C2',
                                                      'lymphocyte_C73'='Remove_lowRNA'
                                                      
)
lym.data[['subCellType']] = unlist(lapply(as.character(colData(lym.data)$assigned_cell_type), function(a){
  if(grepl('_',a)){return( strsplit(a,'_')[[1]][1])}else{return(a)}
}))

cols.cell = cols.default[1:length(unique(colData(lym.data)$subCellType))]
#cols.cell[6]='#fbc543'
names(cols.cell)=unique(colData(lym.data)$subCellType)[order(unique(colData(lym.data)$subCellType))]
# cols.clusters = generate_subColors(cols.cell ,unique(lym.data[['assigned_cell_type']]))
# names(cols.clusters)= as.character(unique(lym.data[['assigned_cell_type']]))




plot_cells(lym.data,group_cells_by = 'subCellType',color_cells_by = 'subCellType',
           cell_stroke = 0,cell_size = 0.5,label_cell_groups = F,
           group_label_size = 4)+scale_color_manual(values = cols.cell )+
  theme_void()

saveRDS(lym.data,file = './variables_v2/lym.data.Rds')

marker.t = read.csv(file = './tables_v2/lymphocyte assigned cell type marker_v2.csv',row.names = 1)

####Re-analysis after remove####
marker_test_res.all = read.csv(file = './tables_v2/markers of all sub-clusters.csv')
marker_test_res.all$cell_group = paste(marker_test_res.all$cell_group,'M')

lym.data = readRDS(file = './variables_v2/lym.data.Rds')

lym.data = lym.data[,lym.data[['subCellType']] != 'Remove']

cc=table(lym.data[['From']],lym.data[['cluster.lymphocyte']])/apply(table(lym.data[['From']],lym.data[['cluster.lymphocyte']]),1,sum)
heatmap(cc)

table(lym.data[['subCellType']],lym.data[['majorCellType']])

B.data = lym.data[,lym.data[['subCellType']] %in% c('B','Early pro B','GC B','IgG low plasma','Memory B','Naive B','Plasma') ]
cc=table(B.data[['From']],B.data[['cluster.lymphocyte']])/apply(table(B.data[['From']],B.data[['cluster.lymphocyte']]),1,sum)
heatmap(cc)
colData(B.data)$assigned_cell_type = B.data[['cluster.lymphocyte']] # 

colData(B.data)$assigned_cell_type <- dplyr::recode(colData(B.data)$assigned_cell_type,
                                                      
                                                      'lymphocyte_C4'='Naive B_C1',#CCR7,  SELL, 
                                                      'lymphocyte_C21'='Naive B_C2',#IgG,FCER2
                                                      'lymphocyte_C51'='Activated B_C7',#NR4A1
                                                      'lymphocyte_C43'='GC B_C1',#RGS13, AICDA
                                                      'lymphocyte_C40'='Naive B_C3',#IGHD, SELL
                                                      'lymphocyte_C27'='Naive B_C4',## SELL
                                                      'lymphocyte_C30'='Plasma_C2',
                                                      'lymphocyte_C10'='Plasma_C1',#
                                                      'lymphocyte_C37'='Plasma_C3',#SSR4,FKBP11,DERL3
                                                      'lymphocyte_C50'='Plasma_C5',
                                                      'lymphocyte_C48'='Plasma_C4',
                                                      'lymphocyte_C17'='B_C1',#
                                                      'lymphocyte_C71'='IgG low plasma_C1',
                                                      'lymphocyte_C62'='ZBTB20+ B_C1',#PAX5, EBF1
                                                      'lymphocyte_C22'='B_C2',#
                                                      'lymphocyte_C47'='B_C3',
                                                      'lymphocyte_C35'='B_C4',
                                                      'lymphocyte_C46'='B_C5',
                                                      'lymphocyte_C36'='Naive B_C5'
                                                    
)
B.data[['subCellType']] = unlist(lapply(as.character(colData(B.data)$assigned_cell_type), function(a){
  if(grepl('_',a)){return( strsplit(a,'_')[[1]][1])}else{return(a)}
}))


cols.cell = cols.default[1:length(unique(colData(B.data)$subCellType))]
names(cols.cell)=levels(colData(B.data)$subCellType)
scater::plotReducedDim(B.data,dimred = 'UMAP',colour_by = 'subCellType',point_size = 1,text_by = 'subCellType',text_size = 2)+
  ggplot2::scale_color_manual(values = cols.cell)+
  theme_void()
saveRDS(B.data,file = './variables/B.data.Rds')
scater::plotReducedDim(B.data,dimred = 'UMAP',colour_by = 'cluster.lymphocyte',point_size = 1,text_by = 'subCellType',text_size = 2)+
  theme_void()

T.data = lym.data[,lym.data[['subCellType']] %in% c('B','Early pro B','GC B','IgG low plasma','Memory B','Naive B','Plasma') == F ]
cc=table(T.data[['From']],T.data[['cluster.lymphocyte']])/apply(table(T.data[['From']],T.data[['cluster.lymphocyte']]),1,sum)
heatmap(cc)

colData(T.data)$assigned_cell_type = T.data[['cluster.lymphocyte']] #

colData(T.data)$assigned_cell_type <- dplyr::recode(colData(T.data)$assigned_cell_type,
                                                    'lymphocyte_C5'='ALDOA+ T_C1',
                                                    "lymphocyte_C6"="ALDOA+ T_C2",
                                                    'lymphocyte_C38'='CD8 Exhausted T_C1',#LAG3 IFNG
                                                    'lymphocyte_C13'='NK_C2',
                                                    "lymphocyte_C20"="NK_C3",
                                                    'lymphocyte_C12'='NK_C1',
                                                     "lymphocyte_C19"="CD4 Effector Memory_C1",#S100A4,IL32,ANXA1,CCL5
                                                    'lymphocyte_C61'='NK_C6',
                                                     "lymphocyte_C25"="CD4 Effector Memory_C2",#RORA,IFNG,CXCR6
                                                   
                                                     'lymphocyte_C70'='ALDOA+ T_C3',
                                                    'lymphocyte_C72'='Exhausted T_C1',#TIGIT,HAVCR2
                                                    'lymphocyte_C63'='¦Ã¦Ä T_C1',#TRDC
                                                    'lymphocyte_C64'='Proliferation lymphocyte_C2',
                                                    
                                                    "lymphocyte_C24"="CD4 Effector Memory_C3",#CD69,IFNG,RORA
                                                    "lymphocyte_C39"="NK_C5",#
                                                    'lymphocyte_C41'='NK_C4',#
                                                    "lymphocyte_C23"="CD8 Cytotoxic T_C1",#
                                                    
                                                    'lymphocyte_C45'='CD4 Naive T_C3',
                                                    'lymphocyte_C11'='CD4 Naive T_C4',#
                                                    
                                                    'lymphocyte_C26'='NK_C7',#
                                                    'lymphocyte_C18'='CD8 Exhausted T_C2',#LAG3,RGS1,IFNG
                                                    'lymphocyte_C2'='CD8 Cytotoxic T_C2',
                                                    
                                                    "lymphocyte_C1"="CD8 Exhausted T_C3",#LAG3, TIGIT, PDCD1, HAVCR2
                                                    'lymphocyte_C15'='CXCL13+ T_C1',#CXCL13, CTLA4, PDCD1, TOX2,TIGIT
                                                    'lymphocyte_C3'='Proliferation lymphocyte_C1',
                                                    
                                                    'lymphocyte_C55'='CD4 Resident Effector Memory_C3',#RORA,IL7R,CD69
                                                    'lymphocyte_C42'='CD8 Cytotoxic T_C3',
                                                    'lymphocyte_C33'='CD8 Cytotoxic T_C4',
                                                    'lymphocyte_C14'='CD4 Resident Effector Memory_C1',#CD69,CD40LG,IL7R,S100A4
                                                    
                                                    'lymphocyte_C7'='CD4 Naive T_C1',#
                                                    'lymphocyte_C9'='CD4 Naive T_C2',
                                                    'lymphocyte_C53'='ARHGAP15+ T',#
                                                    'lymphocyte_C31'='CXCL13+ T_C2',
                                                    'lymphocyte_C8'='Treg_C1',
                                                    'lymphocyte_C16'='Treg_C2',
                                                    'lymphocyte_C28'='Treg_C3',
                                                     "lymphocyte_C34"="CD4 Resident Effector Memory_C2",#CD69,CD40LG,IL7R,S100A4
                                                    'lymphocyte_C49'='CXCL13+ T_C3'
                                                     
                                                 
)

marker.genes = c('CD3E','CD3D','CD4','CD8A','CD8B',
                 'FCGR3A','XCL1','KLRF1',
                 'S100A4','GPR183',
                 'IL7R','SELL','TCF7','CCR7','LEF1',
                 'RORA','CD69',
                 'LAG3','TIGIT','PDCD1','HAVCR2','CTLA4',
                 'IL2','GZMA','GNLY','PRF1','GZMB','GZMK','IFNG','NKG7',
                 'IL2RA','FOXP3','TGFB1',
                 'CXCL13','CXCR5','PDCD1',
                 'MKI67',
                 'NCAM1','ARHGAP15','ALDOA',
                 'TRDC')

plot_genes_by_group(T.data,marker.genes,
                    group_cells_by="assigned_cell_type",
                    ordering_type="cluster_row_col",
                    max.size=4)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))


scater::plotReducedDim(T.data,dimred = 'UMAP',colour_by = 'assigned_cell_type',point_size = 1,text_by = 'assigned_cell_type',text_size = 2)

T.data[['subCellType']] = unlist(lapply(as.character(colData(T.data)$assigned_cell_type), function(a){
  if(grepl('_',a)){return( strsplit(a,'_')[[1]][1])}else{return(a)}
}))
cols.cell = cols.default[1:length(unique(colData(T.data)$subCellType))]
names(cols.cell)=levels(colData(T.data)$subCellType)
scater::plotReducedDim(T.data,dimred = 'UMAP',colour_by = 'subCellType',point_size = 1,text_by = 'subCellType',text_size = 2)+
  ggplot2::scale_color_manual(values = cols.cell)+
  theme_void()
T.data = T.data[,T.data[['subCellType']] != 'Proliferation lymphocyte']
saveRDS(T.data,file = './variables_v2/T.data.Rds')
