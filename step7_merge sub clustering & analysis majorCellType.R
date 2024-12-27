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

####7.1 merge sub clustering results####

cds.all = readRDS(file = './variables_v2/inter_cds.Rds')

Myeloid.data = readRDS(file = './variables_v2/myeloid.data.Rds')
lym.data = readRDS(file = './variables_v2/lym.data.Rds')
stromal.data <- readRDS(file = './variables_v2/stromal.data.Rds')
epi.data = readRDS(file = './variables_v2/epi.data.Rds')


cds.all[['subClusters']] = cds.all[['majorCellType']] 

meta.data = as.data.frame(colData(cds.all))
meta.data[colnames(epi.data),'subClusters'] = epi.data[['cluster.epi']]
meta.data[colnames(Myeloid.data),'subClusters'] = Myeloid.data[['cluster.Myeloid']]
meta.data[colnames(stromal.data),'subClusters'] = stromal.data[['cluster.stromal']]
meta.data[colnames(lym.data),'subClusters'] = lym.data[['cluster.lymphocyte']]

cds.all[['subClusters']] = meta.data$subClusters


table(cds.all[['subClusters']])

cds.all = cds.all[,cds.all[['subClusters']] != 'Epi'] #only one Epi in the lymph normal
marker_test_res <- top_markers(cds.all, group_cells_by="subClusters", genes_to_test_per_group=100,
                               reference_cells=1000, cores=16)
write.csv(marker_test_res,file = './tables_v2/markers of all sub-clusters.csv')
marker_test_res$cell_group = paste(marker_test_res$cell_group,'M')
write.csv(as.data.frame(colData(cds.all)),file = './tables_v2/meta data all_temp.csv')

meta.data.show = cds.all@colData
meta.data.show = as.data.frame(meta.data.show)

# 
ggplot2::ggplot(meta.data.show,ggplot2::aes(x=subClusters,y=nFeature_RNA,fill=subClusters))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')


colData(cds.all)$assigned_cell_type_v1 <- dplyr::recode(colData(cds.all)$subClusters,
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
                                                 "Epi_C18"="Remove_lowRNA",
                                                 
                                                 "Myeloid_C1"="MP_C1",
                                                 'Myeloid_C2'='FABP4+ MP_C1',
                                                 'Myeloid_C3'='STAB1+SPP1+ MP_C1',
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
                                                 "Myeloid_C24"="FABP4+ MP_C2",
                                                 
                                                 "Stromal_C1"="mCAF_C1",
                                                 "Stromal_C2"="Veins",#
                                                 "Stromal_C3"="iCAF",
                                                 "Stromal_C4"="vCAF",#NOTCH3, COL18A1
                                                 "Stromal_C5"="Arteries",#HEY1,IGFBP3
                                                 "Stromal_C6"="Remove_twoCellTypes",#C7
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
                                                 'Stromal_C19'='Remove_lowNumber',
                                                 
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
                                                 'lymphocyte_C11'='CD4 resident effector memory_C1',#LTB IL7R
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

cds.all[['subCellType']] = unlist(lapply(as.character(colData(cds.all)$assigned_cell_type_v1), function(a){
  if(grepl('_',a)){return( strsplit(a,'_')[[1]][1])}else{return(a)}
}))




write.csv(as.data.frame(colData(cds.all)),file = './tables_v2/meta data all_subCellTypes.csv')
#cds.all = cds.all[,cds.all[['subClusters']] != 'Epi'] #only one Epi in the lymph normal

#meta.data.all = read.csv(file = './tables_v2/meta data all_subCellTypes.csv')

#### update based on sub clustering####
meta.data.all = as.data.frame(colData(cds.all))
meta.data.all[meta.data.all$subCluster =='Stromal_C11','majorCellType' ] = 'Mast'
cds.all[['majorCellType']] = meta.data.all$majorCellType
cds.all[['subCellType']] = meta.data.all$subCellType
cds.all[['assigned_cell_type_v1']] = meta.data.all$assigned_cell_type_v1


cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559","#32709f",
                 '#b8a3c4','#60428d','#9f5633','grey','#c92e2a','#c79795','#d9c17c',"#d0896b",'#5d94b4',"black","blue",'cyan')

cds = cds.all[,cds.all[['subCellType']] != 'Remove' & cds.all[['majorCellType']] != 'Remove']
cols.cl =cols.default[1:length(unique(colData(cds)$majorCellType))]
names(cols.cl)=c('T','B','Plasma','Epi','Mast','Myeloid','Fibro','Endo','NK')

scater::plotReducedDim(cds,dimred = 'UMAP',colour_by = 'majorCellType',point_size = 1)+
  ggplot2::scale_color_manual(values = cols.cl)+
  theme_void()+theme(legend.position = 'none')

scater::plotReducedDim(cds[,sample(colnames(cds),2000)],dimred = 'UMAP',colour_by = 'majorCellType',point_size = 1)+
  ggplot2::scale_color_manual(values = cols.cl)+
  theme_void()

cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')
plot_cells(cds,  color_cells_by="From",
           cell_size=0.1,
           cell_stroke=0,
           rasterize = T,
           label_cell_groups = F)+
  ggplot2::scale_color_manual(values = cols.From)+
  facet_wrap(~From, ncol = 3)+
  theme_void()+
  theme( legend.position = 'none')
cols.sources = jcolors('pal8')
names(cols.sources)=unique(cds[['Source']])
plot_cells(cds, color_cells_by="Source", label_cell_groups=FALSE,cell_size = 0.1,cell_stroke = 0)+
  facet_wrap(~Source,ncol = 4)+
  scale_color_manual(values = cols.sources)+
  theme_void()+theme(legend.position = 'none')
####7.3 statistics####
meta.data.show = cds@colData
meta.data.show = as.data.frame(meta.data.show)
meta.data.show$Disease2 = substr(meta.data.show$Disease,1,4)


#write.csv(meta.data.show,file = './tables_v2/all cds meta data.csv')


ggplot2::ggplot(meta.data.show,ggplot2::aes(x=From,fill=majorCellType))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values = cols.cl)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))




meta.data.show$Patient = factor(meta.data.show$Patient,
                                levels = paste0('P',0:18))
ggplot2::ggplot(meta.data.show[!is.na(meta.data.show$Patient),],ggplot2::aes(x=orig.ident,fill=majorCellType))+
  facet_grid(~Patient,scales = 'free',space = 'free')+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values = cols.cl)+
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))


marker.genes.classical = c('EPCAM','KRT17','KRT19',
                           'CD79A','MS4A1',"CD19",
                           'JCHAIN','MZB1','IGKC','IGHM',
                           'CD3D','CD3E','TRBC2',
                           'APOE','MARCO','CD68','LYZ',
                           'COL1A1','COL1A2','ACTA2',
                           'CLDN5','PECAM1','VWF',
                           'TPSAB1','TPSB2',
                           'NKG7','KLRB1','KLRD1'
)
plot_genes_by_group(cds,
                    marker.genes.classical,
                    group_cells_by="majorCellType",
                    ordering_type="cluster_row_col",
                    max.size=3)+scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))

saveRDS(cds,file = './variables_v2/inter_cds_refine.Rds')
cds = readRDS(file = './variables_v2/inter_cds_refine.Rds')
