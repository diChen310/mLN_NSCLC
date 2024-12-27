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
library(ggsci)
####11.1 Read previous results####
cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559","#32709f",
                 '#b8a3c4','#60428d','#9f5633','grey','#c92e2a','#c79795','#d9c17c',"#d0896b",'#5d94b4',"black","blue",'cyan')

stromal.data = readRDS(file = './variables_v2/stromal.data.refine.Rds')
marker.stromal = read.csv(file = './tables_v2/Stromal cluster marker.r2.csv',row.names = 1)
marker.genes = c('COL1A1','COL3A1','ACTA2',
                 'MMP11','COL1A2','MYH11',#mCAF
                 'C3','CFD','CXCL12','CXCL14','DPP4',#iCAF
                 'NOTCH3','COL18A1',#vCAF
                 'MCAM' ,'RGS5',#pericyte
                 'GAPDH','ENO1','VEGFA','MME',#tCAF
                 'LYVE1','CCL21',#LEC
                 'PECAM1','VWF','CLDN5',#BEC,
                 'CA4','HPGD','ACKR1','HEY1','C7',
                 'SOSTDC1','EDNRB',
                 'CST1','WIF1',
                 'HEY1','COX4I2','PROX1',
                 'MKI67','TOP2A','CD34','KDR'
)


plot_genes_by_group(stromal.data,marker.genes,
                    group_cells_by="subCellType",
                    ordering_type="cluster_row_col",
                    max.size=4)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))+
  theme(axis.text.x = element_text(angle = 90))

plot_genes_violin(stromal.data[c('CST1','CST3'),] ,
                  group_cells_by = 'cluster.stromal')+ scale_fill_manual(values = cols.cl)


cols.cell = cols.default[1:length(unique(colData(stromal.data)$subCellType))]
names(cols.cell)=unique(colData(stromal.data)$subCellType)
scater::plotReducedDim(stromal.data,dimred = 'UMAP',colour_by = 'subCellType',point_size = 1)+
ggplot2::scale_color_manual(values = cols.cell)+
theme_void()

####11.2 statistics on the cell sub types####

#####subCellType barplot in LUAD and LUSC####

meta.data.show = stromal.data@colData
meta.data.show = as.data.frame(meta.data.show)
meta.data.show$Disease2 = substr(meta.data.show$Disease,1,4)
meta.data.show$subCellType = as.character(meta.data.show$subCellType)


ggplot2::ggplot(meta.data.show[meta.data.show$subCellType %in% c('LEC','Capillaries','Veins','Arteries','Aerocyte','Endo progenitor'),],ggplot2::aes(x=From,fill=subCellType))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+facet_wrap(~Disease2,ncol = 1)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values =  cols.cell)+
  theme(axis.text  = element_text(size = 8),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))


ggplot2::ggplot(meta.data.show[meta.data.show$subCellType %in% c('LEC','Capillaries','Veins','Arteries','Aerocyte','Endo progenitor','Proliferation stromal') == F,],ggplot2::aes(x=From,fill=subCellType))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+facet_wrap(~Disease2,ncol = 1)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values =  cols.cell)+
  theme(axis.text  = element_text(size = 8),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))

##### clusters difference#####
meta.data.show.endo = meta.data.show[meta.data.show$subCellType %in% c('LEC','Capillaries','Veins','Arteries','Aerocyte','Endo progenitor') == F,]
com2orig = table(meta.data.show.endo$orig.ident,meta.data.show.endo$subCellType)

f.data<-as.matrix(com2orig)
f.data = f.data/apply(f.data, 1, sum)
#f.data = t(f.data)
#rownames(ids.all)=ids.all$id
meta.data.u = meta.data.show
meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
rownames(meta.data.u)=meta.data.u$orig.ident

com2orig.df = melt(f.data,measure.vars = colnames(f.data))
com2orig.df$patient = meta.data.u[as.character(com2orig.df$Var1),'PatientID']
com2orig.df$From = meta.data.u[as.character(com2orig.df$Var1),'From']
cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')

ggplot(com2orig.df,aes(x=From,y=value,fill=From,color=From))+geom_violin(draw_quantiles = 0.5,scale = 'width')+
  facet_wrap(~Var2,scales = 'free',ncol = 5)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 0))+
  scale_fill_manual(values = alpha(cols.From,0.7))+
  scale_color_manual(values = cols.From)



#####Subtype markers#####

marker_test_res.subCellType = top_markers(stromal.data, group_cells_by="subCellType", genes_to_test_per_group=200,
                                          reference_cells=1000, cores=32)



filter.vCAF = filter(marker_test_res.subCellType,marker_test_p_value < 0.01 & 
                       substr(gene_id,1,2) %in% c('MT','RP') == F & 
                       cell_group == 'vCAF' &
                       gene_id != 'MALAT1') %>% 
  arrange(cell_group,marker_test_p_value)
filter.vCAF$label = ''
filter.vCAF$label[1:25]=filter.vCAF$gene_id[1:25]
ggplot(filter.vCAF,aes(color=log2(mean_expression+1),y=-log10(marker_test_p_value),x=specificity))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label = label),size=2,max.overlaps = 1000)+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours = hcl.colors(100))
# plot_genes_by_group(stromal.data,as.character(filter.from$gene_id),
#                     group_cells_by="From",
#                     ordering_type="maximal_on_diag",
#                     max.size=4)+
#   scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))
# plot_cells(stromal.data,genes = c('IGFBP7','COL18A1','NOTCH3','RGS5','COL4A1'),cell_size = 0.1,cell_stroke = 0,scale_to_range = F,label_cell_groups = F)+
#   scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))
filter.CST1 = filter(marker_test_res.subCellType,marker_test_p_value < 0.01 & 
                       substr(gene_id,1,2) %in% c('MT','RP') == F & 
                       cell_group == 'CST1+ Fibro' &
                       gene_id != 'MALAT1') %>% 
  arrange(cell_group,-specificity)
filter.CST1$label = ''
filter.CST1$label[1:25]=filter.CST1$gene_id[1:25]
filter.CST1$label[filter.CST1$gene_id %in% c('MMP2','ZEB1','IGFBP3','CD302')]
ggplot(filter.CST1,aes(color=log2(mean_expression+1),y=-log10(marker_test_p_value),x=specificity))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label = label),size=2,max.overlaps = 1000)+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours = hcl.colors(100))

write.csv(marker_test_res.subCellType,file = './tables_v2/Fibro subtype markers.csv')

####11.3 cell type based From difference#####
cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')

marker.diff.from.endo = top_markers(stromal.data[,stromal.data[['subCellType']] %in% c('LEC','Capillaries','Veins','Arteries','Aerocyte','Endo progenitor')], group_cells_by="From", genes_to_test_per_group=300,
                               reference_cells=1000, cores=24)
marker.diff.from.fibro = top_markers(stromal.data[,stromal.data[['subCellType']] %in% c('LEC','Capillaries','Veins','Arteries','Aerocyte','Endo progenitor') == F], group_cells_by="From", genes_to_test_per_group=300,
                                    reference_cells=1000, cores=24)

endo.data = stromal.data[,stromal.data[['subCellType']] %in% c('LEC','Capillaries','Veins','Arteries','Aerocyte','Endo progenitor')]
endo.sc = CreateSeuratObject(counts = assay(endo.data))
endo.sc = AddMetaData(endo.sc,as.data.frame(colData(endo.data)))
Idents(endo.sc)=endo.sc$From
fibro.data = stromal.data[,stromal.data[['subCellType']] %in% c('LEC','Capillaries','Veins','Arteries','Aerocyte','Endo progenitor','Proliferation stromal') == F]
fibro.sc = CreateSeuratObject(counts = assay(fibro.data))
fibro.sc = AddMetaData(fibro.sc,as.data.frame(colData(fibro.data)))
Idents(fibro.sc)=fibro.sc$From

endo.mLN.diff = FindMarkers(endo.sc,ident.1 = 'mLN',only.pos = F)
endo.mLNPT.diff = FindMarkers(endo.sc,ident.1 = 'mLN+ PT',only.pos = F)
fibro.mLN.diff = FindMarkers(fibro.sc,ident.1 = 'mLN',only.pos = F)
fibro.mLNPT.diff = FindMarkers(fibro.sc,ident.1 = 'mLN+ PT',ident.2 = 'mLN- PT',only.pos = F)

saveRDS(endo.mLN.diff,file = './variables_v2/endo.mLN.diff.Rds')
saveRDS(endo.mLNPT.diff,file = './variables_v2/endo.mLNPT.diff.Rds')
saveRDS(fibro.mLN.diff,file = './variables_v2/fibro.mLN.diff.Rds')
saveRDS(fibro.mLNPT.diff,file = './variables_v2/fibro.mLNPT.diff.Rds')
marker.diff.from.endo.top = filter(marker.diff.from.endo,mean_expression > 1 & marker_test_p_value<0.01) %>% group_by(cell_group) %>% top_n(10,specificity)
marker.diff.from.fibro.top = filter(marker.diff.from.fibro,mean_expression > 1 & marker_test_p_value<0.01)%>% group_by(cell_group) %>% top_n(10,specificity)

plot_genes_violin(fibro.data[c('TXNRD1','RGS5','PLEKHA7','AMBRA1',
                                 'IGFBP7','COL18A1','TGFB1','COL8A1'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot_genes_violin(fibro.data[c('CST1','CCL19','TGFB1','TGFBI'),fibro.data[['From']]],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))
plot_genes_violin(fibro.data[c('CST1','SMOC2','ISLR','WISP2'),fibro.data[['From']]],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))


plot_genes_violin(fibro.data[c('CST1'),fibro.data[['From']] %in% c('mLN+ PT','mLN','mLN- PT')],
                  group_cells_by = 'From',normalize = T)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))


# plot_genes_violin(fibro.data[c('CCL19','CXCL9'),],
#                   group_cells_by = 'From',
#                   ncol = 4)+
#   scale_fill_manual(values = cols.From)+
#   theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot_genes_violin(endo.data[c('COL4A1','SPARC','AMBRA1','HSPA1A',
                              'PLEKHA7','TXNRD1','PLVAP','MGP'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot_genes_violin(stromal.data[c('CST1'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))



endo.mLN.markers = intersect(rownames(filter(endo.mLN.diff,p_val_adj<0.01 & avg_log2FC > 0.5)),
                             rownames(filter(endo.mLNPT.diff,p_val_adj<0.01 & avg_log2FC > 0.5)))
fibro.mLN.markers = intersect(rownames(filter(fibro.mLN.diff,p_val_adj<0.01 & avg_log2FC > 0.5)),
                             rownames(filter(fibro.mLNPT.diff,p_val_adj<0.01 & avg_log2FC > 0.5)))

intersect(endo.mLN.markers,fibro.mLN.markers)
#####pathway enrichment#####

#####go enrichment#####
library(clusterProfiler)
go.items = read.gmt(file('./c5.all.v2024.1.Hs.symbols.gmt'))

go.endo.mLN = enricher(endo.mLN.markers,TERM2GENE = go.items)@result

go.fibro.mLN = enricher(fibro.mLN.markers,TERM2GENE = go.items)@result

hallmark.items = read.gmt(file('./MSigDB_Hallmark.gmt'))

h.endo.mLN = enricher(endo.mLN.markers,TERM2GENE = hallmark.items)@result
h.fibro.mLN = enricher(fibro.mLN.markers,TERM2GENE = hallmark.items)@result

h.endo.mLN$CellType = 'Endo'
h.fibro.mLN$CellType = 'Fibroblast'

merge.hallmarker = rbind(h.endo.mLN,h.fibro.mLN)

merge.hallmarker = filter(merge.hallmarker,p.adjust<0.05) %>% group_by(CellType) %>% arrange(p.adjust)

merge.hallmarker$ID = factor(merge.hallmarker$ID,levels = rev(merge.hallmarker$ID))
merge.hallmarker$CellType = factor(merge.hallmarker$CellType,levels = c('Fibroblast','Endo'))
ggplot(merge.hallmarker,aes(x=-log10(p.adjust),fill = -log10(p.adjust),y=ID))+geom_bar(stat = 'identity',width = 0.6)+
  theme_cowplot(font_size = 8)+
  facet_grid(CellType~'',space='free',scales = 'free_y')+
  scale_fill_gradientn(colors = hcl.colors(100))+
  theme(strip.background = element_blank(),
        strip.text = element_text(angle = 90))
write.csv(merge.hallmarker,file = './tables_v2/stromal From diff hallmark.csv')


cst1.path = enricher(filter.CST1$gene_id,TERM2GENE = hallmark.items)@result

cst1.hallmarker = filter(cst1.path,p.adjust<0.01) %>% arrange(p.adjust)

cst1.hallmarker$ID = factor(cst1.hallmarker$ID,levels = rev(cst1.hallmarker$ID))
ggplot(cst1.hallmarker,aes(x=-log10(p.adjust),fill = -log10(p.adjust),y=ID))+geom_bar(stat = 'identity',width = 0.6)+
  theme_cowplot(font_size = 8)+
  scale_fill_gradientn(colors = hcl.colors(100))+
  theme(strip.background = element_blank(),
        strip.text = element_text(angle = 90))

write.csv(cst1.path,file = './tables_v2/stromal CST1+ fibro markers pathway.csv')


#### 11.4 EMT and Angiogenesis correlations####
hallmark.items = read.gmt(file('./MSigDB_Hallmark.gmt'))
EMT.genes = unique(hallmark.items[hallmark.items$term == 'EPITHELIAL_MESENCHYMAL_TRANSITION',2])
ang.genes = unique(hallmark.items[hallmark.items$term == 'ANGIOGENESIS',2])
cds.GSVA = readRDS(file = './variables_v2/cds.hallmark.gsva.Rds')
stromal.st = CreateSeuratObject(counts = assay(stromal.data))
stromal.st = AddMetaData(stromal.st,as.data.frame(colData(stromal.data)))
stromal.st = NormalizeData(stromal.st)
stromal.st = AddModuleScore(stromal.st,features = list(EMT = EMT.genes),name = 'EMT_Score')
stromal.st = AddModuleScore(stromal.st,features = list(ANG = ang.genes),name = 'ANGIOGENESIS_Score')



stromal.data[['EMT']] = stromal.st$EMT_Score1
stromal.data[['ANGIOGENESIS']] = stromal.st$ANGIOGENESIS_Score1

plot_cells(stromal.data,color_cells_by = 'EMT',cell_size = 0.25,cell_stroke = 0,label_cell_groups = F)

meta.data.plot = as.data.frame(colData(stromal.data))

ggplot(meta.data.plot,aes(x=subCellType,y=EMT,fill=subCellType))+geom_boxplot(outlier.size = 0.2)+
  scale_fill_manual(values = cols.cell)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggplot(meta.data.plot[meta.data.plot$subCellType %in% 
                        c('LEC','Capillaries','Veins','Arteries',
                          'Aerocyte','Endo progenitor','Proliferation stromal') == F,],
       aes(x=reorder(subCellType,EMT,FUN=median),y=EMT,fill=subCellType))+geom_boxplot(outlier.size = 0.2)+
  scale_fill_manual(values = cols.cell)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggplot(meta.data.plot[meta.data.plot$subCellType %in% 
                        c('LEC','Capillaries','Veins','Arteries',
                          'Aerocyte','Endo progenitor'),],
       aes(x=reorder(subCellType,ANGIOGENESIS,FUN=median),y=ANGIOGENESIS,fill=subCellType))+
  geom_boxplot(outlier.size = 0.2)+
  scale_fill_manual(values = cols.cell)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))



saveRDS(stromal.data,file = './variables_v2/stromal.data.refine.Rds')

####5. Spatial####
st.m = readRDS(file='./variables_v2/st.m_bayesSpace.Rds')
names(st.m@images) = unique(st.m$orig.ident)

st.f = st.m[,st.m$Source == 'PLM']

st.m.fibro = st.m[,st.m$cluster.name == 'Fibroblast']

Idents(st.m.fibro)=st.m.fibro$Source

fibro.st.diff = FindMarkers(st.m.fibro,ident.1 = 'PLM')
write.csv(fibro.st.diff,file = './tables_v2/fibro.st.diff.Source.csv')
SpatialFeaturePlot(st.m,
                   features =c('CST1'),slot = 'data',images = c('P6L','P6T'),
                   ncol = 8,
                   alpha = 0.7,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0)
SpatialFeaturePlot(st.m.fibro,
                   features =c('IGHG4'),slot = 'data',images = c('P6L','P6T'),
                   ncol = 8,
                   alpha = 0.7,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0)
SpatialDimPlot(st.m.fibro,group.by = 'cluster.name',
                   images = c('P6L','P6T'), 
                   ncol = 4,
                   alpha = 0.8,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 1,
               cols = cols.cl)
VlnPlot(st.f.fibro[,st.f.fibro$orig.ident %in% c('P6L6','P6T'),],group.by = 'orig.ident',features = c('COL4A1','IGFBP7'),pt.size = 0 )

SpatialFeaturePlot(st.f,
                   features =c('CST1'),
                   images = c('P3T','P3L','P6T','P6L'), 
                   ncol = 4,
                   alpha = 1,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0.6)
SpatialFeaturePlot(st.f,
                   features =c('CST1'),
                   images = c('P5T','P5L','P7T','P7L'), 
                   ncol = 4,
                   alpha = 1,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0.6)
SpatialFeaturePlot(st.m,
                   features =c('CST1'),
                   images = c( "P10_T1", "P11_T1","P15_T2", "P16_T1" ,
                               "P17_T2","P19_T1","P24_T1", "P25_T2"), 
                   ncol = 4,
                   alpha = 1,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0.6)


SpatialFeaturePlot(st.f.fibro,
                   features =c('CST1'),
                   images = c('P3T','P3L','P6T','P6L'), 
                   ncol = 4,
                   alpha = 1,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0)



#####CST1 analysis#####

cds.sm.fibro = new_cell_data_set(as(as.matrix(st.m.fibro@assays$Spatial@counts[ rownames(st.m.fibro),]), "sparseMatrix"),
                                 cell_metadata = st.m.fibro@meta.data,
                                 gene_metadata = data.frame(id = rownames(st.m.fibro),gene_short_name=rownames(st.m.fibro),
                                                            row.names = rownames(st.m.fibro))
)
cds.sm.fibro[['From']] = factor(cds.sm.fibro[['From']],levels = c('mLN','mLN+ PT','mLN- PT'))
plot_genes_violin(cds.sm.fibro['CST1',],group_cells_by = 'From',normalize = T)+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(2,3),c(1,3)))+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

st.m$CST1_Level = ifelse(st.m@assays$Spatial@data['CST1',] > median(st.m@assays$Spatial@data['CST1',]),'High','Low')
Idents(st.m) = st.m$CST1_Level

st.m.cst1.diff = FindMarkers(st.m,ident.1 = 'High')
write.csv(st.m.cst1.diff,file = './tables_v2/st.m.cst1.diff.csv')

VlnPlot(st.m,features = c('CST1','IGHG4'),group.by = 'CST1_Level',pt.size = 0,adjust = 5)



up.CST1.genes = rownames(dplyr::filter(st.m.cst1.diff,p_val_adj<0.01 & avg_log2FC > 0.5))
library(clusterProfiler)
go.items = read.gmt(file('./c5.all.v2024.1.Hs.symbols.gmt'))
hallmark.items = read.gmt(file('./MSigDB_Hallmark.gmt'))

up.CST1.hallmark = enricher(up.CST1.genes,TERM2GENE = hallmark.items)@result
write.csv(up.CST1.hallmark,file = './tables_v2/CST1.up.hallmark.csv')
up.CST1.hallmark.filter = dplyr::filter(up.CST1.hallmark,p.adjust < 0.05)

ggplot(up.CST1.hallmark.filter,aes(x=-log10(p.adjust),y=Count,color=-log10(p.adjust),size = -log10(p.adjust)))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label = ID,size = 3))+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours = hcl.colors(100))
#####Distance with epi#####
samples = unique(st.f$orig.ident)

unique(st.f$orig.ident)
fibro.spots = c()
fibro.distance.min = c()
fibro.distance.median = c()
orig.flags = c()

for(sample in samples){
  
  st.f.i.fibro = colnames(st.f[,st.f$orig.ident == sample & st.f$cluster.name == 'Fibroblast'])
  st.f.i.epi = colnames(st.f[,st.f$orig.ident == sample & st.f$cluster.name == 'Epithelial'])
  
  st.f.i.distance = epiDistance(st.f,sample,st.f.i.fibro,st.f.i.epi)
  
  
  
  fibro.spots = append(fibro.spots,st.f.i.fibro)
  fibro.distance.min = append(fibro.distance.min,
                              apply(st.f.i.distance,1,min))
  fibro.distance.median = append(fibro.distance.median,
                                 apply(st.f.i.distance, 1, median))
  orig.flags = append(orig.flags,rep(sample,length(st.f.i.fibro)))
  print(sample)
}

fibro.distance.res = data.frame(row.names = fibro.spots,
                                distance.min = fibro.distance.min,
                                distance.median = fibro.distance.median,
                                sample = orig.flags)
median(fibro.distance.res$distance.min)
max(fibro.distance.res$distance.min)

fibro.distance.res$distance.tag = ifelse(fibro.distance.res$distance.min> 3,'Far','Near')

st.f.fibro = st.f[,st.f$cluster.name == 'Fibroblast']

st.f.fibro$distance.tag = fibro.distance.res[colnames(st.f.fibro),'distance.tag']

st.f.fibro$distanceToEpi = fibro.distance.res[colnames(st.f.fibro),'distance.min']

Idents(st.f.fibro)=st.f.fibro$distance.tag
distance.diff = FindMarkers(st.f.fibro,ident.1 = 'Near',logfc.threshold = 0)
st.f.fibro = NormalizeData(st.f.fibro)
SpatialDimPlot(st.f.fibro,group.by = 'distance.tag',
               images = c('P3L','P5L','P6L','P7L','P3T','P5T','P6T','P7T'), 
               pt.size.factor = 1.5,
               stroke = 0,
               image.alpha = 0.5,
               ncol=4)
SpatialDimPlot(st.f.fibro,group.by = 'distance.tag',
               images = c('P6T'), 
               pt.size.factor = 1.5,
               stroke = 0,
               image.alpha = 0.5,
               ncol=1)
cor.test(st.f.fibro$distanceToEpi,st.f.fibro@assays$Spatial@data['CST1',],
         method = 'spearman')


cds.st.fibro = new_cell_data_set(as(as.matrix(st.f.fibro@assays$Spatial@counts[ rownames(st.f.fibro),]), "sparseMatrix"),
                                 cell_metadata = st.f.fibro@meta.data,
                                 gene_metadata = data.frame(id = rownames(st.f.fibro),gene_short_name=rownames(st.f.fibro),
                                                            row.names = rownames(st.f.fibro))
)
cds.st.fibro[['From']] = factor(cds.st.fibro[['From']],levels = c('mLN','mLN+ PT','mLN- PT'))
plot_genes_violin(cds.st.fibro['CST1',],group_cells_by = 'distance.tag',normalize = T)+
  ggpubr::stat_compare_means()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())
plot_genes_violin(cds.st.fibro[c('TGFB1','TGFBI'),],
                  group_cells_by = 'distance.tag',normalize = T,
                  ncol = 4)+
  ggpubr::stat_compare_means()+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())


saveRDS(st.f.fibro,file = './variables_v2/st.f.fibro.Rds')

st.f.fibro = readRDS(file = './variables_v2/st.f.fibro.Rds')


#####NicheNet#####


cds.sm.fibro[['From']] = factor(cds.sm.fibro[['From']],levels = c('mLN','mLN+ PT','mLN- PT'))
marker_test_res.subCellType = top_markers(cds.sm.fibro, 
                                          group_cells_by="From", 
                                          genes_to_test_per_group=500,
                                          reference_cells=1000, cores=32)
#####Figure 4I#####
mlnpt.res = marker_test_res.subCellType[marker_test_res.subCellType$cell_group == 'mLN+ PT',]


library(nichenetr)
sc.t = st.m[,st.m$Source == 'PLM']
Idents(sc.t)=sc.t$cluster.name

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(file = './data/ligand_target_matrix_nsga2r_final.rds')
weighted_networks <- readRDS(file = "./data/weighted_networks_nsga2r_final.rds")

receiver = "Fibroblast"
expressed_genes_receiver <- get_expressed_genes(receiver, sc.t, pct = 0.1)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()


sender_celltypes <- c('Plasma','Plasma+Myeloid','Fibroblast','Myeloid','Epithelial')

list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, sc.t, 0.25)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 8702
length(potential_ligands)
## [1] 597
length(potential_ligands_focused)
## [1]365

geneset_oi <- as.character(as.data.frame(arrange(mlnpt.res,specificity) %>% top_n(30,specificity))[,'gene_id'])


background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 4567
length(geneset_oi)
## [1] 50
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = order(aupr_corrected,decreasing = T))
ligand_activities


best_upstream_ligands <- ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% tidyr::drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.01) 


order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
vis_ligand_target.row.mean = apply(vis_ligand_target, 1, mean)
vis_ligand_target.row.mean = sort(vis_ligand_target.row.mean,decreasing = T)
vis_ligand_target.high = vis_ligand_target[names(vis_ligand_target.row.mean),]
vis_ligand_target.high = vis_ligand_target.high[,apply(vis_ligand_target.high, 2, sum) >0]
saveRDS(vis_ligand_target,file = './variables_v2/vis_ligand_target_fibro.csv')
make_heatmap_ggplot(vis_ligand_target.high, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  theme(axis.text =  element_text(size = 8))+
  scale_fill_gradientn(colours = rev(hcl.colors(256,'Blue-Red3'))[140:256])
