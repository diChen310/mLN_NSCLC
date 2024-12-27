.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')
library(Seurat)
library(monocle3)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggsci)
#### Read previous results####
cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559","#32709f",
                 '#b8a3c4','#60428d','#9f5633','grey','#c92e2a','#c79795','#d9c17c',"#d0896b",'#5d94b4',"black","blue",'cyan')

cds.all = readRDS(file = './variables_v2/inter_cds.Rds')

Myeloid.data = readRDS(file = './variables_v2/myeloid.data.refine.Rds')
table(Myeloid.data[['subCellType']])
B.data = readRDS(file = './variables_v2/B.data.Rds')
table(B.data[['subCellType']])

T.data = readRDS(file='./variables_v2/T.data.Rds')
table(T.data[['subCellType']])

stromal.data <- readRDS(file = './variables_v2/stromal.data.refine.Rds')
table(stromal.data[['subCellType']])

epi.data = readRDS(file = './variables_v2/epi.data.Rds')
epi.data = epi.data[,epi.data[['subCellType']] != 'Remove']

lym.data = readRDS(file = './variables_v2/lym.data.Rds')
lym.proliferation = lym.data[,lym.data[['subCellType']] == 'Proliferation lymphocyte']
table(lym.proliferation[['subCellType']])
table(epi.data[['subCellType']])



cds.all[['subCellType']] = 'Remove'

meta.data = as.data.frame(colData(cds.all))
meta.data[colnames(epi.data),'subCellType'] = epi.data[['subCellType']]
meta.data[colnames(Myeloid.data),'subCellType'] = as.character(Myeloid.data[['subCellType']])
meta.data[colnames(stromal.data),'subCellType'] = stromal.data[['subCellType']]
meta.data[colnames(T.data),'subCellType'] = T.data[['subCellType']]
meta.data[colnames(B.data),'subCellType'] = as.character(B.data[['subCellType']])
meta.data[colnames(lym.proliferation),'subCellType'] = lym.proliferation[['subCellType']]
meta.data[meta.data$majorCellType == 'Mast','subCellType'] = 'Mast'
table(meta.data$subCellType)



cds.all[['subCellType']] = meta.data$subCellType

cds.all[['majorCellType2']] = recode(cds.all[['subCellType']],
                                     "CD4 Naive T"  ='T',
                                     "CD8 Exhausted T" ='T',
                                     "AT2" ='Epi',
                                     "B" ='B',
                                     "CD4 Resident Effector Memory" = 'T',
                                     "CD8 Cytotoxic T"  ='T',
                                     "APOE+ MP"  ='Myeloid',
                                     "cDC2"   ='Myeloid',
                                     "ALDOA+ T"='T',
                                     "Gran"  ='Myeloid',     
                                     "CD4 Effector Memory"   ='T',    
                                     "Club"     ='Epi'       ,
                                     "Ciliated" ='Epi',
                                     "SPP1+ MP" ='Myeloid',
                                     "Treg"='T', 
                                     "LEC" ='Endo',
                                     "Naive B"  ='B',
                                     "¦Ã¦Ä T"='T',
                                     "Activated B"   = 'B',               
                                     "Veins"  ='Endo',
                                     "cDC1" ='Myeloid', 
                                     "mCAF" ='Fibro',
                                     "LAMP3+ DC"   ='Myeloid',
                                     "FABP4+ MP" ='Myeloid',  
                                     "PIP+ Epi"   ='Epi',                 
                                     "Basal"='Epi',
                                     "CXCL13+ T"  ='T',
                                     "CD1C- DC"   ='Myeloid',
                                     "ZBTB20+ B"  ='B',
                                     "GC B"='B',
                                     "pDC" ='Myeloid', 
                                     "Arteries" ='Endo',
                                     "Exhausted T"='T',
                                     "SUCNR1+ Epi"  ='Epi',
                                     "AT1"='Epi',
                                     "iCAF" ='Fibro',
                                     "Capillaries" ='Endo',
                                     "CST1+ Fibro"='Fibro',
                                     "SMC"='Fibro', 
                                     "ARHGAP15+ T" ='T',
                                     "Aerocyte"='Endo',
                                     "tCAF" ='Fibro',
                                     "vCAF" ='Fibro', 
                                     "Endo progenitor"='Endo',
                                     "Pericyte"='Fibro',
                                     "Neu" ='Myeloid',
                                     "IgG low plasma" ='Plasma')

cds.all = cds.all[,cds.all[['subCellType']] != 'Remove'] # 14787 848337

cols.cl =cols.default[1:length(unique(colData(cds.all)$majorCellType))]
names(cols.cl)=c('T','B','Plasma','Epi','Mast','Myeloid','Fibro','Endo','NK')

scater::plotReducedDim(cds.all,dimred = 'UMAP',colour_by = 'majorCellType',point_size = 1)+
  ggplot2::scale_color_manual(values = cols.cl)+
  theme_void()+theme(legend.position = 'none')

saveRDS(cds.all,file = './variables_v2/inter_cds_refine.Rds')
cds.all = readRDS(file = './variables_v2/inter_cds_refine.Rds')##########################re-analysis need!!!!!!!!!!!!!!!!!!!!!!!

cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')
plot_cells(cds.all,  color_cells_by="From",
           cell_size=0.1,
           cell_stroke=0,
           rasterize = T,
           label_cell_groups = F)+
  ggplot2::scale_color_manual(values = cols.From)+
  facet_wrap(~From, ncol = 3)+
  theme_void()+
  theme( legend.position = 'none')
cols.sources = jcolors('pal8')
names(cols.sources)=unique(cds.all[['Source']])
plot_cells(cds.all, color_cells_by="Source", label_cell_groups=FALSE,cell_size = 0.1,cell_stroke = 0)+
  facet_wrap(~Source,ncol = 4)+
  scale_color_manual(values = cols.sources)+
  theme_void()+theme(legend.position = 'none')


meta.data.show = cds.all@colData
meta.data.show = as.data.frame(meta.data.show)
meta.data.show$Disease2 = substr(meta.data.show$Disease,1,4)
write.csv(meta.data.show,file = './tables_v2/meta data all_subCellTypes.csv')
meta.data.show = read.csv(file = './tables_v2/meta data all_subCellTypes.csv',row.names = 1)

ggplot2::ggplot(meta.data.show,ggplot2::aes(x=From,fill=majorCellType))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values = cols.cl)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

com2orig = table(meta.data.show$orig.ident,meta.data.show$majorCellType)
f.data<-as.matrix(com2orig)
f.data = f.data/apply(f.data, 1, sum)
#f.data = t(f.data)
#rownames(ids.all)=ids.all$id
meta.data.u = meta.data.show
meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
rownames(meta.data.u)=meta.data.u$orig.ident

com2orig.df = melt(f.data,measure.vars = colnames(f.data))
com2orig.df = com2orig.df[com2orig.df$Var2 != 'Unassigned',]##zero values

com2orig.df$Var2 = as.character(com2orig.df$Var2)
com2orig.df$From = factor(meta.data.u[as.character(com2orig.df$Var1),'From'],levels = 
                            c('mLN+ N','mLN','mLN+ PT','mLN- N','mLN- PT','nLN'))
com2orig.df$Disease = substr(meta.data.u[as.character(com2orig.df$Var1),'Disease'],1,4)
ggplot(com2orig.df,aes(x=From,y=log2(value+1),fill=From))+
  ggpubr::stat_compare_means(comparisons = list(c(2,5),c(3,5)))+
  geom_boxplot()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))+
  scale_fill_manual(values = cols.From)


ggplot(com2orig.df[com2orig.df$Var2 == 'Plasma',],aes(x=From,y=value,fill=From))+
  stat_ydensity(draw_quantiles = 0.5,scale = 'width',adjust = 5)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))+
  scale_fill_manual(values = cols.From)



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
plot_genes_by_group(cds.all,
                    marker.genes.classical,
                    group_cells_by="majorCellType",
                    ordering_type="cluster_row_col",
                    max.size=3)+scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))

meta.data.u = meta.data.show
meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
rownames(meta.data.u)=meta.data.u$orig.ident

com2orig.L = table(meta.data.show$orig.ident,meta.data.show$subCellType)

f.data<-as.matrix(com2orig.L)
f.data = f.data/apply(f.data, 1, sum)

pheatmap(f.data,show_rownames = F,
         annotation_row = data.frame(Disease = meta.data.u[rownames(f.data),'Disease2'],
                                     From = meta.data.u[rownames(f.data),'From'],
                                     Stage = meta.data.u[rownames(f.data),'Stage'],
                                     row.names = rownames(f.data)),
         clustering_distance_cols = 'correlation')


corrs.cell = cor(f.data,method = 'spearman')

pheatmap(corrs.cell,color = colorRampPalette(rev(hcl.colors(100,'Spectral')))(100))
library(corrplot)
corrplot(corrs.cell,method = 'pie',type = 'upper',col=rev(COL2('RdBu',200)),
         order='hclust',
         hclust.method='ward',
         diag=F,pch=2,
         tl.cex = .5
)

####PLM based correlations####

cds = cds.all[,cds.all[['Source']] == 'PLM']#14787 334956##########################re-analysis need!!!!!!!!!!!!!!!!!!!!!!!


meta.data.show.L = meta.data.show[meta.data.show$Source == 'PLM',]
com2orig.L = table(meta.data.show.L$orig.ident,meta.data.show.L$subCellType)

f.data<-as.matrix(com2orig.L)
f.data = f.data/apply(f.data, 1, sum)



corrs.cell = cor(f.data,method = 'spearman')

pheatmap(corrs.cell,color = colorRampPalette(rev(hcl.colors(100,'Spectral')))(100))


##Epi: malignant or not
epi.data.mal = readRDS(file = './variables_v2/epi.data.mal.Rds')
epi.data.mal.LUSC = epi.data.mal[,epi.data.mal[['PatientID']] %in% c('P0','P1','P5','P9','P15')]
epi.data.mal.LUAD = epi.data.mal[,epi.data.mal[['PatientID']] %in% c('P0','P1','P5','P9','P15') ==F]

meta.data.cds = as.data.frame(colData(cds))

meta.data.cds$subCellType_v2 = meta.data.cds$subCellType
meta.data.cds[colnames(epi.data.mal.LUAD),'subCellType_v2'] = 'Malignant Epi'
meta.data.cds[colnames(epi.data.mal.LUSC),'subCellType_v2'] = 'Malignant Epi'
meta.data.cds[meta.data.cds$majorCellType == 'Epi' & (rownames(meta.data.cds) %in% colnames(epi.data.mal) == F),'subCellType_v2'] = 'Non-malignant Epi'
write.csv(meta.data.cds,file = './tables_v2/meta.data.cds with malignant or not Epi.csv') ###Source = PLM only
meta.data.cds = read.csv(file = './tables_v2/meta.data.cds with malignant or not Epi.csv',row.names = 1)##########################re-analysis need!!!!!!!!!!!!!!!!!!!!!!!
cds[['subCellType_v2']] = meta.data.cds$subCellType_v2##########################re-analysis need!!!!!!!!!!!!!!!!!!!!!!!

com2orig.2 = table(meta.data.cds$orig.ident,meta.data.cds$subCellType_v2)

f.data<-as.matrix(com2orig.2)
f.data = f.data/apply(f.data, 1, sum)

corrs.cell = cor(f.data,method = 'spearman')

pheatmap(corrs.cell,color = colorRampPalette(rev(hcl.colors(100,'Spectral')))(100))
write.csv(corrs.cell,file = './tables_v2/subCellType spearman correlations PLM.csv')
corr.plot = data.frame(Malignant=f.data[,'Malignant Epi'],
                       SPP1_MP = f.data[,'SPP1+ MP'])
ggscatter(corr.plot, y = "Malignant", x = "SPP1_MP",color = 'SPP1_MP',size=1,
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
)+
  stat_cor()+
  scale_color_gradientn(colors = hcl.colors(100))+
  theme_cowplot(font_size = 8)

library(corrplot)

corrs.cell.t = cor(f.data[,c(colnames(f.data)[grepl(' MP',colnames(f.data))],'Malignant Epi','Non-malignant Epi')],metho='spearman')
corrplot(corrs.cell.t,method = 'pie',type = 'lower',col=rev(COL2('RdBu',200)),
         order='hclust',
         hclust.method='ward',
         diag=F,pch=2,
         tl.cex = .5
)


corrs.cell.t = cor(f.data[,c('Malignant Epi','CST1+ Fibro','iCAF','mCAF','Pericyte','SMC','vCAF','tCAF',
                             'Aerocyte','Arteries','Capillaries','Endo progenitor','LEC','Veins')],metho='spearman')
corrplot(corrs.cell.t,method = 'pie',type = 'lower',col=rev(COL2('RdBu',200)),
         order='hclust',
         hclust.method='ward',
         diag=F,pch=2,
         tl.cex = .5
)
# meta.data.cds.mLN= meta.data.cds[meta.data.cds$From == 'mLN',]
# com2orig.3 = table(meta.data.cds.mLN$orig.ident,meta.data.cds.mLN$subCellType_v2)
# 
# f.data<-as.matrix(com2orig.3)
# f.data = f.data/apply(f.data, 1, sum)
# 
# corrs.cell.mLN = cor(f.data,method = 'spearman')
# 
# pheatmap(corrs.cell.mLN,color = colorRampPalette(rev(hcl.colors(100,'Spectral')))(100))
# 
write.csv(corrs.cell.mLN,file = './tables_v2/subCellType spearman correlations mLN.csv')

####cellChat analysis####
cds.sub = cds[,!grepl('Proliferation',cds[['subCellType']]) & cds[['From']] != 'mLN+ N']

subCellType.count = table(cds.sub[['subCellType_v2']])
subCellType.count = subCellType.count[order(subCellType.count)]

input.subCellType.not = names(subCellType.count[subCellType.count<100])

input.subCellType = names(subCellType.count)[names(subCellType.count) %in% input.subCellType.not == F]

cds.sub = cds[,cds[['subCellType_v2']] %in% input.subCellType & cds[['From']] != 'mLN+ N']

dim(cds.sub)
saveRDS(cds.sub,file = './variables_v2/cds.PLM.trimForCellChat.Rds')

library(CellChat)
library(patchwork)
library(monocle3)
library(pheatmap)

cds.PLM = readRDS(file = './variables_v2/cds.PLM.trimForCellChat.Rds')
options(stringsAsFactors = FALSE)
#####Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# # use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

future::plan("multiprocess", workers = 16) # do parallel


cds.PLM.T = cds.PLM[,cds.PLM[['From']] == "mLN" ]

data.input = cds.PLM.T@assays@data$counts
meta = as.data.frame(colData(cds.PLM.T))
#
cellchat.TI = createCellChat(object = data.input,meta,group.by = 'subCellType_v2')

levels(cellchat.TI@idents)

# [1] "Activated B"                  "ALDOA+ T"                     "APOE+ MP"                    
# [4] "ARHGAP15+ T"                  "Arteries"                     "Capillaries"                 
# [7] "CD1C- DC"                     "CD4 Effector Memory"          "CD4 Naive T"                 
# [10] "CD4 Resident Effector Memory" "CD8 Cytotoxic T"              "CD8 Exhausted T"             
# [13] "cDC1"                         "cDC2"                         "CST1+ Fibro"                 
# [16] "CXCL13+ T"                    "Endo progenitor"              "Exhausted T"                 
# [19] "FABP4+ MP"                    "GC B"                         "Gran"                        
# [22] "iCAF"                         "LAMP3+ DC"                    "LEC"                         
# [25] "Malignant Epi"                "Mast"                         "mCAF"                        
# [28] "Memory B"                     "Naive B"                      "Neu"                         
# [31] "NK"                           "Non-malignant Epi"            "pDC"                         
# [34] "Plasma"                       "SMC"                          "SPP1+ MP"                    
# [37] "Treg"                         "vCAF"                         "Veins"                       
# [40] "ZBTB20+ B"                    "¦Ã¦Ä T"                             


options(future.globals.maxSize = 1000 * 1024^2)

####Preprocessing the expression data for cell-cell communication analysis
# set the used database in the object
cellchat.TI@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat.TI <- CellChat::subsetData(cellchat.TI) # This step is necessary even if using the whole database
cellchat.TI <- identifyOverExpressedGenes(cellchat.TI)
cellchat.TI <- identifyOverExpressedInteractions(cellchat.TI)
# project gene expression data onto PPI network (optional)
cellchat.TI <- projectData(cellchat.TI, PPI.human)

####Compute the communication probability and infer cellular communication network
cellchat.TI <- computeCommunProb(cellchat.TI)
cellchat.TI <- computeCommunProbPathway(cellchat.TI)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.TI <- filterCommunication(cellchat.TI, min.cells = 20)

cellchat.TI <- aggregateNet(cellchat.TI)
groupSize <- as.numeric(table(cellchat.TI@idents))

pheatmap(cellchat.TI@net[['count']])
pheatmap(cellchat.TI@net[['weight']])


netVisual_bubble(cellchat.TI, sources.use =  c(17,34), targets.use =c(17,34), remove.isolate = T)
netVisual_bubble(cellchat.TI, sources.use =  c(26,36), targets.use =c(26,36), remove.isolate = T)
netVisual_bubble(cellchat.TI, sources.use =  c(17),  remove.isolate = T)
netVisual_bubble(cellchat.TI, targets.use =  c(17),  remove.isolate = T)
netVisual_bubble(cellchat.TI, sources.use =  c(17),  remove.isolate = T,signaling = c('COLLAGEN','CXCL','CCL','FN1'))

netVisual_chord_gene(cellchat.TI, targets.use =  c(17),  sources.use =  c(23,28),signaling = c('COLLAGEN','CXCL','CCL','FN1'))
netVisual_chord_gene(cellchat.TI, targets.use =  c(6,33),  sources.use =  c(17))

netVisual_bubble(cellchat.TI, sources.use =  c(17,34), targets.use =c(17,34), remove.isolate = T)



saveRDS(cellchat.TI,file = './variables_v2/cellchat.mLNPT.Rds')

saveRDS(cellchat.TI,file = './variables_v2/cellchat.mLN.Rds')


cellchat.TI.pt = readRDS(file = './variables_v2/cellchat.mLNPT.Rds')
cellchat.TI = readRDS(file = './variables_v2/cellchat.mLN.Rds')

#####mLN+ PT#####
#####Plasma
netVisual_bubble(cellchat.TI.pt,sources.use = c(34))
netVisual_bubble(cellchat.TI.pt,targets.use = c(34))
netVisual_bubble(cellchat.TI.pt,targets.use = c(34),signaling = c('MIF','CXCL'))
netVisual_circle(cellchat.TI.pt@net$weight, targets.use =  c(34), sources.use = c(3,5,7,14,15,16,18,20,22,23,24,25,26,27,28,38,39),
                 vertex.weight =1, vertex.label.cex = 0.5,remove.isolate = T,
                 weight.scale = T, label.edge= F, title.name =  "Interaction weights/strength")

#####SPP1+ MP
netVisual_bubble(cellchat.TI.pt,sources.use = c(36))
netVisual_bubble(cellchat.TI.pt,targets.use = c(36))
netVisual_chord_gene(cellchat.TI.pt,sources.use = c(36),targets.use = c(25))

#####Malignant and Fibro
netVisual_bubble(cellchat.TI.pt, targets.use =  c(25))
netVisual_chord_gene(cellchat.TI.pt, targets.use =  c(25),  sources.use =  c(38,27),signaling = c('COLLAGEN','CXCL','CCL','FN1'))

####CXCL13+ T
netVisual_bubble(cellchat.TI.pt, sources.use =  c(16))
netVisual_bubble(cellchat.TI.pt, sources.use =  c(16),signaling = c('CXCL','CCL'))

####mCAF

#####mLN####

#####Malignant and Fibro
CellChat::netVisual_bubble(cellchat.TI, targets.use =  c(25))
CellChat::netVisual_chord_gene(cellchat.TI.pt, targets.use =  c(25),  sources.use =  c(38,27),signaling = c('COLLAGEN','CXCL','CCL','FN1'))


####NicheNet analysis####

library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
cds.t = cds[,cds[['From']] == 'mLN+ PT']
sc.t = CreateSeuratObject(counts = assay(cds.t))
sc.t = AddMetaData(sc.t,as.data.frame(colData(cds.t)))
Idents(sc.t)=sc.t$subCellType_v2

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(file = './data/ligand_target_matrix_nsga2r_final.rds')
weighted_networks <- readRDS(file = "./data/weighted_networks_nsga2r_final.rds")

receiver = "Malignant Epi"
expressed_genes_receiver <- get_expressed_genes(receiver, sc.t, pct = 0.1)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()


sender_celltypes <- unique(sc.t$subCellType_v2)[unique(sc.t$subCellType_v2) %in% receiver ==F]

list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, sc.t, 0.1)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 8702
length(potential_ligands)
## [1] 597
length(potential_ligands_focused)
## [1]365


cyto.genes.luad = names(cyto.LUAD$cytoGenes[abs(cyto.LUAD$cytoGenes)>0.4])
cyto.genes.lusc = names(cyto.LUSC$cytoGenes[abs(cyto.LUSC$cytoGenes)>0.4])


#geneset_oi <- intersect(cyto.genes.luad,cyto.genes.lusc)

geneset_oi <- c('IGHA1','IGHG4','IGLC2','IGHG1')

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 5917
length(geneset_oi)
## [1] 46




ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = order(aupr_corrected,decreasing = T))
ligand_activities


best_upstream_ligands <- ligand_activities %>% top_n(100, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 1500) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.0) 


order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

plot_genes_by_group(cds,
                    order_ligands,
                    group_cells_by="majorCellType",
                    ordering_type="cluster_row_col",
                    max.size=3)+scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))


#####receiver, malignant cells#####
receiver = "Malignant Epi"
expressed_genes_receiver <- get_expressed_genes(receiver, sc.t, pct = 0.1)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()


sender_celltypes <- unique(sc.t$subCellType_v2)[unique(sc.t$subCellType_v2) %in% receiver ==F]

list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, sc.t, 0.1)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 8702
length(potential_ligands)
## [1] 597
length(potential_ligands_focused)
## [1]365

cyto.genes.luad = names(cyto.LUAD$cytoGenes)[1:50]
cyto.genes.lusc = names(cyto.LUSC$cytoGenes)[1:50]


both.cyto.genes = union(cyto.genes.luad,cyto.genes.lusc)#46 genes

geneset_oi <- unique(both.cyto.genes)


background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 5917
length(geneset_oi)
## [1] 46




ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities


best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 


order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")


