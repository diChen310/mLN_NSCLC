.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')

####Read GSE131907 data#####

library(Seurat)
ref.geo = read.delim(file='./variables/GSE131907_Lung_Cancer_cell_annotation.txt')

table(ref.geo$Sample_Origin)
# mBrain    mLN    nLN  nLung     PE   tL/B  tLung 
# 29060  21479  37446  42995  20304  12073  45149 
geo.data = readRDS(file='./variables/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds')



####Prepare our data####
cds.all = readRDS(file = './variables_v2/inter_cds_refine.Rds')
sel.cells = c()

for(cell.i in unique(colData(cds.all)$subCellType)){
  cells.i = colnames(cds.all[,cds.all[['subCellType']] == cell.i])
  cells.i.s = sample(cells.i,length(cells.i)*0.25)
  sel.cells = append(sel.cells,cells.i.s)
}
length(sel.cells)# 212062
####Single R####
library(SingleR)
pred.res.i = SingleR(test=as.matrix(geo.data),
                     ref=as.matrix(assay(cds.all[,sel.cells])), 
                     labels=cds.all[,sel.cells][['subCellType']], 
                     de.method="wilcox")
print(Sys.time())
saveRDS(pred.res.i,file = './variables_v2/singleR gse131907 prediction.Rds')

ref.geo$pred.type = pred.res.i$pruned.labels

library(ggplot2)
library(cowplot)
library(ggsci)
library(Seurat)
#####

rm(pred.res.i)


geo.sc = CreateSeuratObject(counts = geo.data)
geo.sc = AddMetaData(geo.sc,ref.geo)
  
  
####Stemness Scores####

cyto.LUAD = readRDS(file = './variables/cytoTRACE_LUAD.Rds')
cyto.genes.luad = names(cyto.LUAD$cytoGenes[cyto.LUAD$cytoGenes>0.549])

#geo.sc = NormalizeData(geo.sc)
geo.sc = AddModuleScore(geo.sc,features = list(cyto.genes.luad),name = 'CSC_Score')

ref.geo$CSC_Score = geo.sc$CSC_Score1[as.character(ref.geo$Index)]
plot.data =ref.geo
plot.data$Sample_Origin = factor(plot.data$Sample_Origin,levels = c('mBrain','mLN','nLN','nLung','PE','tL/B','tLung'))
ggplot(plot.data[plot.data$Cell_type == 'Epithelial cells',],aes(x=Sample_Origin,y=CSC_Score,fill=Sample_Origin))+
  geom_violin(adjust=3,scale = 'width')+
  theme_cowplot()+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank())

plot.data.1 = plot.data[plot.data$Cell_type == 'Epithelial cells',c('CSC_Score','Sample_Origin')]
marker_test_res = read.csv(file = './tables_v2/all_refine_markers.csv')



####FTH1####

ref.geo$FTH1 = as.numeric(geo.data['FTH1',as.character(ref.geo$Index)])
plot.data =ref.geo
plot.data$Sample_Origin = factor(plot.data$Sample_Origin,levels = c('mBrain','mLN','nLN','nLung','PE','tL/B','tLung'))
ggplot(plot.data[plot.data$Cell_type == 'Epithelial cells',],aes(x=Sample_Origin,y=FTH1,fill=Sample_Origin))+
  geom_violin(adjust=3,scale = 'width')+
  theme_cowplot(font_size = 8)+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank())+
  ylab('Log2TPM_FTH1')
write.csv(plot.data[plot.data$Cell_type == 'Epithelial cells',c('Sample_Origin','FTH1')],file = './tables_v3/add Fig 7i.csv')

ref.geo$IGHG4 = as.numeric(geo.data['IGHG4',as.character(ref.geo$Index)])
ref.geo$IGHG3 = as.numeric(geo.data['IGHG3',as.character(ref.geo$Index)])
ref.geo$ALDOA = as.numeric(geo.data['ALDOA',as.character(ref.geo$Index)])

plot.data =ref.geo
plot.data$Sample_Origin = factor(plot.data$Sample_Origin,levels = c('mBrain','mLN','nLN','nLung','PE','tL/B','tLung'))
ggplot(plot.data[plot.data$Cell_type == 'T lymphocytes',],aes(x=Sample_Origin,y=ALDOA,fill=Sample_Origin))+
  geom_violin(adjust=5,scale = 'width')+
  theme_cowplot(font_size = 8)+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank())

####IG scores####

genes = c('IGHG1','IGHG4','IGHA1','IGLC2')
geo.sc = AddModuleScore(geo.sc,features = list(genes),name = 'Ig_Score')

ref.geo$Ig_Score = as.numeric(geo.data['IGHG4',ref.geo$Index])
plot.data =ref.geo

ggplot(plot.data[plot.data$Cell_type == 'Epithelial cells',],aes(x=Sample_Origin,y=log2(Ig_Score+1),fill=Sample_Origin))+
  geom_violin(adjust=5,scale = 'width')+
  theme_cowplot()+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank())
plot.data.2 = plot.data[plot.data$Cell_type == 'Epithelial cells',c('Ig_Score','Sample_Origin')]

####CST1 marker genes####
CST1.markers = filter(marker_test_res,marker_test_p_value<0.05 & cell_group == 'CST1+ Fibro') %>% top_n(20,specificity)
CST1.markers = as.character(CST1.markers$gene_id)

geo.sc = AddModuleScore(geo.sc,features = list(CST1.markers),name = 'CST1_Score')

ref.geo$CST1_Score = geo.sc$CST1_Score1[as.character(ref.geo$Index)]
plot.data =ref.geo

ggplot(plot.data[plot.data$Cell_type == 'Fibroblasts',],aes(x=Sample_Origin,y=CST1_Score,fill=Sample_Origin))+
  geom_violin(adjust=3,scale = 'width')+
  theme_cowplot()+
  stat_summary(fun = mean, geom = "point", size = 0.5,
                          color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank())
plot.data.3 = plot.data[plot.data$Cell_type == 'Fibroblasts',c('CST1_Score','Sample_Origin')]
####CXCL13 marker genes####
CXCL13.markers = filter(marker_test_res,marker_test_p_value<0.05 & cell_group == 'CXCL13+ T') %>% top_n(20,specificity)
CXCL13.markers = as.character(CXCL13.markers$gene_id)

geo.sc = AddModuleScore(geo.sc,features = list(CXCL13.markers),name = 'CXCL13T_Score')

ref.geo$CXCL13T_Score = geo.sc$CXCL13T_Score1[as.character(ref.geo$Index)]
plot.data =ref.geo


ggplot(plot.data[plot.data$Cell_type == 'T lymphocytes',],
       aes(x=Sample_Origin,y=CXCL13T_Score,fill=Sample_Origin))+
geom_violin(adjust = 3,scale = 'width')+
theme_cowplot(font_size = 8)+
scale_fill_manual(values = ggsci::pal_npg()(7))+
stat_summary(fun = mean, geom = "point", size = 0.5,
color = "black")+
theme(axis.text.x = element_text(angle = 45,hjust = 1),
strip.background = element_blank(),
panel.grid = element_blank())
plot.data.4 = plot.data[plot.data$Cell_type == 'T lymphocytes',c('CXCL13T_Score','Sample_Origin')]
####IFI30 markers#####
IFI30.diff = read.csv(file = './tables_v2/IFI30.diff.csv',row.names = 1)
IFI30.top = filter(IFI30.diff,p_val_adj<0.01 & avg_log2FC >0) %>% top_n(10,avg_log2FC)
IFI30.top = as.character(rownames(IFI30.top))



geo.sc = AddModuleScore(geo.sc,features = list(c(IFI30.top,'SPP1')),name = 'IFI30MP_Score')

ref.geo$IFI30MP_Score = geo.sc$IFI30MP_Score1[as.character(ref.geo$Index)]
plot.data =ref.geo


ggplot(plot.data[plot.data$Cell_type == 'Myeloid cells',],
       aes(x=Sample_Origin,y=IFI30MP_Score,fill=Sample_Origin))+
  geom_violin(adjust = 3,scale = 'width')+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        panel.grid = element_blank())
plot.data.5 = plot.data[plot.data$Cell_type == 'Myeloid cells',c('IFI30MP_Score','Sample_Origin')]

####Merge####

plot.data.1$Source = 'Epithelial_CSC'
colnames(plot.data.1) = c('Score','Sample_Origin','Source')

plot.data.3$Source = 'Fibroblast_CST1'
colnames(plot.data.3) = c('Score','Sample_Origin','Source')

plot.data.4$Source = 'CXCL13+ T'
colnames(plot.data.4) = c('Score','Sample_Origin','Source')

plot.data.5$Source = 'SPP1+IFI30+ MP'
colnames(plot.data.5) = c('Score','Sample_Origin','Source')

plot.data.m = rbind(plot.data.1,plot.data.4,plot.data.5)
plot.data.m$Sample_Origin = factor(plot.data.m$Sample_Origin,levels = c('mBrain','mLN','nLung','PE','tL/B','tLung','nLN'))
plot.data.m$Source = factor(plot.data.m$Source,levels = c('Epithelial_CSC','SPP1+IFI30+ MP','CXCL13+ T'))
ggplot(plot.data.m,
       aes(x=Sample_Origin,y=Score,fill=Sample_Origin))+
  facet_wrap(~Source,scales = 'free_y',ncol = 1)+
  geom_violin(adjust = 2,scale = 'width')+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        panel.grid = element_blank())

signature.gene = marker_test_res.sig[,c('gene_id','cell_group','specificity','marker_test_p_value')]
write.csv(signature.gene,file = './tables_v3/Signature subCellType markers.csv')

sig.gene = data.frame(Gene = c(both.cyto.genes,CXCL13.markers,IFI30.top,'SPP1'),
                      Type = c(rep('CSC',length(both.cyto.genes)),
                               rep('CXCL13+ T',length(CXCL13.markers)),
                               rep('IFI30+SPP1+',length(IFI30.top)+1)))
write.csv(sig.gene,file = './tables_v3/GSE131907 signature genes used.csv')
