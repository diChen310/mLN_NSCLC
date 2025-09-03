.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')
library(Seurat)
#
cds = readRDS(file = './step3_anno_cds_GSE123904.Rds')
geo.data = readRDS(file='./variables_v4/step0_1_GSE123904_sc.m.Rds')

cyto.genes = read.csv(file = './tables_v3/fig 2d-luad.csv')
cyto.genes = as.character(cyto.genes$X)[1:25]

geo.data$majorCellType = cds[['majorCellType']]
geo.data$cluster = cds[['cluster']]





geo.sc = NormalizeData(geo.data)
geo.sc = ScaleData(geo.sc)
geo.sc = geo.sc[,geo.sc$orig.ident != 'LX685_N']
geo.sc = AddModuleScore(geo.sc,features = list(cyto.genes),name = 'CSC_Score')

plot.data =geo.sc@meta.data
plot.data1 = plot.data[plot.data$majorCellType == 'Epi',]

plot.data1$From = factor(plot.data1$From,levels = c('METASTASIS','PRIMARY','NORMAL'))
a=ggplot(plot.data1,aes(x=From,y=CSC_Score1,fill=From))+
  geom_violin(adjust=3,scale = 'width')+
  coord_flip()+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
  theme_cowplot(font_size = 8)+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7)[c(5,4,3)])+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank(),
        legend.position = 'none')

a
geo.epi = geo.sc[,geo.sc$majorCellType == 'Epi']
Idents(geo.epi)=factor(geo.epi$From,levels = c('METASTASIS','PRIMARY','NORMAL'))
diff.N1=FindMarkers(geo.epi,ident.1 = 'METASTASIS',ident.2 = 'PRIMARY',only.pos = T,logfc.threshold = 0.01)

intersect(cyto.genes,rownames(diff.N1))

b=DotPlot(geo.epi,
          features = c('ACTB','PKM','FTH1','GAPDH','LDHA','ENO1'),
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 0),
        text = element_text(size = 8))
Idents(geo.epi)=geo.epi$orig.ident
DotPlot(geo.epi,
        features = c('IGHG1','IGHG4','IGHA1','IGLC2'),
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))
a|b


marker_test_res = read.csv(file = './tables_v2/all_refine_markers.csv')

####CXCL13+ T####

CXCL13.markers = filter(marker_test_res,marker_test_p_value<0.05 & cell_group == 'CXCL13+ T') %>% top_n(20,specificity)
CXCL13.markers = as.character(CXCL13.markers$gene_id)
ex.markers = c("CXCL13", "HAVCR2", "PDCD1", "TIGIT", "LAG3", "CTLA4", "LAYN", "RBPJ", "VCAM1", "GZMB", "TOX", "MYO7A")

geo.sc = AddModuleScore(geo.sc,features = list(CXCL13.markers),name = 'CXCL13T_Score')
geo.sc = AddModuleScore(geo.sc,features = list(cc=ex.markers),name = 'exScore')
geo.sc$CXCL13_exp = geo.sc@assays$RNA@counts['CXCL13',]

plot.data2 =geo.sc@meta.data
plot.data2 = plot.data2[plot.data2$majorCellType == 'T' ,]
plot.data2$From = factor(plot.data2$From,levels = c('METASTASIS','PRIMARY','NORMAL'))

a=ggplot(plot.data2,aes(x=From,y=CXCL13T_Score1,fill=From))+
  geom_violin(adjust=3,scale = 'width')+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
  theme_cowplot(font_size = 8)+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7)[c(5,4,3)])+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank(),
        legend.position = 'none')+
  coord_flip()
a

geo.t = geo.sc[,geo.sc$majorCellType == 'T']
Idents(geo.t)=factor(geo.t$From,levels = c('METASTASIS','PRIMARY','NORMAL'))
diff.t.n1 = FindMarkers(geo.t,ident.1 = 'METASTASIS',ident.2 = 'PRIMARY',only.pos = T,logfc.threshold = 0.01)
intersect(CXCL13.markers,rownames(diff.t.n1))
# DoHeatmap(geo.t,features = c('CXCL13','CTLA4','PDCD1','TIGIT','TOX2','BATF','TNFRSF18'),disp.max = 0.2,
#           disp.min = -0.2,
#           group.by = 'N.stage')+
#   scale_fill_gradientn(colours = rev(hcl.colors(100,'RdBu')))
b=DotPlot(geo.t,
          features = c('CXCL13','CTLA4','MAF','TIGIT','SIRPG','TNFRSF18'),
          cols = c('blue','red'),
)+
  scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 0),
        text = element_text(size = 8))
a | b

####IFI30 SPP1 MP####
IFI30.diff = read.csv(file = './tables_v2/IFI30.diff.csv',row.names = 1)
IFI30.top = filter(IFI30.diff,p_val_adj<0.01 & avg_log2FC >0) %>% top_n(10,avg_log2FC)
IFI30.top = as.character(rownames(IFI30.top))



geo.sc = AddModuleScore(geo.sc,features = list(c(IFI30.top,'SPP1')),name = 'IFI30MP_Score')

plot.data3 = geo.sc@meta.data

plot.data3 = plot.data3[plot.data3$majorCellType == 'Myeloid',]
plot.data3$From = factor(plot.data3$From,levels = c('METASTASIS','PRIMARY','NORMAL'))

a=ggplot(plot.data3,
         aes(x=From,y=IFI30MP_Score1,fill=From))+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
  geom_violin(adjust = 3,scale = 'width')+
  coord_flip()+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = ggsci::pal_npg()(7)[c(5,4,3)])+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')
a
geo.myeloid = geo.sc[,geo.sc$majorCellType == 'Myeloid' ,]

Idents(geo.myeloid)=factor(geo.myeloid$From,
                           levels = c('METASTASIS','PRIMARY','NORMAL'))
diff.m.n1 = FindMarkers(geo.myeloid,ident.1 = 'METASTASIS',ident.2 = 'PRIMARY',only.pos = T,
                        logfc.threshold = 0.1)
intersect(IFI30.top,rownames(diff.m.n1))


# Idents(geo.myeloid)=geo.myeloid$orig.ident
# DotPlot(geo.myeloid,
#         features = c('SPP1','NME2','MIF','IGKC','GABARAP','IFI30'),
# )
b=DotPlot(geo.myeloid,
          features = c('SPP1','NME2','MIF','IGKC','GABARAP','IFI30'),
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 0),
        text = element_text(size = 8))

a|b

####Correlations####

library(dplyr)

plot.data1.sum =  plot.data1 %>%
  group_by(orig.ident) %>%          
  summarize(mean_value = mean(CSC_Score1)) 
plot.data2.sum = plot.data2 %>%
  group_by(orig.ident) %>%
  summarize(mean_value=mean(CXCL13T_Score1)) 

plot.data3.sum = plot.data3 %>%
  group_by(orig.ident) %>%
  summarize(mean_value=mean(IFI30MP_Score1))

sum.data = data.frame(CSC=plot.data1.sum$mean_value,
                      CXCL13T=plot.data2.sum$mean_value,
                      IFI30MP=plot.data3.sum$mean_value)
sum.data = read.csv(file = './GSE127465_sumData.csv',row.names = 1)
sum.data = rbind(sum.data,sum.data.2)
corr.sum = cor(as.matrix(sum.data),method = 'spearman')
library(corrplot)

corrplot(corr.sum,method = 'pie',
         ,col = rev(COL2('RdYlBu', 100)),
         type = 'upper',
         is.corr = T,addCoef.col = 'black')
write.csv(sum.data,file = './GSE127465&GSE123904&GSE131907_sumData.csv')
sum.data = read.csv(file = './GSE127465&GSE123904&GSE131907_sumData.csv',row.names = 1)

save.image(file = 'GSE123904Validation.RData')

