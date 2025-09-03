
ref.geo = read.delim(file='./variables/GSE131907_Lung_Cancer_cell_annotation.txt')

table(ref.geo$Sample_Origin)
# mBrain    mLN    nLN  nLung     PE   tL/B  tLung 
# 29060  21479  37446  42995  20304  12073  45149 
geo.data = readRDS(file='./data/GSE131907_Lung_Cancer_raw_UMI_matrix.rds')



geo.sc = CreateSeuratObject(counts = geo.data)
geo.sc = AddMetaData(geo.sc,ref.geo)


cyto.genes = read.csv(file = './tables_v3/fig 2d-luad.csv')
cyto.genes = as.character(cyto.genes$X)[1:25]


geo.sc = NormalizeData(geo.sc)
geo.sc = ScaleData(geo.sc)

####CSC####
geo.sc = AddModuleScore(geo.sc,features = list(cyto.genes),name = 'CSC_Score')

ref.geo$CSC_Score = geo.sc$CSC_Score1[as.character(ref.geo$Index)]
ref.geo$FTH1 = geo.sc@assays$RNA@data['FTH1',as.character(ref.geo$Index)]
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


ggplot(plot.data[plot.data$Cell_type == 'Epithelial cells',],aes(x=Sample_Origin,y=log2(FTH1+1),fill=Sample_Origin))+
  geom_violin(adjust=5,scale = 'width')+
  theme_cowplot()+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank())
write.csv(plot.data[plot.data$Cell_type == 'Epithelial cells',c('Index','Cell_type','Sample_Origin','FTH1')],file = './tables_v3/fig 7f_new.csv')

plot.data.1 = plot.data[plot.data$Cell_type == 'Epithelial cells',c('CSC_Score','Sample_Origin','Sample')]
rownames(ref.geo)=ref.geo$Index
geo.sc$Cell_type = ref.geo[colnames(geo.sc),'Cell_type']
geo.sc$Cell_subtype = ref.geo[colnames(geo.sc),'Cell_subtype']
geo.sc$Sample_Origin = ref.geo[colnames(geo.sc),'Sample_Origin']
geo.epi = geo.sc[,geo.sc$Cell_type == 'Epithelial cells',]
Idents(geo.epi)=factor(geo.epi$Sample_Origin,levels = c('mBrain','mLN','nLN','nLung','PE','tL/B','tLung'))
b.csc = DotPlot(geo.epi,
        features = c('FTH1','GAPDH','PKM','LDHA','ENO1','ANXA2','CD44','SOX2'),
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))




####CXCL13 T####

marker_test_res = read.csv(file = './tables_v2/all_refine_markers.csv')
CXCL13.markers = filter(marker_test_res,marker_test_p_value<0.01 & cell_group == 'CXCL13+ T') %>% top_n(20,specificity)
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
plot.data.2 = plot.data[plot.data$Cell_type == 'T lymphocytes',c('CXCL13T_Score','Sample_Origin','Sample')]

geo.t = geo.sc[,geo.sc$Cell_type == 'T lymphocytes']
Idents(geo.t)=factor(geo.t$Sample_Origin,levels = c('mBrain','mLN','nLN','nLung','PE','tL/B','tLung'))
b.t=DotPlot(geo.t,
          features = c('CXCL13','CTLA4','BTLA','TIGIT','SIRPG','TNFRSF18','ICA1','TOX2'),
          cols = c('blue','red'),
)+
  scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))
b.t

####IFI30 MP####

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
plot.data.3 = plot.data[plot.data$Cell_type == 'Myeloid cells',c('IFI30MP_Score','Sample_Origin','Sample')]

geo.myeloid = geo.sc[,geo.sc$Cell_subtype == 'mo-Mac' ,]

Idents(geo.myeloid)=factor(geo.myeloid$Sample_Origin,levels = c('mBrain','mLN','nLN','nLung','PE','tL/B','tLung'))
b.m=DotPlot(geo.myeloid,
          features = c('SPP1','MIF','GABARAP','IFI30','CD274','MMP12','TNFSF13'),
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))
b.m

plot.data.1$Source = 'Epithelial_CSC'
colnames(plot.data.1) = c('Score','Sample_Origin','Source')

plot.data.2$Source = 'CXCL13+ T'
colnames(plot.data.2) = c('Score','Sample_Origin','Source')

plot.data.3$Source = 'SPP1+IFI30+ MP'
colnames(plot.data.3) = c('Score','Sample_Origin','Source')

plot.data.m = rbind(plot.data.1,plot.data.2,plot.data.3)
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
write.csv(plot.data.m,file = './tables_v3/add Fig 7a.csv')

plot.data.m = read.csv(file = './tables_v3/add Fig 7a.csv')

b.csc | b.m |b.t

write.csv(rbind(b.csc$data,b.m$data,b.t$data),file = './tables_v3/Fig S7_validate dotplot GSE131907')




