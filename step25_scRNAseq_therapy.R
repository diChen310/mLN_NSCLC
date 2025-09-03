sc.data = readRDS(file = '/share/dichen/data/zendo_14728962/NSCLC-singlecell-all_new.rds')
marker_test_res = read.csv(file = './tables_v2/all_refine_markers.csv')
CXCL13.markers = filter(marker_test_res,marker_test_p_value<0.01 & cell_group == 'CXCL13+ T') %>% top_n(20,specificity)
CXCL13.markers = as.character(CXCL13.markers$gene_id)

sc.data = NormalizeData(sc.data)
sc.data = AddModuleScore(sc.data,features = list(c('BTLA','CD2', "CD40LG" , "CTLA4" ,   "CXCL13"  , "GNG4"   , "ICA1"   ,  "KLRB1"  ,   "LINC01281", "MAF"   , 
                                                   "PDCD1"  ,   "SIRPG"  ,  "SLA"   ,    "TIGIT"  ,   "TNFRSF18" , "TOX2")),name = 'CXCL13T_Score',nbin=50)


plot.data = sc.data@meta.data

ggplot(plot.data[plot.data$lineage == 'T',],aes(x=Group,y=CXCL13T_Score1,fill=Group))+
  geom_boxplot()+  theme_cowplot()+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(2,3),c(3,4),c(1,3),c(2,4),c(1,4)))+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank())

VlnPlot(sc.data[,sc.data$lineage == 'T'],group.by = 'Group',features = c('CXCL13','BTLA','CD2','CTLA4'),pt.size = 0)

write.csv(plot.data[plot.data$lineage == 'T',],file = './tables_v3/revise_YanDataset_CXCL13TScores.csv')
####IFI30 MP####

IFI30.diff = read.csv(file = './tables_v2/IFI30.diff.csv',row.names = 1)
IFI30.top = filter(IFI30.diff,p_val_adj<0.01 & avg_log2FC >0) %>% top_n(10,avg_log2FC)
IFI30.top = as.character(rownames(IFI30.top))


sc.data = AddModuleScore(sc.data,features = list(c(IFI30.top,'SPP1')),name = 'IFI30MP_Score',nbin = 50)

plot.data = sc.data@meta.data

ggplot(plot.data[plot.data$major %in% c('Macrophage'),],aes(x=Group,y=IFI30MP_Score1,fill=Group))+
  geom_boxplot()+  theme_cowplot()+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(2,3),c(3,4),c(1,3),c(2,4),c(1,4)))+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank())
write.csv(plot.data[plot.data$major %in% c('Macrophage'),],file = './tables_v3/revise_YanDataset_IFI30MPScores.csv')