
#saveRDS(sc.m,file='./step0_1_GSE127465_sc.m.Rds')
library(GEOquery)
data.info = getGEO(GEO = 'GSE127465')

data.ph = phenoData(data.info[[1]])

data.ph = data.ph@data
library(data.table)
library(scDblFinder)
meta.data = fread('/share/dichen/data/GSE127465/barcodes.tsv.gz')
meta.data = as.data.frame(meta.data)
rownames(meta.data)=paste0(meta.data$Library,'_',meta.data$Barcode)

sc.input.folder = '/share/dichen/data/GSE127465/'
sc.files = list.files(path = sc.input.folder,recursive = F)
sc.files = sc.files[grepl('human',sc.files)]

sc.objs = c()
cell_ids = c()

for(sc.file in sc.files){
  mat.i = fread(paste0(sc.input.folder,sc.file))
  id.i = strsplit(sc.file,'_')[[1]][3]
  expr.d = as(t(mat.i[,-1]), "sparseMatrix")
  colnames(expr.d)=mat.i$barcode
  sc.obj = CreateSeuratObject(expr.d, names.field = 1,min.cells = 3)
  sc.obj$orig.ident = id.i
  #sc.obj$From = type.i
  cell_ids = append(cell_ids,id.i)
  sc.objs = append(sc.objs,sc.obj)
}

print(cell_ids)
sc.m = merge(sc.objs[1][[1]],sc.objs[-1],add.cell.ids = cell_ids)
rownames(data.ph)= unlist(lapply(data.ph$title, function(a){
  id.i = gsub(']','',strsplit(a,'_')[[1]][4])
}))
sc.m$From = data.ph[sc.m$orig.ident,'source_name_ch1']
sc.m$stage = data.ph[sc.m$orig.ident,'characteristics_ch1.5']
sc.m$stage = gsub('tumor stage, ajcc 8th edition:','',sc.m$stage)
sc.m$cellType.pre = meta.data[colnames(sc.m),'Major cell type']
sc.m$tissue = ifelse(grepl('b',sc.m$orig.ident),'blood','tumor')
sc.m$disease = ifelse(data.ph[sc.m$orig.ident,'histology-reduced (who categories based on diagnosis reported in surgical pathology report):ch1']=='Adeno','LUAD','LUSC')
print(unique(sc.m$orig.ident))


sc.m = PercentageFeatureSet(sc.m, "^HB[^(P)]", col.name = "percent_hb")
sc.m = PercentageFeatureSet(sc.m, "^MT-", col.name = "percent_mito")

sc.m = NormalizeData(sc.m)
selected_genes <- rownames(sc.m)[ Matrix::rowSums(sc.m@assays$RNA@counts > 0)> 30]
selected_genes.e <- selected_genes[!grepl('\\.',selected_genes)] #19063

sc.m.2 = FindVariableFeatures(sc.m[selected_genes.e,],nfeatures=5000)

selected_genes.e = VariableFeatures(sc.m.2)
sc.m = sc.m[,sc.m$nFeature_RNA > 250]
sc.m = sc.m[,sc.m$percent_mito < 20]

expr.d = sc.m@assays$RNA@counts[selected_genes.e,]
dim(expr.d)#  5000 51093
gene_annotation = data.frame(id = rownames(expr.d),gene_short_name=rownames(expr.d),num_cells_expresssed = Matrix::rowSums(expr.d!=0))

####3.2 merge as one monocle object####
cds <- new_cell_data_set(as(expr.d, "sparseMatrix"),
                         cell_metadata = sc.m@meta.data,
                         gene_metadata = gene_annotation)

gc()
rm(expr.d)
rm(sc.m.2)
gc()


meta.data.show = cds@colData
meta.data.show = as.data.frame(meta.data.show)
# 
# ggplot2::ggplot(meta.data.show,ggplot2::aes(x=cluster,y=nFeature_RNA,fill=cluster))+
#   ggplot2::geom_violin(scale = 'width')+
#   theme_cowplot()+
#   theme(axis.text.x = ggplot2::element_text(angle=90))
# 
# ggplot2::ggplot(meta.data.show,ggplot2::aes(x=cluster,y=percent_mito,fill=cluster))+
#   ggplot2::geom_violin(scale = 'width')+
#   theme_cowplot()+
#   theme(axis.text.x = ggplot2::element_text(angle=90))
ggplot2::ggplot(meta.data.show,ggplot2::aes(x=orig.ident,y=nFeature_RNA,fill=orig.ident))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')

####3.3 preprocess of cds####

cds[['cellType.pre']] = sc.m$cellType.pre


cds <- preprocess_cds(cds, num_dim = 50,method="PCA",norm_method="log") 
plot_pc_variance_explained(cds)



####3.4 remove batch effects and reduction####
cds <- align_cds(cds, num_dim = 30, alignment_group = "orig.ident")
cds <- reduce_dimension(cds,cores = 20)



plot_cells(cds, color_cells_by="orig.ident",label_cell_groups = F)

plot_cells(cds, color_cells_by="cellType.pre",label_cell_groups = F)



plot_cells(cds,genes = c('CD8A','CD8B','CD4','IL7R','CD3E','CD3D'),cell_size = 0.1,cell_stroke = 0,scale_to_range = F)+
  scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))


plot_cells(cds,genes = c('COL1A1','CD79A','LYZ','EPCAM'),cell_size = 0.1,cell_stroke = 0,scale_to_range = F)+
  ggplot2::scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))

####3.5 clustering####
cds <- cluster_cells(cds, resolution=1e-4)
cols.cl = colorRampPalette(brewer.pal(8, "Accent"))(length(unique(clusters(cds))))

plot_cells(cds,group_label_size = 4,cell_size = 0.1,cell_stroke = 0)+ggplot2::scale_color_manual(values = cols.cl)



table(colData(cds)$orig.ident,clusters(cds))
colData(cds)$cluster = paste0('C',clusters(cds))


####3.6 find markers####
marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group=50,
                               reference_cells=1000, cores=12)
marker_test_res$cell_group = paste('C',marker_test_res$cell_group)

marker_test_res$cell_group = paste(marker_test_res$cell_group,'M')
unique(cds[['From']])


####3.7 annotation####
meta.data.show$cluster = colData(cds)$cluster
ggplot2::ggplot(meta.data.show,ggplot2::aes(x=reorder(cluster,nFeature_RNA,median),y=nFeature_RNA,fill=cluster))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')




colData(cds)$assigned_cell_type <- as.character(clusters(cds))


colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1"="T",
                                                 "2"="B",
                                                 "3"="Myeloid",
                                                 "4"="T",
                                                 '5'='Myeloid',
                                                 '6'='Myeloid',
                                                 "7"="NK",
                                                 "8"="Myeloid",
                                                 "9"="T",#CXCL13 or Treg
                                                 "10"="Myeloid",#
                                                 "11"="Epi",
                                                 "12"="Myeloid",#
                                                 '13'='Epi',
                                                 "14"="Myeloid",#
                                                 "15"="Plasma",
                                                 "16"="Fibro",
                                                 "17"="Mast",#
                                                 "18"="T",#CXCL13
                                                 "19"="Epi",#
                                                 "20"="Plasma",
                                                 "21"="Endo",#
                                                 "22"="Epi",
                                                 "23"="Myeloid",
                                                 "24"="Myeloid",
                                                 "25"="Myeloid",
                                                 "26"="HBB",
                                                 '27'='Epi',
                                                 '28'='Myeloid',
                                                 '29'='Epi'
                                                 
                                                 
)

colData(cds)$majorCellType = unlist(lapply(colData(cds)$assigned_cell_type,function(a){
  strsplit(a,'_')[[1]][1]
}))

cols.cl =cols.default[1:length(unique(colData(cds)$majorCellType))]
names(cols.cl)=c('T','B','Plasma','Myeloid','Epi','Mast','Endo','Fibro','HBB','NK')
plot_cells(cds,  color_cells_by="majorCellType",
           cell_size=0.1,
           cell_stroke=0,
           label_cell_groups = F)+ggplot2::scale_color_manual(values = cols.cl)
plot_cells(cds,  color_cells_by="majorCellType",
           cell_size=0.1,
           cell_stroke=0,
           label_cell_groups = F)+ggplot2::scale_color_manual(values = cols.cl)+
  facet_grid(~From)




colData(cds)$cellType_cluster <- paste0(colData(cds)$majorCellType,'_C',as.character(clusters(cds)))
saveRDS(cds,file = './step3_anno_cds_GSE127465.Rds')



cc.t = table(cds[['cellType_cluster']],cds[['orig.ident']])

pheatmap::pheatmap(scale(cc.t),clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation')

cc.t.cor = cor(t(cc.t),method = 'spearman')

pheatmap::pheatmap(cc.t.cor)

####Validation of CSC_scores####

cds = readRDS(file = './step3_anno_cds_GSE127465.Rds')
geo.data = sc.m
geo.data$majorCellType = cds[['majorCellType']]
geo.data$cluster = cds[['cluster']]

geo.data$N.stage = ifelse(grepl('Nx',geo.data$stage),'Nx',ifelse(grepl('N1',geo.data$stage),'N1','N0'))
cyto.genes = read.csv(file = './tables_v3/fig 2d-luad.csv')
cyto.genes.2 = read.csv(file = './tables_v3/fig 2d-lusc.csv')

cyto.genes.2 = as.character(cyto.genes.2$X)[2:26]

cyto.genes = union(as.character(cyto.genes$X)[1:25],cyto.genes.2)


geo.sc = NormalizeData(geo.data)
geo.sc = ScaleData(geo.sc)

geo.sc = AddModuleScore(geo.sc,features = list(cyto.genes),name = 'CSC_Score')

plot.data =geo.sc@meta.data
plot.data1 = plot.data[plot.data$majorCellType == 'Epi' & plot.data$From == 'lung tumor' &
                         plot.data$N.stage %in% c('N0','N1'),]

plot.data1$N.stage = factor(plot.data1$N.stage,levels = c('N1','N0'))
a=ggplot(plot.data1,aes(x=N.stage,y=CSC_Score1,fill=N.stage))+
  geom_violin(adjust=3,scale = 'width')+
  coord_flip()+
  ggpubr::stat_compare_means(comparisons = list(c(1,2)))+
  theme_cowplot(font_size = 8)+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank(),
        legend.position = 'none')


geo.epi = geo.sc[,geo.sc$majorCellType == 'Epi' & geo.sc$From =='lung tumor' &
                   geo.sc$N.stage %in% c('N1','N0')]
Idents(geo.epi)=geo.epi$N.stage
diff.N1=FindMarkers(geo.epi,ident.1 = 'N1',ident.2 = 'N0',only.pos = T)

intersect(cyto.genes,rownames(diff.N1)[1:50])

b=DotPlot(geo.epi[,geo.epi$N.stage != 'Nx'],
          features = c('AKR1C1','ALDH3A1','FTH1','GAPDH','ALDH1A1','GPX2')
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))
a|b

Idents(geo.epi)=geo.epi$orig.ident
DotPlot(geo.epi[,geo.epi$N.stage != 'Nx'],
        features = c('IGHG1','IGHG4','IGHA1','IGLC2')
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))



####CXCL13+ T####
marker_test_res = read.csv(file = './tables_v2/all_refine_markers.csv')

CXCL13.markers = filter(marker_test_res,marker_test_q_value<0.01 & cell_group == 'CXCL13+ T') %>% top_n(20,specificity)
CXCL13.markers = as.character(CXCL13.markers$gene_id)
ex.markers = c("CXCL13", "HAVCR2", "PDCD1", "TIGIT", "LAG3", "CTLA4", "LAYN", "RBPJ", "VCAM1", "GZMB", "TOX", "MYO7A")

geo.sc = AddModuleScore(geo.sc,features = list(CXCL13.markers),name = 'CXCL13T_Score')
geo.sc = AddModuleScore(geo.sc,features = list(cc=ex.markers),name = 'exScore')
geo.sc$CXCL13_exp = geo.sc@assays$RNA@counts['CXCL13',]

plot.data2 =geo.sc@meta.data
plot.data2 = plot.data2[plot.data2$majorCellType == 'T' & 
                          plot.data2$From =='lung tumor' & 
                          plot.data2$N.stage %in% c('N0','N1'),]
plot.data2=plot.data2[!is.na(plot.data2$N.stage),]
plot.data2$N.stage = factor(plot.data2$N.stage,levels = c('N1','N0'))

a=ggplot(plot.data2,aes(x=N.stage,y=CXCL13T_Score1,fill=N.stage))+
  geom_violin(adjust=3,scale = 'width')+
  ggpubr::stat_compare_means(comparisons = list(c(1,2)))+
  theme_cowplot(font_size = 8)+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank(), 
        panel.grid = element_blank(),
        legend.position = 'none')+
  coord_flip()
a
#ggplot(plot.data2,aes(x=exScore1,y=log2(CXCL13_exp+1),color=From))+geom_point()

geo.t = geo.sc[,geo.sc$majorCellType == 'T' & geo.sc$From == 'lung tumor' & geo.sc$N.stage != 'Nx']
Idents(geo.t)=geo.t$N.stage
diff.t.n1 = FindMarkers(geo.t,ident.1 = 'N1',ident.2 = 'N0',only.pos = T,logfc.threshold = 0.1)
intersect(CXCL13.markers,rownames(diff.t.n1))
# DoHeatmap(geo.t,features = c('CXCL13','CTLA4','PDCD1','TIGIT','TOX2','BATF','TNFRSF18'),disp.max = 0.2,
#           disp.min = -0.2,
#           group.by = 'N.stage')+
#   scale_fill_gradientn(colours = rev(hcl.colors(100,'RdBu')))
b=DotPlot(geo.t,
          features = c('CXCL13','CTLA4','PDCD1','TIGIT','BATF','TNFRSF18')
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))

a | b

####IFI30 SPP1 MP####
IFI30.diff = read.csv(file = './tables_v2/IFI30.diff.csv',row.names = 1)
IFI30.top = filter(IFI30.diff,p_val_adj<0.01 & avg_log2FC >0) %>% top_n(10,avg_log2FC)
IFI30.top = as.character(rownames(IFI30.top))



geo.sc = AddModuleScore(geo.sc,features = list(c(IFI30.top,'SPP1')),name = 'IFI30MP_Score')

plot.data3 = geo.sc@meta.data

plot.data3 = plot.data3[plot.data3$majorCellType == 'Myeloid' & 
                          plot.data3$From == 'lung tumor' &
                          plot.data3$N.stage %in% c('N0','N1'),]
plot.data3$N.stage = factor(plot.data3$N.stage,levels = c('N1','N0'))

a=ggplot(plot.data3,
         aes(x=N.stage,y=IFI30MP_Score1,fill=N.stage))+
  ggpubr::stat_compare_means(comparisons = list(c(1,2)))+
  geom_violin(adjust = 3,scale = 'width')+
  coord_flip()+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = ggsci::pal_npg()(7))+
  stat_summary(fun = mean, geom = "point", size = 0.5,
               color = "black")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none')
a
geo.myeloid = geo.sc[,geo.sc$majorCellType == 'Myeloid' & geo.sc$From == 'lung tumor' &
                       geo.sc$N.stage %in% c('N0','N1'),]

Idents(geo.myeloid)=geo.myeloid$N.stage
diff.m.n1 = FindMarkers(geo.myeloid,ident.1 = 'N1',ident.2 = 'N0',only.pos = T,
                        logfc.threshold = 0.1)
intersect(IFI30.top,rownames(diff.m.n1))
DoHeatmap(geo.myeloid,features = c('SPP1','IFI30','MIF','RNASEK','IGKC'),
          group.by = 'N.stage')+
  scale_fill_gradientn(colours = rev(hcl.colors(100,'RdBu')))
b=DotPlot(geo.myeloid[,geo.myeloid$N.stage !='Nx'],
          features = c('SPP1','IFI30','MIF','RNASEK','IGKC','GABARAP'),
          cols = c('blue','red'),
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))
a|b

save.image(file = 'GSE127475Validation.RData')
