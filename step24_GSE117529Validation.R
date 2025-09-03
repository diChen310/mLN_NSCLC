
library(Seurat)
library(monocle3)

#####ReadData#####

geo.dir.1 = '/share/dichen/data/GSE117529/B16F0/'
geo.dir.2 = '/share/dichen/data/GSE117529/LN6/'
geo.dir.3 = '/share/dichen/data/GSE117529/NoTumor/'

sc.obj.1 = Read10X(geo.dir.1)
sc.obj.1 = CreateSeuratObject(sc.obj.1,min.cells = 3)
sc.obj.1$orig.ident='Parental'


sc.obj.2 = Read10X(geo.dir.2)
sc.obj.2 = CreateSeuratObject(sc.obj.2,min.cells = 3)
sc.obj.2$orig.ident='LN6'

sc.obj.3 = Read10X(geo.dir.3)
sc.obj.3 = CreateSeuratObject(sc.obj.3,min.cells = 3)
sc.obj.3$orig.ident='noTumor'


sc.m = merge(sc.obj.1,c(sc.obj.2,sc.obj.3),
             add.cell.ids = c('Parental','LN6','noTumor'))
sc.m = PercentageFeatureSet(sc.m, "^Hb[^(p)]", col.name = "percent_hb")
sc.m = PercentageFeatureSet(sc.m, "^mt-", col.name = "percent_mito")

Idents(sc.m)=sc.m$orig.ident
P1=VlnPlot(sc.m, features = "nCount_RNA", pt.size = 0,log = F) + NoLegend()
P2=VlnPlot(sc.m, features = "nFeature_RNA", pt.size = 0,log = F) + NoLegend()
P3=VlnPlot(sc.m, features = "percent_mito", pt.size = 0)+ NoLegend()
P4=VlnPlot(sc.m, features = "percent_hb", pt.size = 0)+ NoLegend()

P1 / P2 
P3+ggplot2::geom_hline(yintercept = c(10,15,20),lty=3,color='red')


sc.m = sc.m[,sc.m$nFeature_RNA > 200]
#sc.m = sc.m[,sc.m$percent_mito < 25]
sc.m = sc.m[,sc.m$percent_mito < 10]
dim(sc.m)#18078 28534


####Annotation####

sc.m = NormalizeData(sc.m)
selected_genes <- rownames(sc.m)[ Matrix::rowSums(sc.m@assays$RNA@counts > 0)> 10]

sc.m.2 = FindVariableFeatures(sc.m[selected_genes,],nfeatures=5000)

selected_genes.e = VariableFeatures(sc.m.2)

expr.d = sc.m@assays$RNA@counts[selected_genes.e,]
dim(expr.d)#  5000 28534
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

cds <- preprocess_cds(cds, num_dim = 50,method="PCA",norm_method="log") 
plot_pc_variance_explained(cds)



####3.4 remove batch effects and reduction####
cds <- align_cds(cds, num_dim = 30, alignment_group = "orig.ident")
cds <- reduce_dimension(cds,cores = 20)
plot_cells(cds, color_cells_by="orig.ident",label_cell_groups = F)

plot_cells(cds,genes = c('Cd8a','Cd8b','Cd4','Il7r','Cd3e','Cd3d'),cell_size = 0.1,cell_stroke = 0,scale_to_range = F)+
  scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))


plot_cells(cds,genes = c('Col1a1','Cd79a','Lyz','Epcam'),cell_size = 0.1,cell_stroke = 0,scale_to_range = F)+
  ggplot2::scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))

plot_cells(cds,genes = c('Cxcl13','Ifi30','Spp1','Ighg1'),cell_size = 0.1,cell_stroke = 0,scale_to_range = T)+
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


####3.7 annotation####
meta.data.show$cluster = colData(cds)$cluster
ggplot2::ggplot(meta.data.show,ggplot2::aes(x=reorder(cluster,nFeature_RNA,median),y=nFeature_RNA,fill=cluster))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')




colData(cds)$assigned_cell_type <- as.character(clusters(cds))


colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1"="B",
                                                 "2"="B",
                                                 "3"="T",
                                                 "4"="T",
                                                 '5'='B',
                                                 '6'='T',#Treg
                                                 "7"="T",
                                                 "8"="Epi",
                                                 "9"="T",#
                                                 "10"="Myeloid",#
                                                 "11"="Proliferation",
                                                 "12"="T",#
                                                 '13'='Myeloid',
                                                 "14"="Myeloid",#
                                                 "15"="NK",
                                                 "16"="Myeloid",
                                                 "17"="T",#
                                                 "18"="B",#
                                                 "19"="Myeloid",#
                                                 "20"="Poor",
                                                 "21"="Poor",#
                                                 "22"="Myeloid"
                                                 
                                                 
)

colData(cds)$majorCellType = unlist(lapply(colData(cds)$assigned_cell_type,function(a){
  strsplit(a,'_')[[1]][1]
}))

cols.cl =colorRampPalette(brewer.pal(8, "Accent"))(length(unique(cds[['majorCellType']])))

names(cols.cl)=c('T','B','Myeloid','Epi','NK','Proliferation','Poor')
plot_cells(cds,  color_cells_by="majorCellType",
           cell_size=0.1,
           cell_stroke=0,
           label_cell_groups = F)+ggplot2::scale_color_manual(values = cols.cl)

saveRDS(cds,file = './step24_anno_cds_GSE117529.Rds')


####Validation####
geo.sc = NormalizeData(sc.m)
geo.sc = ScaleData(geo.sc)
geo.sc$majorCellType = cds[['majorCellType']]
rownames(geo.sc)=toupper(rownames(geo.sc))
cyto.genes = read.csv(file = './tables_v3/fig 2d-luad.csv')
cyto.genes.2 = read.csv(file = './tables_v3/fig 2d-lusc.csv')

cyto.genes.2 = as.character(cyto.genes.2$X)[2:26]

cyto.genes = union(as.character(cyto.genes$X)[1:25],cyto.genes.2)


geo.sc = AddModuleScore(geo.sc,features = list(cyto.genes),name = 'CSC_Score')

#####CXCL13#####

marker_test_res = read.csv(file = './tables_v2/all_refine_markers.csv')

####CXCL13+ T####

CXCL13.markers = filter(marker_test_res,marker_test_p_value<0.05 & cell_group == 'CXCL13+ T') %>% top_n(20,specificity)
CXCL13.markers = as.character(CXCL13.markers$gene_id)
CXCL13.markers = c(
  "Btla",      "Cd2" ,      "Cd3d" ,     "Cd3g"   ,   "Cd40lg"  ,  "Ctla4"   ,
  "Cxcl13" , "Gng4"  ,    "Ica1"     , "Klrb1"  ,   "Linc01281" ,"Maf" ,
  "Pdcd1"  ,   "Sirpg","Sla"  ,     "Tigit"   ,  "Tnfrsf18",  "Tox2"   ,
  "Trbc1"  ,   "Trbc2" 
)
ex.markers = c("Cxcl13", "Havcr2", "Pdcd1", "Tigit", "Lag3", "Ctla4", "Layb", "Rbpj", "Vcam1", "Gzmb", "Tox",
               "Myo7a")

geo.sc = AddModuleScore(geo.sc,features = list(CXCL13.markers),name = 'CXCL13T_Score')
geo.sc = AddModuleScore(geo.sc,features = list(cc=ex.markers),name = 'exScore')

plot.data2 =geo.sc@meta.data
plot.data2 = plot.data2[plot.data2$majorCellType == 'T',]

plot.data2$orig.ident = factor(plot.data2$orig.ident,levels = c('LN6','Parental','noTumor'))

a=ggplot(plot.data2,aes(x=orig.ident,y=CXCL13T_Score1,fill=orig.ident))+
  geom_violin(adjust=3,scale = 'width')+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
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

geo.t = geo.sc[,geo.sc$majorCellType == 'T']
Idents(geo.t)=geo.t$orig.ident

b=DotPlot(geo.t,
          features = c('Cxcl13','Ctla4','Pdcd1','Tigit','Batf','Tnfrsf18')
)+scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))

a | b

b
write.csv(b$data,file = './tables_v3/figS10L_1_new.csv')
#####IFI30#####


IFI30.diff = read.csv(file = './tables_v2/IFI30.diff.csv',row.names = 1)
IFI30.top = filter(IFI30.diff,p_val_adj<0.01 & avg_log2FC >0) %>% top_n(10,avg_log2FC)
IFI30.top = as.character(rownames(IFI30.top))
IFI30.top = c("Calhm6",  "Igkc",    "Tnfsf13", "Rnasek",  "Mif",
              "Apoc2" ,  "Mmp12" ,  "Nme2","Gabarap", "Ifi30" )


geo.sc = AddModuleScore(geo.sc,features = list(c(IFI30.top,'Spp1')),name = 'IFI30MP_Score')

plot.data3 = geo.sc@meta.data

plot.data3 = plot.data3[plot.data3$majorCellType == 'Myeloid' ,]
plot.data3$orig.ident = factor(plot.data3$orig.ident,levels = c('LN6','Parental','noTumor'))

a=ggplot(plot.data3,
         aes(x=orig.ident,y=IFI30MP_Score1,fill=orig.ident))+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
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
geo.myeloid = geo.sc[,geo.sc$majorCellType == 'Myeloid',]

Idents(geo.myeloid)=geo.myeloid$orig.ident
diff.m.n1 = FindMarkers(geo.myeloid,ident.1 = 'LN6',ident.2 = 'noTumor',only.pos = T,
                        logfc.threshold = 0.1)
intersect(IFI30.top,rownames(diff.m.n1))
DoHeatmap(geo.myeloid,features = c('Spp1','Mif','Gabarap','Igkc'),
          group.by = 'orig.ident')+
  scale_fill_gradientn(colours = rev(hcl.colors(100,'RdBu')))
b=DotPlot(geo.myeloid,
          features = c('Spp1','Mif','Gabarap','Igkc','Cd274'))+
  scale_size_continuous(range = c(1,5))+
  scale_color_gradientn(colours = rev(hcl.colors(100,'RdBu')))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 8),
        axis.text.y = element_text(size = 8),
        text = element_text(size = 8))
a|b
b
write.csv(b$data,file = './tables_v3/figS10L_2_new.csv')
save.image('GSE117529.RData')

