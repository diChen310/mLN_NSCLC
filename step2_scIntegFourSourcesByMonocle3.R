####Step2, scIntegration by monocle3 pipeline,2024-09-13####
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
checkOverlap = function(genes1,genes2,N){
  
  ab = length(intersect(genes1,genes2))
  ab_= length(genes2)-ab
  a_b = length(genes1) - ab
  other = N-ab-ab_-a_b
  
  table.con<-matrix(c(ab,ab_,a_b,other),nrow = 2)
  res<-fisher.test(table.con)
  p.value<-res$p.value
  return(p.value)
  
}

####default colors####
cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559","#32709f",
                 '#b8a3c4','#60428d','#9f5633','grey','#c92e2a','#c79795','#d9c17c',"#d0896b",'#5d94b4',"black","blue",'cyan')
# cols.default = c('#e69f00','#54b6e9','#009e73','#f0e442','#0072b2','#d55e00',
#                  '#cc79a7','#666666','#ad7700','#1c91d4','#007756','#d5c711',
#                  '#005685', '#a04700','#b14380','#4d4d4d','#ffbe2d','#80c7ef',
#                  '#00f6b3','#f4eb71')
####2.1 read and merge the pre-processed data####
sc.chenlu = readRDS(file = './variables_v2/chenluData.Rds')
sc.chenlu = sc.chenlu[,sc.chenlu$N ==  'N0']
sc.our = readRDS(file='./variables_v2/sc.our_newName.Rds')
sc.ln = readRDS(file = './variables_v2/GSE131907_nLN_newName.Rds')
sc.ana = readRDS(file = './variables_v2/sc.ana_newName.Rds')
sc.ana = sc.ana[,sc.ana$N0 == 'N0']

cell_ids = unlist(lapply(sc.our, function(a){unique(a$orig.ident)}))
#cell_ids[1:11]='' ###The first 11 with orig.ident already
sc.PLM = merge(sc.our[1][[1]],sc.our[-1],add.cell.ids = cell_ids)


sc.chenlu$orig.ident = sc.chenlu$SampleID
Idents(sc.PLM)=sc.PLM$orig.ident
sc.PLM$PatientID = unlist(lapply(sc.PLM$orig.ident,function(a){
  strsplit(a,'\\.')[[1]][1]
}))
sc.ln$PatientID = paste0('GEO',sc.ln$orig.ident)
sc.PLM$Source = 'PLM'
sc.chenlu$Source = 'chenlulab'
sc.ln$Source = 'normalLN'

sc.ana$Source = 'Ana'
sc.ana = sc.ana[,sc.ana$majorCellType != 'Unassigned']
sc.ana$orig.ident = sc.ana$patient_sample
sc.ana$PatientID=sc.ana$patient
####Check the names between different sources####
names.intersect = intersect(intersect(rownames(sc.chenlu),rownames(sc.ln)),intersect(rownames(sc.ana),rownames(sc.PLM)))#14799

sc.m = merge(sc.chenlu[names.intersect,],c(sc.PLM[names.intersect,],sc.ln[names.intersect,],sc.ana[names.intersect,]))
print(unique(sc.m$orig.ident))
rm(sc.chenlu)
rm(sc.PLM)
rm(sc.ln)
rm(sc.ana)
rm(sc.our)
gc()
sc.m@active.assay = 'RNA'

dim(sc.m)# 14799 1051186

#sc.m = NormalizeData(sc.m)
selected_genes <- rownames(sc.m)[ Matrix::rowSums(sc.m@assays$RNA@counts > 0)> 500] #14787


####3.2 QC together####
sc.m = PercentageFeatureSet(sc.m, "^HB[^(P)]", col.name = "percent_hb")
sc.m = PercentageFeatureSet(sc.m, "^MT-", col.name = "percent_mito")
meta.data.show = sc.m@meta.data
meta.data.show = as.data.frame(meta.data.show)
# 
ggplot2::ggplot(meta.data.show,ggplot2::aes(x=Source,y=nFeature_RNA,fill=Source))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')
ggplot2::ggplot(meta.data.show,ggplot2::aes(x=Source,y=nCount_RNA,fill=Source))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')

ggplot2::ggplot(meta.data.show,ggplot2::aes(x=Source,y=percent_mito,fill=Source))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')

sc.m = sc.m[,sc.m$percent_mito < 25]
sc.m = sc.m[,sc.m$nCount_RNA >1000 & sc.m$nFeature_RNA > 300 ]
dim(sc.m)
# 14799 903633
expr.d = sc.m@assays$RNA@counts[selected_genes,]
dim(expr.d)#  14787 903633
gene_annotation = data.frame(id = rownames(expr.d),gene_short_name=rownames(expr.d),num_cells_expresssed = Matrix::rowSums(expr.d!=0))

####3.3 begin analysis by monocle3
cds <- new_cell_data_set(as(expr.d, "sparseMatrix"),
                         cell_metadata = sc.m@meta.data,
                         gene_metadata = gene_annotation)

rm(sc.our)
gc()
rm(expr.d)

gc()


####3.3 preprocess of cds####


colData(cds)$orig.ident[colData(cds)$orig.ident == 'P0.no']='P0.T'#revise tissue type
saveRDS(cds,file = './variables_v2/inter_cds.Rds')
rm(sc.m)
gc()

cds <- preprocess_cds(cds, num_dim = 50,method="PCA",norm_method="log") 
plot_pc_variance_explained(cds)
print(Sys.time())
# "2024-09-13 16:34:43 +08"

####3.4 remove batch effects and reduction####
cds <- align_cds(cds, num_dim = 20, alignment_group = "PatientID")
# "2024-09-13 19:54:56 +08"

cds <- reduce_dimension(cds,cores = 20)
# "2024-09-13 20:12:21 +08"

cols.sources = jcolors('pal8')
names(cols.sources)=unique(cds[['Source']])
plot_cells(cds, color_cells_by="Source", label_cell_groups=FALSE,cell_size = 0.1,cell_stroke = 0)+
  facet_wrap(~Source,ncol = 4)+
  scale_color_manual(values = cols.sources)+
  theme(legend.position = 'none')+
  theme_void()


####3.5 clustering####
cds <- cluster_cells(cds, resolution=1e-6)
print(Sys.time())
# "2024-09-13 20:34:24 +08"1e-6, 27 clusters
cols.cl = colorRampPalette(brewer.pal(12, "Accent"))(length(unique(clusters(cds))))

plot_cells(cds[,sample(colnames(cds),100000)],
           group_label_size = 4,cell_size = 0.1,cell_stroke = 0)+
  ggplot2::scale_color_manual(values = cols.cl)
table(clusters(cds))


colData(cds)$cluster.main = paste('C',clusters(cds))


####3.6 find markers####
marker_test_res <- top_markers(cds, group_cells_by="cluster", genes_to_test_per_group=100,
                               reference_cells=1000, cores=24)
marker_test_res$cell_group = paste0('C',marker_test_res$cell_group)
write.csv(marker_test_res,file = './tables_v2/markers of main clusters.csv')
marker_test_res$cell_group = paste(marker_test_res$cell_group,'M')


####3.7 provide some patient information####

colData(cds)$From = ''
colData(cds)$From[cds[['Source']] == 'Ana']='mLN- PT'
colData(cds)$From[cds[['Source']] == 'normalLN']='nLN'

colData(cds)$From[cds[['Source']] == 'chenlulab' & cds[['Disease']] %in% c('LUAD','LUSC')]='mLN- PT'
colData(cds)$From[cds[['Source']] == 'chenlulab' & cds[['Disease']] %in% c('LUAD_Normal','LUSC_Normal')]='mLN- N'

colData(cds)$From[cds[['Source']] == 'PLM' & grepl('T',cds[['orig.ident']])]='mLN+ PT'
colData(cds)$From[cds[['Source']] == 'PLM' & grepl('L',cds[['orig.ident']])]='mLN'
colData(cds)$From[cds[['Source']] == 'PLM' & grepl('N',cds[['orig.ident']])]='mLN+ N'

colData(cds)$From[cds[['Source']] == 'PLM' & cds[['orig.ident']] == 'P0.3A']='mLN'


colData(cds)$Disease[cds[['patient']] %in% c('Patient 4','Patient 11',"Patient 17","Patient 19")] = 'LUSC'
colData(cds)$Disease[cds[['patient']] %in% c('Patient 4','Patient 11',"Patient 17","Patient 19") == F & cds[['Source']] == 'Ana'] = 'LUAD'


PLM.infor = read.csv(file = './documents/Patient Info.csv')
PLM.infor$Stage2 = gsub('A|B','',PLM.infor$Stage)
GEO.info = read.csv(file = './documents/Geo Info.csv')
rownames(GEO.info)=GEO.info$orig.ident
rownames(PLM.infor)=PLM.infor$Patient_ID
colData(cds)$Disease[colData(cds)$Source == 'PLM'] = PLM.infor[colData(cds)$PatientID[colData(cds)$Source == 'PLM'],'Disease']
colData(cds)$Disease[colData(cds)$orig.ident %in% c('P3.N','P6.N','P7.N')]='LUAD_Normal'
colData(cds)$Disease[colData(cds)$orig.ident %in% c('P5.N')]='LUSC_Normal'
colData(cds)$Disease[colData(cds)$Source == 'normalLN'] = 'LUAD_Normal'

table(cds[['Source']],cds[['From']])

#               mLN mLN- N mLN- PT mLN+ N mLN+ PT    nLN
# Ana            0      0  320005      0       0      0
# chenlulab      0  78451   94125      0       0      0
# normalLN       0      0       0      0       0  37067
# PLM       156076      0       0  33552  184260      0


table(cds[['Source']],cds[['Disease']])
#           LUAD LUAD_Normal   LUSC LUSC_Normal
# Ana       188352           0 131653           0
# chenlulab  34153       41919  59972       36532
# normalLN       0       37067      0           0
# PLM       263432       31549  77001        2003


cols.diseases = jcolors('pal5')[1:4]
names(cols.diseases)=unique(cds[['Disease']])
plot_cells(cds, color_cells_by="Disease", label_cell_groups=FALSE,cell_size = 0.1,cell_stroke = 0)+
  facet_wrap(~Disease,ncol = 4)+
  scale_color_manual(values = cols.diseases)+
  theme(legend.position = 'none')+
  theme_void()

gc()
meta.data.show = cds@colData
meta.data.show = as.data.frame(meta.data.show)

# 
ggplot2::ggplot(meta.data.show,ggplot2::aes(x=cluster.main,y=nFeature_RNA,fill=cluster.main))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')
ggplot2::ggplot(meta.data.show,ggplot2::aes(x=cluster.main,y=scS,fill=cluster.main))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')

ggplot2::ggplot(meta.data.show,ggplot2::aes(x=cluster.main,y=percent_mito,fill=cluster.main))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')

ggplot2::ggplot(meta.data.show,ggplot2::aes(x=orig.ident,y=nFeature_RNA,fill=orig.ident))+
  ggplot2::geom_violin(scale = 'width')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8),
        legend.position = 'none')

ggplot2::ggplot(meta.data.show,ggplot2::aes(x=orig.ident,fill=cluster.main))+
  ggplot2::geom_bar(position = 'fill')+facet_wrap(~Source,ncol = 1,scales = 'free')+theme_cowplot()+
  theme(axis.text.x = ggplot2::element_text(angle=90,size = 8))+
  scale_fill_manual(values = cols.cl)

####3.8 annotation####

colData(cds)$assigned_cell_type <- as.character(clusters(cds))


colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1"="Myeloid",
                                                 "2"="T",
                                                 "3"="T",
                                                 "4"="B",
                                                 '5'='T',
                                                 '6'='Epi',
                                                 "7"="Myeloid",
                                                 "8"="Myeloid",
                                                 "9"="Epi",
                                                 "10"="Plasma",
                                                 "11"="NK",
                                                 "12"="Myeloid",
                                                 '13'='Epi',
                                                 "14"="Myeloid",
                                                 "15"="Epi",
                                                 "16"="Epi",
                                                 "17"="B",
                                                 "18"="Mast",
                                                 "19"="Fibro",
                                                 "20"="Epi",
                                                 "21"="Myeloid",
                                                 "22"="Endo",
                                                 "23"="Epi",
                                                 "24"="Myeloid",
                                                 "25"="Epi",
                                                 "26"="Remove",#high nFeature
                                                 "27"="Remove"#high MT, nFeature
                                                 )

colData(cds)$majorCellType =colData(cds)$assigned_cell_type

cols.cl =cols.default[1:length(unique(colData(cds)$majorCellType))]
names(cols.cl)=c('T','B','Plasma','Epi','Mast','Myeloid','Fibro','Endo','NK','Remove')
plot_cells(cds,  color_cells_by="majorCellType",
           cell_size=0.1,
           cell_stroke=0,
           label_cell_groups = F)+
  ggplot2::scale_color_manual(values = cols.cl)+
  theme_void()

scater::plotReducedDim(cds,dimred = 'UMAP',colour_by = 'majorCellType',point_size = 0.2)+
  ggplot2::scale_color_manual(values = cols.cl)+
  theme_void()
saveRDS(cds,file = './variables_v2/inter_cds.Rds')


cds = readRDS(file='./variables_v2/inter_cds.Rds')

# marker_test_res <- top_markers(cds, group_cells_by="majorCellType", genes_to_test_per_group=100,
#                                reference_cells=1000, cores=36)
# write.csv(marker_test_res,file = './tables_v2/cds_all_markers_cellType.csv')
#marker_test_res.all = read.csv(file = './variables/cds_all_markers.csv')

####3.9 plot and statistics####

cds[['From']] = factor(cds[['From']],levels = c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN'))

cols.From = as.character(jcolors('pal6')[1:6])
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
# 
# cds = cds[,cds[['majorCellType']] != 'Remove']
# 
# cols.cl =cols.default[1:length(unique(colData(cds)$majorCellType))]
# names(cols.cl)=c('T','B','Plasma','Epi','Mast','Myeloid','Fibro','Endo','NK')
# plot_cells(cds,  color_cells_by="majorCellType",
#            cell_size=0.1,
#            cell_stroke=0,
#            label_cell_groups = F)+
#   ggplot2::scale_color_manual(values = cols.cl)+
#   theme_void()
# 
# marker.genes.classical = c('EPCAM','KRT17','KRT19',
#                            'CD79A','MS4A1',"CD19",
#                            'JCHAIN','MZB1','IGKC','IGHM',
#                            'CD3D','CD3E','TRBC2',
#                            'APOE','MARCO','CD68','LYZ',
#                            'COL1A1','COL1A2','ACTA2',
#                            'CLDN5','PECAM1','VWF',
#                            'TPSAB1','TPSB2',
#                            'NKG7','KLRB1','KLRD1'
#                            )
# plot_genes_by_group(cds,
#                     marker.genes.classical,
#                     group_cells_by="majorCellType",
#                     ordering_type="cluster_row_col",
#                     max.size=3)+scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))
# 
# plot_cells(cds,genes = c('EPCAM','LYZ','NKG7','MS4A1','JCHAIN','CD3D','DCN','PECAM1','TPSAB1','MKI67'),
#            cell_size = 0.1,cell_stroke = 0,scale_to_range = F,
#            label_cell_groups = F,
#            rasterize = T)+
#   scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))+
#   theme_void()
# 
# 
# #colData(cds)$cluster = clusters(cds)
# 
# meta.data.show = cds@colData
# meta.data.show = as.data.frame(meta.data.show)
# meta.data.show$Disease2 = substr(meta.data.show$Disease,1,4)
# 
# 
# write.csv(meta.data.show,file = './tables_v2/all cds meta data.csv')
# 
# 
# ggplot2::ggplot(meta.data.show,ggplot2::aes(x=From,fill=majorCellType))+
#   ggplot2::geom_bar(position = 'fill',width = 0.6)+
#   theme_cowplot()+
#   ggplot2::scale_fill_manual(values = cols.cl)+
#   theme(axis.text.x = element_text(angle = 90,hjust = 1))
# 
# 
# 
# 
# meta.data.show$Patient = factor(meta.data.show$Patient,
#                                 levels = paste0('P',0:18))
# ggplot2::ggplot(meta.data.show[!is.na(meta.data.show$Patient),],ggplot2::aes(x=orig.ident,fill=majorCellType))+
#   facet_grid(~Patient,scales = 'free',space = 'free')+
#   ggplot2::geom_bar(position = 'fill',width = 0.6)+
#   theme_cowplot()+
#   ggplot2::scale_fill_manual(values = cols.cl)+
#   theme(axis.text = element_text(size = 8),
#         axis.text.x = element_text(angle = 90,hjust = 1),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 8))
# 
# 
# com2orig = table(cds[['orig.ident']],cds[['cluster']])
# 
# f.data<-as.matrix(com2orig)
# f.data = f.data/apply(f.data, 1, sum)
# #f.data = t(f.data)
# #rownames(ids.all)=ids.all$id
# meta.data.u = meta.data.show
# meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
# rownames(meta.data.u)=meta.data.u$orig.ident
# 
# com2orig = table(meta.data.show$orig.ident,meta.data.show$majorCellType)
# f.data<-as.matrix(com2orig)
# f.data = f.data/apply(f.data, 1, sum)
# com2orig.df = melt(f.data,measure.vars = colnames(f.data))
# com2orig.df$patient = meta.data.u[as.character(com2orig.df$Var1),'PatientID']
# com2orig.df$From = meta.data.u[as.character(com2orig.df$Var1),'From']
# 
# ggplot(com2orig.df,aes(x=From,y=value,fill=From,color=From))+geom_boxplot(outlier.size = 1,outlier.alpha = 0.7)+
#   facet_wrap(~Var2,scales = 'free',ncol = 10)+
#   theme_cowplot()+
#   theme(axis.text = element_text(size=8),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 8),
#         axis.text.x = element_text(size = 0))+
#   scale_fill_manual(values = alpha(brewer.pal(7,'Set1'),0.7))+
#   scale_color_manual(values = brewer.pal(7,'Set1'))
# 
# 
# 
# 
# meta.data.u = meta.data.show
# meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
# rownames(meta.data.u)=meta.data.u$orig.ident
# 
# L.orig.idents = meta.data.u[meta.data.u$From == 'mLN','orig.ident']
# nLN.orig.idents = meta.data.u[meta.data.u$From == 'nLN','orig.ident'] 
# T.orig.idents = meta.data.u[meta.data.u$From == 'mLN+ PT','orig.ident']
# Tumor.orig.idents = meta.data.u[meta.data.u$From == 'mLN- PT','orig.ident']
# Normal.orig.idents = meta.data.u[meta.data.u$From == 'mLN- N','orig.ident']
# 
# N.orig.idents = meta.data.u[meta.data.u$From == 'mLN+ N','orig.ident']
# 
# L2nLN.diff.p = unlist(lapply(1:ncol(f.data), function(a){
#   wilcox.test(as.double(f.data[L.orig.idents,a]),
#               as.double(f.data[nLN.orig.idents,a]))$p.value
# }))
# 
# L2nLN.diff = unlist(lapply(1:ncol(f.data), function(a){
#   mean(as.double(f.data[L.orig.idents,a]))- mean(as.double(f.data[nLN.orig.idents,a]))
# }))
# 
# L2nLN.res = data.frame(cellType = colnames(f.data),pvalue = L2nLN.diff.p,diff = L2nLN.diff)
# 
# T2T.diff.p = unlist(lapply(1:ncol(f.data), function(a){
#   wilcox.test(as.double(f.data[T.orig.idents,a]),
#               as.double(f.data[Tumor.orig.idents,a]))$p.value
# }))
# 
# T2T.diff = unlist(lapply(1:ncol(f.data), function(a){
#   mean(as.double(f.data[T.orig.idents,a]))- mean(as.double(f.data[Tumor.orig.idents,a]))
# }))
# 
# 
# 
# 
# T2T.res = data.frame(cellType = colnames(f.data),pvalue = T2T.diff.p,diff = T2T.diff)
# 
# T2N.diff.p = unlist(lapply(1:ncol(f.data), function(a){
#   wilcox.test(as.double(f.data[T.orig.idents,a]),
#               as.double(f.data[N.orig.idents,a]))$p.value
# }))
# 
# T2N.diff = unlist(lapply(1:ncol(f.data), function(a){
#   mean(as.double(f.data[T.orig.idents,a]))- mean(as.double(f.data[N.orig.idents,a]))
# }))
# 
# 
# T2N.2.diff.p = unlist(lapply(1:ncol(f.data), function(a){
#   wilcox.test(as.double(f.data[Tumor.orig.idents,a]),
#               as.double(f.data[Normal.orig.idents,a]))$p.value
# }))
# 
# T2N.2.diff = unlist(lapply(1:ncol(f.data), function(a){
#   mean(as.double(f.data[Tumor.orig.idents,a]))- mean(as.double(f.data[Normal.orig.idents,a]))
# }))
# 
# mLN2N.diff = unlist(lapply(1:ncol(f.data), function(a){
#   mean(as.double(f.data[L.orig.idents,a]))- mean(as.double(f.data[N.orig.idents,a]))
# }))
# 
# 
# 
# comp.plot = data.frame(cellType = colnames(f.data),FC_mLN2nLN = L2nLN.diff, FC_mLNPT2PT = T2T.diff)
# ggplot(comp.plot,aes(x=FC_mLN2nLN,y=FC_mLNPT2PT))+geom_point(aes(color=cellType))+ggrepel::geom_label_repel(aes(label = cellType,color=cellType))+
#   theme_cowplot()+scale_color_manual(values=cols.cl)+
#   geom_hline(yintercept = 0,lty=3)+
#   geom_vline(xintercept = 0,lty=3)
# write.csv(comp.plot,file = './tables_v2//Fig1_compare 1.csv')
# #saveRDS(cds,file = './variables_v2/inter_cds.Rds')
# 
# #####To note!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!############################
# 
# #The meta.data , some features were not useful
