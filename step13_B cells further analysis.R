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
####1 Read previous results####
cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559","#32709f",
                 '#b8a3c4','#60428d','#9f5633','grey','#c92e2a','#c79795','#d9c17c',"#d0896b",'#5d94b4',"black","blue",'cyan')

B.data = readRDS(file = './variables_v2/B.data.Rds')
colData(B.data)$assigned_cell_type = B.data[['cluster.lymphocyte']] # 

colData(B.data)$assigned_cell_type <- dplyr::recode(colData(B.data)$assigned_cell_type,
                                                    
                                                    'lymphocyte_C4'='Activated B_C1',#CD69
                                                    'lymphocyte_C21'='Naive B_C1',#IgG,FCER2
                                                    'lymphocyte_C51'='Activated B_C2',#NR4A1,CD69,CD83
                                                    'lymphocyte_C43'='GC B_C1',#RGS13, AICDA
                                                    'lymphocyte_C40'='Naive B_C2',#IGHD, SELL
                                                    'lymphocyte_C27'='Memory B_C1',## BANK1, LTB, LY86,BLK
                                                    'lymphocyte_C30'='Plasma_C2',
                                                    'lymphocyte_C10'='Plasma_C1',#
                                                    'lymphocyte_C37'='Plasma_C3',#SSR4,FKBP11,DERL3
                                                    'lymphocyte_C50'='Plasma_C5',
                                                    'lymphocyte_C48'='Plasma_C4',
                                                    'lymphocyte_C17'='Memory B_C2',#BANK1, LTB
                                                    'lymphocyte_C71'='IgG low plasma_C1',
                                                    'lymphocyte_C62'='ZBTB20+ B_C1',#PAX5, EBF1
                                                    'lymphocyte_C22'='Activated B_C3',#CD69,CD83
                                                    'lymphocyte_C47'='Activated B_C4',#CD69,CD83
                                                    'lymphocyte_C35'='Memory B_C3',#BANK1, LTB, LY86,IFI44L,COTL1
                                                    'lymphocyte_C46'='Activated B_C5',#CD69, CD83
                                                    'lymphocyte_C36'='Naive B_C3'
                                                    
)
marker.genes = c('IGHD','IGHM','CCR7','SELL','TCL1A','FCMR','IL4R','FCER2',
                 'IGHA1','IGHA2','IGHG1','IGHG2','IGHG4','IGHGP','IGLC2',
                 'CD27','LY86','BLK','LTB','BANK1','CHMP1B',
                 'CD69','JUN','CD83',
                 'MS4A1','CD79A','CD79B','CD19',
                 'AICDA','RGS13',
                 "ZBTB20",'ANKRD44','BANK1',
                 'JCHAIN','MZB1','XBP','IGHG1','IGKC'
                 
)
plot_genes_by_group(B.data,marker.genes,
                    group_cells_by="subCellType",
                    ordering_type="cluster_row_col",norm_method = 'size_only',
                    max.size=4,scale_max = 10)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))+
  theme(axis.text.x = element_text(angle = 90))
B.data[['subCellType']] = factor(B.data[['subCellType']],levels = c('Naive B','GC B','Memory B','Activated B',
                                                                    'Plasma','IgG low plasma','ZBTB20+ B'))


cols.cell = cols.default[1:length(unique(colData(B.data)$subCellType))]
names(cols.cell)=levels(colData(B.data)$subCellType)
scater::plotReducedDim(B.data,dimred = 'UMAP',colour_by = 'subCellType',point_size = 1,text_by = 'subCellType',text_size = 2)+
  ggplot2::scale_color_manual(values = cols.cell)+
  theme_void()

####2. statistics on the cell sub types####

cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')

meta.data.show = B.data@colData
meta.data.show = as.data.frame(meta.data.show)
meta.data.show$Disease2 = substr(meta.data.show$Disease,1,4)
meta.data.show$subCellType = as.character(meta.data.show$subCellType)


ggplot2::ggplot(meta.data.show,ggplot2::aes(x=From,fill=subCellType))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+facet_wrap(~Disease2,ncol = 1)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values =  cols.cell)+
  theme(axis.text  = element_text(size = 8),
        axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))

com2orig = table(meta.data.show$orig.ident,meta.data.show$subCellType) 

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

ggplot(com2orig.df,aes(x=From,y=value,fill=From))+geom_boxplot(outlier.size = 0.1)+
  facet_wrap(~Var2,scales = 'free',ncol = 11)+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 0))+
  scale_fill_manual(values = cols.From)



####3. From difference####

marker_test_res.B = top_markers(B.data,group_cells_by = 'subCellType',genes_to_test_per_group = 100,reference_cells = 1000,cores = 32)
marker_test_res.B = filter(marker_test_res.B,fraction_expressing > 0.5 & marker_test_p_value < 0.01) %>% group_by(cell_group)

activated.B.markers = filter(marker_test_res.B,marker_test_p_value < 0.01 & 
                               fraction_expressing>0.7 &
                       substr(gene_id,1,2) %in% c('MT','RP') == F & 
                       cell_group == 'Activated B' &
                       gene_id != 'MALAT1') %>% 
  arrange(marker_test_p_value)
activated.B.markers$label = ''
activated.B.markers$label[1:10]=activated.B.markers$gene_id[1:10]
ggplot(activated.B.markers,aes(color=log2(mean_expression+1),y=-log10(marker_test_p_value),x=specificity))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label = label),size=2,max.overlaps = 1000)+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours = hcl.colors(100))



B.data.From = top_markers(B.data[,!grepl('lasma',B.data[['subCellType']])],group_cells_by = 'From',genes_to_test_per_group = 100,reference_cells = 1000,cores = 32)

plasma.data.From = top_markers(B.data[,grepl('lasma',B.data[['subCellType']])],group_cells_by = 'From',genes_to_test_per_group = 100,reference_cells = 1000,cores = 32)

	
mLNPT.markers = filter(plasma.data.From,cell_group=='mLN+ PT' & marker_test_p_value <0.01)
mLN.markers = filter(plasma.data.From,cell_group=='mLN' & marker_test_p_value <0.01)
intersect(mLNPT.markers$gene_id,
          mLN.markers$gene_id)
# [1] "B2M"     "CD74"    "FOS"     "FTH1"    "FTL"     "HLA-B"   "HSP90B1" "HSPA1A"  "HSPA5"   "IGHA1"   "IGHG1"   "IGHG3"   "IGHG4"   "IGHGP"   "JUN"     "MT-ATP8" "MT-CO1"  "MT-CO2"  "MT-CO3"  "MT-CYB" 
# [21] "MT-ND3"  "MT-ND4"  "MT-ND4L" "MT-ND5"  "PTMA"    "RGS1"    "RPL13"   "RPL23A"  "RPL30"   "RPL37"   "RPL37A"  "RPL38"   "RPLP1"   "RPS11"   "RPS16"   "RPS19"   "RPS27"   "RPS29"   "SF3B1"   "TPT1"   
# [41] "UBC"    

mLNPT.markers = filter(B.data.From,cell_group=='mLN+ PT' & marker_test_p_value <0.01)
mLN.markers = filter(B.data.From,cell_group=='mLN' & marker_test_p_value <0.01)
intersect(mLNPT.markers$gene_id,
          mLN.markers$gene_id)
# [1] "ACTB"      "AMBRA1"    "ARHGDIB"   "B2M"       "CD52"      "CD53"      "CD69"      "CD74"      "DDX5"      "DUSP1"     "FOS"       "FTH1"      "HNRNPA2B1" "HSP90AB1"  "HSPA8"     "IGHG4"     "IGHM"     
# [18] "JUN"       "KLF6"      "MALAT1"    "MS4A1"     "MT-ATP6"   "MT-ATP8"   "MT-CO1"    "MT-CO2"    "MT-CO3"    "MT-CYB"    "MT-ND1"    "MT-ND2"    "MT-ND4"    "MT-ND5"    "PABPC1"    "PFN1"      "PTPRC"    
# [35] "RGS1"      "RPL13"     "RPL15"     "RPL30"     "RPL31"     "RPLP0"     "RPS19"     "RPS2"      "SF3B1"     "TMSB4X"    "TXNRD1"    "UBB"       "ZFP36L1"  

plot_genes_violin(B.data[,!grepl('lasma',B.data[['subCellType']])][c('CD69','JUN','ALDOA','FOS',
                              'IGHG1','IGHG4','KIAA1551','RGS1'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot_genes_violin(B.data[,grepl('Plasma',B.data[['subCellType']])][c('RGS1','IFI6','ALDOA','IGHV4-34',
                                                                     'IGKV1-5','IGHG4','IGHG3','FKBP11'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))


plot_genes_violin(B.data[,grepl('Plasma',B.data[['subCellType']]) & B.data[['From']] %in% c('mLN','mLN+ PT','mLN- PT')][c('PDCD1','LAG3','HAVCR2'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))


plot_genes_violin(B.data[c('PDCD1','LAG3','HAVCR2'),],
                  group_cells_by = 'subCellType',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot_genes_violin(B.data[c('ZBTB20','AFF3','BANK1','RALGPS2'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plasma.data = B.data[,B.data[['subCellType']] == 'Plasma']

plot_genes_violin(plasma.data[c('PDCD1','HAVCR2','LAG3'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  stat_compare_means(comparisons = list(c(2,5),c(3,5)))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot.data = plot_genes_violin(plasma.data[c('PDCD1','HAVCR2','LAG3'),],
                              group_cells_by = 'From')

plasma.data[['PDCD1']] = plot.data$data$expression[plot.data$data$id == 'PDCD1']
plasma.data[['HAVCR2']] = plot.data$data$expression[plot.data$data$id == 'HAVCR2']
plasma.data[['LAG3']] = plot.data$data$expression[plot.data$data$id == 'LAG3']

plasma.data[['ExhaustionScore']] = plasma.data[['PDCD1']]+plasma.data[['HAVCR2']]+plasma.data[['LAG3']]

meta.data.plasma = as.data.frame(colData(plasma.data))

ggplot(meta.data.plasma,aes(x=From,y=log2(PDCD1+1)))+geom_violin(scale = 'width')+
  stat_compare_means(comparisons = list(c(2,5),c(3,5)))

####4. CytoTRACE analysis####
library(CytoTRACE)
B.data.LUSC = B.data[,B.data[['Disease']] %in% c('LUSC','LUSC_Normal')]
expr.LUSC = as.matrix(B.data.LUSC@assays@data$counts)

# B.data.LUSC = preprocess_cds(B.data.LUSC)
# B.data.LUSC = reduce_dimension(B.data.LUSC,cores = 20)
# 
selected_genes <- rownames(expr.LUSC)[ Matrix::rowSums(B.data.LUSC@assays@data$counts > 0)> 100]
cyto.LUSC <- CytoTRACE(expr.LUSC[selected_genes,],batch = as.character(colData(B.data.LUSC)$PatientID),ncores = 12)
#plotCytoGenes(cyto.LUSC, numOfGenes = 10)

gc()

cyto.genes = cyto.LUSC$cytoGenes


sub.cells = sample(colnames(B.data.LUSC),500)
cyto.genes.temp = cyto.genes[substr(names(cyto.genes),1,2) %in% c('MT','RP')==F]
top.genes = c(cyto.genes.temp[1:25],rev(cyto.genes.temp)[1:25])
plot.mat = log2( as.data.frame(t(expr.LUSC[names(top.genes),sub.cells]))+1)

plot.mat$CytoTRACE = as.numeric(cyto.LUSC[['CytoTRACE']][sub.cells])
plot.mat = plot.mat[order(plot.mat$CytoTRACE),]

pheatmap::pheatmap(t(as.data.frame(lapply(plot.mat,minMax),check.names = F))[c('CytoTRACE',names(top.genes)),],color = rev(hcl.colors(100,'Spectral')),show_colnames = F,cluster_cols = F,
                   cluster_rows = F,gaps_row = c(1,26),fontsize = 8)


B.data.LUSC[['CytoTRACE']] = cyto.LUSC[['CytoTRACE']]



meta.data.B.LUSC = as.data.frame(colData(B.data.LUSC))
ggplot(meta.data.B.LUSC,aes(y=CytoTRACE,x=From,fill=From,color=From))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+scale_fill_manual(values = alpha(cols.From,0.7))+
  scale_color_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle=90))

#####LUAD#####

B.data.LUAD = B.data[,B.data[['Disease']] %in% c('LUAD','LUAD_Normal')]
expr.LUAD = as.matrix(B.data.LUAD@assays@data$counts)

# B.data.LUAD = preprocess_cds(B.data.LUAD)
# B.data.LUAD = reduce_dimension(B.data.LUAD,cores = 20)
# 
selected_genes <- rownames(expr.LUAD)[ Matrix::rowSums(B.data.LUAD@assays@data$counts > 0)> 100]
cyto.LUAD <- CytoTRACE(expr.LUAD[selected_genes,],batch = as.character(colData(B.data.LUAD)$PatientID),ncores = 12)
#plotCytoGenes(cyto.LUAD, numOfGenes = 10)

gc()

cyto.genes = cyto.LUAD$cytoGenes


sub.cells = sample(colnames(B.data.LUAD),500)
cyto.genes.temp = cyto.genes[substr(names(cyto.genes),1,2) %in% c('MT','RP')==F]
top.genes = c(cyto.genes.temp[1:25],rev(cyto.genes.temp)[1:25])
plot.mat = log2( as.data.frame(t(expr.LUAD[names(top.genes),sub.cells]))+1)

plot.mat$CytoTRACE = as.numeric(cyto.LUAD[['CytoTRACE']][sub.cells])
plot.mat = plot.mat[order(plot.mat$CytoTRACE),]

pheatmap::pheatmap(t(as.data.frame(lapply(plot.mat,minMax),check.names = F))[c('CytoTRACE',names(top.genes)),],color = rev(hcl.colors(100,'Spectral')),show_colnames = F,cluster_cols = F,
                   cluster_rows = F,gaps_row = c(1,26),fontsize = 8)


B.data.LUAD[['CytoTRACE']] = cyto.LUAD[['CytoTRACE']]



meta.data.B.LUAD = as.data.frame(colData(B.data.LUAD))
ggplot(meta.data.B.LUAD,aes(y=CytoTRACE,x=From,fill=From,color=From))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+scale_fill_manual(values = alpha(cols.From,0.7))+
  scale_color_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle=90))

meta.data.show$CytoTRACE = 0
meta.data.show[colnames(B.data.LUAD),'CytoTRACE'] =  cyto.LUAD[['CytoTRACE']]
meta.data.show[colnames(B.data.LUSC),'CytoTRACE'] =  cyto.LUSC[['CytoTRACE']]


ggplot(meta.data.show,aes(y=CytoTRACE,x=subCellType,fill=subCellType))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.cell)+
  facet_wrap(~Disease2)+
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_blank())

saveRDS(meta.data.show,file = './variables_v2/meta.data.show.B.Rds')

####5. Plasma cell clustering####


plasma.data <- monocle3::preprocess_cds(plasma.data, num_dim = 50,method="PCA",norm_method="log") 
plasma.data <- monocle3::align_cds(plasma.data, num_dim = 20, alignment_group = "PatientID")
plasma.data <- monocle3::reduce_dimension(plasma.data,cores = 20)
plasma.data = monocle3::cluster_cells(plasma.data, resolution=5e-5)

cols.cl = colorRampPalette(cols.default)(length(unique(monocle3::clusters(plasma.data)))) 
monocle3::plot_cells(plasma.data,group_label_size = 4,cell_size = 1,cell_stroke = 0)+ggplot2::scale_color_manual(values = cols.cl)


plasma.data[['cluster.plasma']] = paste0('Plasma_C',clusters(plasma.data))
marker_test_res.plasma = top_markers(plasma.data,group_cells_by = 'cluster.plasma',genes_to_test_per_group = 100,reference_cells = 1000,cores = 32)

scater::plotReducedDim(plasma.data,dimred = 'UMAP',colour_by = 'cluster.plasma',point_size = 1)+
  ggplot2::scale_color_manual(values = cols.cl)+
  theme_void()

cc=table(plasma.data[['From']],plasma.data[['cluster.plasma']])/apply(table(plasma.data[['From']],plasma.data[['cluster.plasma']]),1,sum)
pheatmap::pheatmap(cc,color = rev(hcl.colors(100,'Spectral')))

marker_test_res.plasma.c1 = filter(marker_test_res.plasma,cell_group == 'Plasma_C1' & marker_test_p_value < 0.01 ) %>% top_n(30,specificity)
marker_test_res.plasma.c2 = filter(marker_test_res.plasma,cell_group == 'Plasma_C2' & marker_test_p_value < 0.01 ) %>% top_n(30,specificity)

plasma.sc = CreateSeuratObject(counts = assay(plasma.data))
plasma.sc = AddMetaData(plasma.sc,as.data.frame(colData(plasma.data)))
Idents(plasma.sc)=plasma.sc$cluster.plasma

plasma.c1c2.diff = FindMarkers(plasma.sc,ident.1 = 'Plasma_C2',ident.2 = 'Plasma_C1')
write.csv(plasma.c1c2.diff,file='./plasma C1C2 diff.csv')
plasma.c1c2.diff = read.csv(file='./plasma C1C2 diff.csv')
# Idents(plasma.sc)=plasma.sc$From
# 
# plasma.mLND.diff = FindMarkers(plasma.sc,ident.1 = 'mLN+ PT',ident.2 = 'mLN- PT')

plasma.c1c2.diff.sig = plasma.c1c2.diff[plasma.c1c2.diff$p_val_adj<0.01,]
plasma.c1c2.diff.sig = plasma.c1c2.diff.sig[order(plasma.c1c2.diff.sig$avg_log2FC),]

plasma.c1c2.diff$label = ''

plasma.c1c2.diff$label[rownames(plasma.c1c2.diff) %in% rownames(plasma.c1c2.diff.sig)[c(1:10,1166:1175)]] = 
  rownames(plasma.c1c2.diff) [rownames(plasma.c1c2.diff) %in% rownames(plasma.c1c2.diff.sig)[c(1:10,1166:1175)]] 
plasma.c1c2.diff$p_val_adj_adj = plasma.c1c2.diff$p_val_adj
plasma.c1c2.diff$p_val_adj_adj[plasma.c1c2.diff$p_val_adj<1e-160] = 1e-160
ggplot(plasma.c1c2.diff,aes(x=avg_log2FC,y=-log10(p_val_adj_adj),color=avg_log2FC))+geom_point()+
  ggrepel::geom_text_repel(aes(label = label),size=2)+
  scale_color_gradientn(colors = rev(hcl.colors(100,'Spectral')))+
  theme_cowplot(font_size = 8)


top.plasma.c1c2.diff = plasma.c1c2.diff.sig[c(1:10,1166:1175),]

top.plasma.c1c2.diff$ID = factor(rownames(top.plasma.c1c2.diff),levels = rownames(top.plasma.c1c2.diff))

ggplot(top.plasma.c1c2.diff,aes(x=avg_log2FC,y=ID,fill = avg_log2FC))+geom_bar(stat = 'identity',width = 0.7)+
  theme_cowplot(font_size = 8)+
  scale_fill_gradientn(colours = rev(hcl.colors(100,'Spectral')))


plasma.data.c1c2 = plasma.data[,plasma.data[['cluster.plasma']] %in% c("Plasma_C1",'Plasma_C2')]
plot_genes_violin(plasma.data.c1c2[c('IGHG2','IGLV2-8','IGLL5','JCHAIN',
                                     'CEP128',"IGHV4-34", 'XIST','OGT',"IFNG-AS1"),],
                  group_cells_by = 'cluster.plasma',
                  ncol = 5,normalize = F)+
  scale_fill_manual(values = cols.cl)+
  ggpubr::stat_compare_means()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,size = 8))

plot_genes_violin(plasma.data.c1c2[c('IGHG1','IGHG3','IGHG4','IGKC',
                                     'IGHM',"IGHA1", 'IGHA2','IGHG2'),],
                  group_cells_by = 'cluster.plasma',
                  ncol = 8,normalize = F)+
  ggpubr::stat_compare_means(aes(label = paste0("p = ", ..p.format..)))+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.cl)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,size = 8),
        strip.background = element_blank())
saveRDS(plasma.data,file = './plasma.data.Rds')
plasma.data = readRDS(file = './plasma.data.Rds')

library(clusterProfiler)
go.items = read.gmt(file('./c5.all.v2024.1.Hs.symbols.gmt'))
hallmark.items = read.gmt(file('./MSigDB_Hallmark.gmt'))
c12.log2FC = plasma.c1c2.diff$avg_log2FC
names(c12.log2FC)=rownames(plasma.c1c2.diff)
c12.log2FC.s = sort(c12.log2FC,decreasing = T)
c12.gsea = GSEA(c12.log2FC.s,minGSSize = 5,TERM2GENE = go.items,eps = 0)
c12.gsea.res = c12.gsea@result
ridgeplot(c12.gsea)
c12.gsea.res.gobp = c12.gsea.res[grepl('GOBP',c12.gsea.res$ID),] 
c12.gsea.res.gobp$Type = ifelse(c12.gsea.res.gobp$NES>0,'Up','Down')
c12.gsea.res.gobp.top = filter(c12.gsea.res.gobp,p.adjust<0.05) %>% group_by(Type) %>% top_n(10,abs(NES))
c12.gsea.res.gobp.top$ID = gsub('GOBP_','',c12.gsea.res.gobp.top$ID)
c12.gsea.res.gobp.top$ID =  factor(c12.gsea.res.gobp.top$ID,levels = c12.gsea.res.gobp.top$ID)
ggplot(c12.gsea.res.gobp.top,aes(x=NES,y=ID,fill=-log10(p.adjust)))+
  geom_bar(stat = 'identity',width = 0.65)+
  scale_fill_gradientn(colours = rev(hcl.colors(100)))+
  theme_cowplot(font_size = 8)

###Pseudotime analysis


