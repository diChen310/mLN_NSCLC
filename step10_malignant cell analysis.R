####Re-analysis of the malignant cells, 2024-09-26####
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

cols.default = c("#93bcd0","#32709f","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559",'#5d94b4',
                 '#b8a3c4','#60428d','#c79795','#c92e2a','#9f5633','#d9c17c',"#d0896b","black",'grey',"blue",'cyan')
####10.1 read epi data####

sc.our = readRDS(file='./variables_v2/sc.our_newName.Rds')
cell_ids = unlist(lapply(sc.our, function(a){unique(a$orig.ident)}))
sc.PLM = merge(sc.our[1][[1]],sc.our[-1],add.cell.ids = cell_ids)

epi.data = readRDS(file = './variables_v2/epi.data.Rds')
epi.data = epi.data[,epi.data[['subCellType']] != 'Remove']

annotations.row = read.csv(file = './tables_v2/all inferCNV annotations of cells.csv',row.names = 1)
malignant.cells.PLM = rownames(annotations.row[annotations.row$CNV.Type == 'Malignant',])#95716
malignant.cells.PLM = intersect(colnames(sc.PLM),malignant.cells.PLM)

expr.d = sc.PLM@assays$RNA@counts[,malignant.cells.PLM]
dim(expr.d)#  54476 95716
gene_annotation = data.frame(id = rownames(expr.d),gene_short_name=rownames(expr.d),num_cells_expresssed = Matrix::rowSums(expr.d!=0))

epi.data.mal = new_cell_data_set(as(expr.d, "sparseMatrix"),
                                 cell_metadata = sc.PLM[,malignant.cells.PLM]@meta.data,
                                 gene_metadata = gene_annotation)
dim(epi.data.mal)
selected_genes <- rownames(epi.data.mal)[gene_annotation$num_cells_expresssed> 500] #19363
epi.data.mal = epi.data.mal[selected_genes,]
dim(epi.data.mal)#19363 95716

epi.data.mal[['PatientID']] = unlist(lapply(epi.data.mal[['orig.ident']],function(a){
  strsplit(a,'\\.')[[1]][1]
}))


# ####10.2 reduction of malignant cells without batch effect correction####
# epi.data.mal = preprocess_cds(epi.data.mal,method="PCA",norm_method="log",num_dim = 20)
# epi.data.mal = reduce_dimension(epi.data.mal,cores = 20)
# epi.data.mal[['PatientID']] = factor(epi.data.mal[['PatientID']],levels = paste0('P',0:18))
# 
# cl.patients = cols.default[1:19]
# names(cl.patients) = paste0('P',0:18)
# 
# scater::plotReducedDim(epi.data.mal,dimred = 'UMAP',colour_by = 'PatientID',point_size = 0.2)+
#   ggplot2::scale_color_manual(values = cl.patients)+
#   theme_void()
# 
# 
# cols.From = as.character(pal_npg()(6))
# names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')
# colData(epi.data.mal)$From[grepl('T',epi.data.mal[['orig.ident']])]='mLN+ PT'
# colData(epi.data.mal)$From[grepl('L',epi.data.mal[['orig.ident']])]='mLN'
# colData(epi.data.mal)$From[grepl('N',epi.data.mal[['orig.ident']])]='mLN+ N'
# 
# colData(epi.data.mal)$From[epi.data.mal[['orig.ident']] == 'P0.3A']='mLN'
# colData(epi.data.mal)$From[epi.data.mal[['orig.ident']] == 'P0.no']='mLN+ PT'
# 
# scater::plotReducedDim(epi.data.mal,dimred = 'UMAP',colour_by = 'From',point_size = 0.2)+
#   ggplot2::scale_color_manual(values = cols.From)+
#   theme_void()


####10.3 reduction of non-malignant cells with batch effect correction####
epi.data.mal = preprocess_cds(epi.data.mal,method="PCA",norm_method="log",num_dim = 20)
epi.data.mal <- align_cds(epi.data.mal, num_dim = 20, alignment_group = "PatientID")
epi.data.mal = reduce_dimension(epi.data.mal,cores = 20)

scater::plotReducedDim(epi.data.mal[,sample(1:ncol(epi.data.mal),3000)],
                       dimred = 'UMAP',colour_by = 'PatientID',point_size = 3)+
  ggplot2::scale_color_manual(values = cl.patients)+
  theme_void()

scater::plotReducedDim(epi.data.mal,dimred = 'UMAP',colour_by = 'From',point_size = 0.2)+
  ggplot2::scale_color_manual(values = cols.From)+
  theme_void()

epi.data.mal = cluster_cells(epi.data.mal,resolution = 5e-5)
epi.data.mal[['cluster.epi.mal']] = paste('MC',clusters(epi.data.mal))
cols.cluster = colorRampPalette(cols.default)(length(unique(epi.data.mal[['cluster.epi.mal']])))
names(cols.cluster) = unique(epi.data.mal[['cluster.epi.mal']])

scater::plotReducedDim(epi.data.mal,dimred = 'UMAP',colour_by = 'cluster.epi.mal',point_size = 0.5)+
  scale_color_manual(values = cols.cluster)+
  theme_void()

table(epi.data.mal[['PatientID']],epi.data.mal[['cluster.epi.mal']])
marker_test_res = top_markers(epi.data.mal, group_cells_by="cluster.epi.mal", genes_to_test_per_group=100,
                              reference_cells=1000, cores=32)
marker_test_res$cell_group = paste(marker_test_res$cell_group,'M')
write.csv(marker_test_res,file = './tables_v2/epi.data.mal.cluster.markers.csv')
epi.data.meta.2 = readRDS(file='./variables_v2/epi.data.meta.2.Rds')

epi.data.mal[['cytoTRACE']] = epi.data.meta.2[colnames(epi.data.mal),'cytoTRACE']
epi.data.mal[['scS']] = epi.data.meta.2[colnames(epi.data.mal),'scS']

saveRDS(epi.data.mal,file = './epi.data.mal.Rds')
epi.data.mal = readRDS(file = './epi.data.mal.Rds')

com2orig = table(epi.data.mal[['orig.ident']],epi.data.mal[['cluster.epi.mal']])
pheatmap(log2(com2orig+1))

meta.data.mal = as.data.frame(colData(epi.data.mal))
meta.data.mal$cluster.epi.mal = factor(meta.data.mal$cluster.epi.mal,levels = paste('MC',1:28))
ggplot(meta.data.mal,aes(x=cluster.epi.mal,fill = PatientID))+
  geom_bar(position = 'fill',width = 0.7)+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.default)+
  theme(axis.text.x = element_text(angle = 90))

ggplot(meta.data.mal,aes(x=cluster.epi.mal,y = cytoTRACE,fill=cluster.epi.mal))+
  geom_boxplot()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))

ggplot(meta.data.mal,aes(x=cluster.epi.mal,y = scS,fill=cluster.epi.mal))+
  geom_boxplot(,outlier.size = 0.1)+
  scale_fill_manual(values = cols.cluster)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))

#####IgG and epi markers relevance#####


sel.markers = c('SFTPB','SFTPA1','SFTPA2','IGHA1','IGHG1','IGHG4','KRT6A','IGLC1','KRT15','KRT5')

plot_genes_by_group(epi.data.mal,sel.markers,
                    group_cells_by="cluster.epi.mal",
                    ordering_type="cluster_row_col",
                    max.size=4)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))+
  theme(axis.text.x = element_text(angle = 90))

epi.data.mal[['IgG_Group']] = ifelse(epi.data.mal[['cluster.epi.mal']] %in%
                                       c('MC21','MC 2','MC 3','MC 9','MC13',"MC 15",'MC 27','MC16'),'IgG+','IgG-')
epi.data.mal[['Disease']] = epi.data.meta.2[colnames(epi.data.mal),'Disease2']

plot_cells(epi.data.mal, genes=c("IGHG4"),cell_size = 1,cell_stroke = 0,scale_to_range = F,
           label_cell_groups = F)+
  scale_color_gradientn(colours = rev(hcl.colors(256,'Spectral')))+
  theme_void()

marker_test_mal.IgG = top_markers(epi.data.mal,group_cells_by = 'IgG_Group', genes_to_test_per_group=100,
                                  reference_cells=1000, cores=32)

meta.data.mal = as.data.frame(colData(epi.data.mal))
meta.data.mal$cluster.epi.mal = factor(meta.data.mal$cluster.epi.mal,levels = paste('MC',1:28))

ggplot(meta.data.mal,aes(x=IgG_Group,y = cytoTRACE,fill=IgG_Group))+
  geom_boxplot(outlier.size = 0.2)+
  ggpubr::stat_compare_means()+
  facet_wrap(~Disease)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

ggplot(meta.data.mal,aes(x=IgG_Group,y = cytoTRACE,fill=IgG_Group))+
  geom_boxplot(outlier.size = 0.2)+
  facet_wrap(Disease~From)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

ggplot(meta.data.mal,aes(x=IgG_Group,fill=PatientID))+
  geom_bar(position = 'fill')+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.default )+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

IgG.markers = filter(marker_test_mal.IgG,cell_group == 'IgG+' & marker_test_p_value <0.01)
path.IgG = enricher(as.data.frame(IgG.markers)$gene_id,TERM2GENE = path2gene.all)@result


go.items = read.gmt(file('./c5.all.v2024.1.Hs.symbols.gmt'))

go.IgG = enricher(as.data.frame(IgG.markers)$gene_id,TERM2GENE = go.items)@result


saveRDS(epi.data.mal,file = './variables_v2/epi.data.mal.Rds')


epi.data.mal.sc = CreateSeuratObject(counts = assay(epi.data.mal))
epi.data.mal.sc = AddMetaData(epi.data.mal.sc,as.data.frame(colData(epi.data.mal)))
Idents(epi.data.mal.sc)=epi.data.mal.sc$IgG_Group
plot_genes_violin(epi.data.mal[c('HLA-DRB1','CST3','HLA-DRA','SLPI'),],group_cells_by = 'IgG_Group',ncol = 4)
plot_genes_violin(epi.data.mal[c('IFNGR1','BCL10','HMOX1','TAP1'),],group_cells_by = 'IgG_Group',ncol = 4,normalize = F)
plot_genes_violin(epi.data.mal[c('CD44','FTH1','FTL','TMSB4X'),],group_cells_by = 'IgG_Group',ncol = 4,normalize = F)
plot_genes_violin(epi.data.mal[c('CD44','FTH1','FTL','TMSB4X'),],group_cells_by = 'IgG_Group',ncol = 4,normalize = F)

saveRDS(epi.data.mal.sc,file = './variables_v2/epi.data.mal.sc.Rds')

Ig.diff.2 = FindMarkers(epi.data.mal.sc,ident.1 = 'IgG+',logfc.threshold = 0,min.pct = 0)
write.csv(Ig.diff.2,file = './tables_v2/Malignant IgGroup diff_seurat version.csv')
Ig.diff.2 = read.csv(file = './tables_v2/Malignant IgGroup diff_seurat version.csv',row.names = 1)
library(clusterProfiler)
go.items = read.gmt(file('./c5.all.v2024.1.Hs.symbols.gmt'))
hallmark.items = read.gmt(file('./MSigDB_Hallmark.gmt'))
c12.log2FC = Ig.diff.2$avg_log2FC
names(c12.log2FC)=rownames(Ig.diff.2)
c12.log2FC.s = sort(c12.log2FC,decreasing = T)
c12.gsea = GSEA(c12.log2FC.s,minGSSize = 5,TERM2GENE = hallmark.items,eps = 0,pvalueCutoff = 2)
c12.gsea.res.hallmark = c12.gsea@result
gseaplot(c12.gsea,geneSetID = 'APOPTOSIS',size=0.5)
ridgeplot(c12.gsea,showCategory = 10)+
  scale_fill_gradientn(colors = rev(hcl.colors(256)))+
  theme_cowplot(font_size = 8)
write.csv(c12.gsea.res.hallmark,file = './tables_v2/GSEA_hallmark_mal IgGroup diff_seurat version.csv')


plot_genes_violin(epi.data.mal[c('CTLA4','CD274','H3F3B','B2M'),],group_cells_by = 'IgG_Group',min_expr = 0.01,
                  ncol = 4,normalize = T)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

plot_genes_violin(epi.data.mal[c('H3F3B','B2M'),],group_cells_by = 'IgG_Group',
                  ncol = 4,normalize = F)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

plot_genes_violin(epi.data.mal[c('FTH1','FTL'),],group_cells_by = 'IgG_Group',ncol = 2,normalize = T)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())



####10.4 difference between mLN and mLN+PT ####
marker_test_LT.diff = top_markers(epi.data.mal, group_cells_by="From", genes_to_test_per_group=100,
                                  reference_cells=1000, cores=32)
mLN.spec = filter(marker_test_LT.diff,cell_group == 'mLN' & marker_test_p_value <0.01) %>% arrange(-specificity)
mLN.spec$label = ''
mLN.spec$label[1:20] = mLN.spec$gene_id[1:20]
ggplot(mLN.spec,aes(x=specificity,y=-log10(marker_test_p_value),color=specificity))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label=label),size=3)+
  scale_color_gradientn(colours = hcl.colors(100))+
  theme_cowplot()
saveRDS(marker_test_LT.diff,file = './tables_v2/marker_test_LT.diff.Rds')


#####LUSC mLN vs mLN+ T#####
epi.data.mal.LUSC = epi.data.mal[,epi.data.mal[['PatientID']] %in% c('P0','P1','P5','P9','P15')]

marker_test_LT.diff.lusc = top_markers(epi.data.mal.LUSC, group_cells_by="From", genes_to_test_per_group=100,
                                       reference_cells=1000, cores=32)
saveRDS(marker_test_LT.diff.lusc,file = './tables_v2/marker_test_LT.diff.lusc.Rds')

mLN.spec.lusc = filter(marker_test_LT.diff.lusc,cell_group == 'mLN' & marker_test_p_value <0.01 & fraction_expressing > 0.8 & substr(gene_id,1,2) %in% c('MT','RP') == F) %>% arrange(-specificity)
mLN.spec.lusc$label = ''
mLN.spec.lusc$label[1:25] = mLN.spec.lusc$gene_id[1:25]


t.spec.lusc = filter(marker_test_LT.diff.lusc,cell_group == 'mLN+ PT' & marker_test_p_value <0.01& fraction_expressing > 0.9 & substr(gene_id,1,2) %in% c('MT','RP') == F) %>% arrange(-specificity)
t.spec.lusc$label = ''
t.spec.lusc$label[1:25] = t.spec.lusc$gene_id[1:25]


merge.spec = data.frame(gene = c(mLN.spec.lusc$gene_id[1:25],
                                 t.spec.lusc$gene_id[1:25]),
                        specificity = c(mLN.spec.lusc$specificity[1:25],
                                        -t.spec.lusc$specificity[1:25]))
merge.spec$gene = factor(merge.spec$gene,levels = rev(merge.spec$gene))
merge.spec$Tumor = 'LUSC'
ggplot(merge.spec,aes(x=specificity,y=gene,fill=specificity))+
  geom_bar(stat = 'identity',width = 0.65)+
  scale_fill_gradientn(colours = rev(hcl.colors(100,'spectral')))+
  theme_cowplot(font_size = 8)

ggplot(merge.spec[1:25,],aes(x=specificity,y=gene,fill=specificity))+
  geom_bar(stat = 'identity',width = 0.65)+
  scale_fill_gradientn(colours = hcl.colors(100))+
  theme_cowplot(font_size = 8)

#####LUAD#####

epi.data.mal.LUAD = epi.data.mal[,epi.data.mal[['PatientID']] %in% c('P0','P1','P5','P9','P15') ==F]

marker_test_LT.diff.luad = top_markers(epi.data.mal.LUAD, group_cells_by="From", genes_to_test_per_group=100,
                                       reference_cells=1000, cores=32)
saveRDS(marker_test_LT.diff.luad,file = './tables_v2/marker_test_LT.diff.luad.Rds')

mLN.spec = filter(marker_test_LT.diff.luad,cell_group == 'mLN' & marker_test_p_value <0.01& fraction_expressing > 0.8 & substr(gene_id,1,2) %in% c('MT','RP') == F) %>% arrange(-specificity)
mLN.spec$label = ''
mLN.spec$label[1:25] = mLN.spec$gene_id[1:25]

t.spec = filter(marker_test_LT.diff.luad,cell_group == 'mLN+ PT' & marker_test_p_value <0.01 & fraction_expressing > 0.75 & substr(gene_id,1,2) %in% c('MT','RP') == F) %>% arrange(-specificity)
t.spec$label = ''
t.spec$label[1:25] = t.spec$gene_id[1:25]


merge.spec.luad = data.frame(gene = c(mLN.spec$gene_id[1:25],
                                      t.spec$gene_id[1:25]),
                             specificity = c(mLN.spec$specificity[1:25],
                                             -t.spec$specificity[1:25]))
merge.spec.luad$gene = factor(merge.spec.luad$gene,levels = rev(merge.spec.luad$gene))
merge.spec.luad$Tumor = 'LUAD'
ggplot(merge.spec.luad,aes(x=specificity,y=gene,fill=specificity))+
  geom_bar(stat = 'identity',width = 0.65)+
  scale_fill_gradientn(colours = rev(hcl.colors(100,'Spectral')))+
  theme_cowplot(font_size = 8)

merge.spec.m = rbind(merge.spec,merge.spec.luad)
write.csv(merge.spec.m,file = './tables_v2/mLN LUAD and LUSC diff marker.csv')


intersect(merge.spec.luad$gene[1:25],merge.spec$gene[1:25])
###"CLDN4" "JUNB"  "DUSP1" "JUN"   "KRT19"

plot_genes_violin(epi.data[c('CLDN4','KRT19',"JUNB" , "DUSP1", "JUN"),],group_cells_by = 'From')+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90))

plot_genes_violin(epi.data[c('CLDN4'),epi.data[['Disease']] %in% c('LUAD','LUAD_Normal')],group_cells_by = 'From')+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90))
plot_genes_violin(epi.data[c('CLDN4'),epi.data[['Disease']] %in% c('LUSC','LUSC_Normal')],group_cells_by = 'From')+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90))


plot_genes_violin(epi.data.mal[c('FTH1'),epi.data.mal[['Disease']] %in% c('LUAD','LUAD_Normal')],group_cells_by = 'From')+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90))
plot_genes_violin(epi.data.mal[c('FTH1'),epi.data.mal[['Disease']] %in% c('LUSC','LUSC_Normal')],group_cells_by = 'From')+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90))

####mBrain####
cols.2 = c(colnames(epi.data.mal),colnames(epi.data[,epi.data[['From']] == 'mLN- PT']))
epi.data.mal.2 = epi.data[,cols.2]

ref.geo = read.delim(file='./variables/GSE131907_Lung_Cancer_cell_annotation.txt')

geo.data = readRDS(file='./variables/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds')

geo.mBrain.cells = ref.geo[ref.geo$Sample_Origin == 'mBrain' & ref.geo$Cell_type.refined =='Epithelial cells','Index']

geo.mBrain.cells.count = geo.data[,intersect(colnames(geo.data),geo.mBrain.cells)]

genes.both = intersect(rownames(epi.data.mal.2),rownames(geo.data))

rownames(ref.geo)=ref.geo$Index

mBrain.epi.meta.data = ref.geo[intersect(colnames(geo.data),geo.mBrain.cells),]
mBrain.epi.meta.data$From = 'mBrain'
mBrain.epi.meta.data$PatientID = mBrain.epi.meta.data$Sample
mBrain.epi.meta.data$orig.ident = mBrain.epi.meta.data$Sample


expr.d = cbind(assay(epi.data.mal.2)[genes.both,],
               geo.mBrain.cells.count[genes.both,])

meta.d = rbind(colData(epi.data.mal.2)[,c('orig.ident','PatientID','From')],
               mBrain.epi.meta.data[,c('orig.ident','PatientID','From')])
gene_annotation = data.frame(id = rownames(expr.d),gene_short_name=rownames(expr.d),num_cells_expresssed = Matrix::rowSums(expr.d!=0))

cds <- new_cell_data_set(as(as.matrix(expr.d), "sparseMatrix"),
                         cell_metadata = meta.d,
                         gene_metadata = gene_annotation)


rm(expr.d)
gc()
cols.From = as.character(pal_npg()(5))
names(cols.From)= c('mBrain','mLN+ PT','mLN','mLN- N','mLN- PT')
saveRDS(cds,file = './variables_v2/cds mBrain.Rds')
plot_genes_violin(cds[c('IGHG4','IGHG1','IGHA1','IGLC2'),],group_cells_by = 'From',ncol=4)+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = cols.From)

