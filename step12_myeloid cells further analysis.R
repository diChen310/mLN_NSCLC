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

Myeloid.data = readRDS(file = './variables_v2/myeloid.data.refine.Rds')

marker.genes = c('FCER1A','CD1C','CLEC10A','CST3','HLA-DRA','CCL17',
                                'XCR1','CLEC9A',
                                'S100A8','S100A9','G0S2','CXCL8',
                                'S100A12','APOBEC3A','FCN1',
                                'LAMP3','HLA-DPA1',
                                'IRF7','IRF8','TCF4','GZMB',
                                'FABP4',
                                'LYZ','MARCO','CD68','APOE','CD14',
                 "GPNMB",'ITGAM',
                 
                                'SPP1','C1QB','CD163',
                                 "TUBB",'STMN1',
                                'LST1','AIF1','FCGR3A'
)


plot_genes_by_group(Myeloid.data,marker.genes,
                    group_cells_by="subCellType",
                    ordering_type="cluster_row_col",
                    max.size=4)+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))+
  theme(axis.text.x = element_text(angle = 90))
Myeloid.data[['subCellType']] = factor(Myeloid.data[['subCellType']],levels = c('cDC1','cDC2','pDC','LAMP3+ DC','CD1C- DC',
                                                                                'APOE+ MP','FABP4+ MP','SPP1+ MP','Gran','Neu',
                                                                                "Proliferation Myeloid"))


cols.cell = cols.default[1:length(unique(colData(Myeloid.data)$subCellType))]
names(cols.cell)=levels(colData(Myeloid.data)$subCellType)
scater::plotReducedDim(Myeloid.data,dimred = 'UMAP',colour_by = 'subCellType',
                       point_size = 1)+
  ggplot2::scale_color_manual(values = cols.cell)+
  theme_void()

Myeloid.data[['subCellType_level1']] = recode(Myeloid.data[['subCellType']] ,
                                              'cDC1'='DC',
                                              'cDC2'='DC',
                                              'pDC'='DC',
                                              'LAMP3+ DC'='DC',
                                              'CD1C- DC'='DC',
                                              'APOE+ MP'='MP',
                                              'SPP1+ MP'='MP',
                                              'FABP4+ MP'='MP') 
scater::plotReducedDim(Myeloid.data,dimred = 'UMAP',colour_by = 'subCellType_level1',point_size = 1,text_by = 'subCellType_level1')+
  theme_void()

####2. statistics on the cell sub types####
cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')

meta.data.show = Myeloid.data@colData
meta.data.show = as.data.frame(meta.data.show)
meta.data.show$Disease2 = substr(meta.data.show$Disease,1,4)
meta.data.show$subCellType = as.character(meta.data.show$subCellType)


ggplot2::ggplot(meta.data.show,ggplot2::aes(x=From,fill=subCellType_level1))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+facet_wrap(~Disease2,ncol = 1)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values =  alpha(pal_aaas()(6),0.7))+
  theme(axis.text  = element_text(size = 8),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))

ggplot2::ggplot(meta.data.show[meta.data.show$subCellType_level1 == 'DC',],ggplot2::aes(x=From,fill=subCellType))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+facet_wrap(~Disease2,ncol = 1)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values =  cols.cell)+
  theme(axis.text  = element_text(size = 8),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))

ggplot2::ggplot(meta.data.show[meta.data.show$subCellType_level1 == 'MP',],ggplot2::aes(x=From,fill=subCellType))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+facet_wrap(~Disease2,ncol = 1)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values =  cols.cell)+
  theme(axis.text  = element_text(size = 8),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))

com2orig = table(meta.data.show$orig.ident,meta.data.show$subCellType)

f.data<-as.matrix(com2orig)
f.data = f.data/apply(f.data, 1, sum)
meta.data.u = meta.data.show
meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
rownames(meta.data.u)=meta.data.u$orig.ident

com2orig.df = melt(f.data,measure.vars = colnames(f.data))
com2orig.df$patient = meta.data.u[as.character(com2orig.df$Var1),'PatientID']
com2orig.df$From = meta.data.u[as.character(com2orig.df$Var1),'From']
com2orig.df$Disease = meta.data.u[as.character(com2orig.df$Var1),'Disease2']



L.orig.idents = meta.data.u[meta.data.u$From == 'mLN','orig.ident']
nLN.orig.idents = meta.data.u[meta.data.u$From == 'nLN','orig.ident'] 
T.orig.idents = meta.data.u[meta.data.u$From == 'mLN+ PT','orig.ident']
Tumor.orig.idents = meta.data.u[meta.data.u$From == 'mLN- PT','orig.ident']
Tumor.orig.idents = intersect(Tumor.orig.idents,rownames(f.data))

L2nLN.diff.p = unlist(lapply(1:ncol(f.data), function(a){
  wilcox.test(as.double(f.data[L.orig.idents,a]),
              as.double(f.data[nLN.orig.idents,a]))$p.value
}))

L2nLN.diff = unlist(lapply(1:ncol(f.data), function(a){
  mean(as.double(f.data[L.orig.idents,a]))- mean(as.double(f.data[nLN.orig.idents,a]))
}))

L2nLN.res = data.frame(cellType = colnames(f.data),pvalue = L2nLN.diff.p,diff = L2nLN.diff)

T2T.diff.p = unlist(lapply(1:ncol(f.data), function(a){
  wilcox.test(as.double(f.data[T.orig.idents,a]),
              as.double(f.data[Tumor.orig.idents,a]))$p.value
}))

T2T.diff = unlist(lapply(1:ncol(f.data), function(a){
  mean(as.double(f.data[T.orig.idents,a]))- mean(as.double(f.data[Tumor.orig.idents,a]))
}))
T2T.res = data.frame(cellType = colnames(f.data),pvalue = T2T.diff.p,diff = T2T.diff)

comp.plot = data.frame(cellType = colnames(f.data),FC_mLN2nLN = L2nLN.diff, FC_mLNPT2PT = T2T.diff)
ggplot(comp.plot,aes(x=FC_mLN2nLN,y=FC_mLNPT2PT))+geom_point(aes(color=cellType))+ggrepel::geom_text_repel(aes(label = cellType,color=cellType),size=2)+
  theme_cowplot(font_size = 6)+scale_color_manual(values=cols.cell)+
  geom_hline(yintercept = 0,lty=3)+
  geom_vline(xintercept = 0,lty=3)
T2T.res$From = 'mLN+ PT vs mLN- PT'
L2nLN.res$From = 'mLN vs nLN'
write.csv(rbind(T2T.res,L2nLN.res),file = './tables_v2/myeloid two comparisions.csv')

####3. Subtype based differences####\

markers_subType = top_markers(Myeloid.data,group_cells_by = 'subCellType',genes_to_test_per_group = 100,reference_cells = 1000,cores = 12)

spp1.data = Myeloid.data[,Myeloid.data[['subCellType']] == 'SPP1+ MP']
DC2.data = Myeloid.data[,Myeloid.data[['subCellType']] == 'cDC2']

scater::plotReducedDim(spp1.data,dimred = 'UMAP',colour_by = 'From',point_size = 1)+
  theme_void()

plot_cells(DC2.data,color_cells_by = 'From',cell_size = 0.1,cell_stroke = 0)+ggplot2::scale_color_manual(values = cols.From)+facet_grid(~From)+
  theme_void()


spp1.diff.from = top_markers(spp1.data,group_cells_by = 'From',genes_to_test_per_group = 100,cores = 12)
write.csv(spp1.diff.from,file = './tables_v2/spp1  from diff res.csv')
mLNn.genes = filter(spp1.diff.from,marker_test_p_value < 0.01 & cell_group == 'mLN+ N')
spp1.diff.filter = spp1.diff.from[spp1.diff.from$gene_id %in% mLNn.genes$gene_id == F,]

spp1.diff.from = read.csv(file = './tables_v2/spp1  from diff res.csv',row.names = 1)
spp1.sc = CreateSeuratObject(counts = assay(spp1.data))
spp1.sc = AddMetaData(spp1.sc,as.data.frame(colData(spp1.data)))
Idents(spp1.sc) = spp1.sc$From
spp1.sc.diff = FindMarkers(spp1.sc,ident.1 = 'mLN+ PT',ident.2 = 'mLN- PT')

plot_genes_violin(spp1.data[c('XIST','H3F3A','TYMP','CCL13',
                             'CCL4','IGHG4','MMP12','IFI30'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))
spp1.data[['SourceFrom']] = paste(spp1.data[['Source']],spp1.data[['From']])


plot_genes_violin(spp1.data[c('IFI30','CXCL9'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot_genes_violin(spp1.data[c('GABARAP','MMP12','MIF','APOC2',
                              'CCL13','XIST','CHIT1','F13A1'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot_genes_violin(spp1.data[c('FCGR2A','FCGR2B','FCGR3A','FCGR3B'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))


DC2.diff.from = top_markers(DC2.data,group_cells_by = 'From',genes_to_test_per_group = 100,cores = 12)

DC2.diff.from.top = filter(DC2.diff.from, fraction_expressing > 0.6 & marker_test_p_value < 0.001) %>% top_n(20,specificity)

dc2.sc = CreateSeuratObject(counts = assay(DC2.data))
dc2.sc = AddMetaData(dc2.sc,as.data.frame(colData(DC2.data)))
Idents(dc2.sc) = dc2.sc$From
dc2.sc.diff = FindMarkers(dc2.sc,ident.1 = 'mLN+ PT',ident.2 = 'mLN- PT')
dc2.sc.diff.2 = FindMarkers(dc2.sc,ident.1 = 'mLN',ident.2 = 'nLN')



plot_genes_violin(DC2.data[c('FCGBP','CCL17','H3F3A','TXNIP',
                             'IGHG4','CCL3','CCL4L2','IFI30'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))
DC2.data[['SourceFrom']] = paste(DC2.data[['Source']],DC2.data[['From']])
plot_genes_violin(DC2.data[c('FCGBP','CCL17','H3F3A','TXNIP'),],
                  group_cells_by = 'SourceFrom',
                  ncol = 4)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot_genes_violin(Myeloid.data[c('IFI30'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

FR.genes = c('FCGR2A','FCGR2B','FCGR3A','FCGR3B')

plot_genes_by_group(spp1.data,FR.genes,
                    group_cells_by="From",
                    ordering_type="cluster_row_col",
                    max.size=4)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))+
  theme(axis.text.x = element_text(angle = 90))
#####Pathway analysis
library(clusterProfiler)
path2gene.all = read.gmt(gmtfile = './MSigDB_Hallmark.gmt')

spp1.markers = filter(markers_subType,cell_group == 'SPP1+ MP' & marker_test_p_value < 0.01)$gene_id
spp1.pathway = enricher(spp1.markers,pvalueCutoff = 1.2,TERM2GENE = path2gene.all)@result


####3. sub-clustering of SPP1-MP####

spp1.data = monocle3::preprocess_cds(spp1.data, num_dim = 50,method="PCA",norm_method="log") 
spp1.data <- monocle3::align_cds(spp1.data, num_dim = 20, alignment_group = "PatientID")
spp1.data <- monocle3::reduce_dimension(spp1.data,cores = 20)
spp1.data = monocle3::cluster_cells(spp1.data, resolution=1e-4)
# 1e-5, cluster number 24
cols.cl = colorRampPalette(cols.default)(length(unique(monocle3::clusters(spp1.data)))) 

monocle3::plot_cells(spp1.data,group_label_size = 4,cell_size = 1,cell_stroke = 0)+ggplot2::scale_color_manual(values = cols.cl)


spp1.data[['cluster.MP_SPP1']]= paste0('SPP1_MP_C',monocle3::clusters(spp1.data))
monocle3::plot_cells(spp1.data,label_cell_groups = F,
                     group_label_size = 4,cell_size = 1,cell_stroke = 0)+
  ggplot2::scale_color_manual(values = cols.cl)+
  theme_void()

monocle3::plot_cells(spp1.data,color_cells_by = 'From',group_label_size = 4,cell_size = 0.1,cell_stroke = 0)
marker_test_res = monocle3::top_markers(spp1.data,
                                        group_cells_by="cluster.MP_SPP1", genes_to_test_per_group=100,
                                        reference_cells=1000, cores=32)
cc=table(spp1.data[['From']],spp1.data[['cluster.MP_SPP1']])/apply(table(spp1.data[['From']],spp1.data[['cluster.MP_SPP1']]),1,sum)
pheatmap::pheatmap(cc,color = rev(hcl.colors(100,'Spectral')))
marker_test_res.filter = filter(marker_test_res, marker_test_p_value < 0.01) %>% group_by(cell_group)%>% top_n(20,mean_expression)
write.csv(marker_test_res.filter,file = './tables_v2/SPP1 sub-clustering markers.csv')
mLN.cluster.markers = filter(marker_test_res,cell_group == 'SPP1_MP_C12' & marker_test_p_value < 0.01) %>% top_n(20,mean_expression)

plot_genes_violin(spp1.data[c('IFI30','MMP12','CCL3','CCL4',
                              'IGHG4','FOS','RGS1','HSPA1A'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())
plot_genes_violin(spp1.data[c('CD274'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())


plot_genes_violin(Myeloid.data[c('CCL18','CCL13','SPP1','CXCL10'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())

plot_genes_violin(spp1.data[c('IFI30'),],
                  group_cells_by = 'cluster.MP_SPP1',
                  ncol = 4)+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.cl)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())

spp1.data[['IFI30_Level']] = ifelse(spp1.data[['cluster.MP_SPP1']] %in% c('SPP1_MP_C11','SPP1_MP_C12','SPP1_MP_C14',
                                                                          
                                                                          'SPP1_MP_C16','SPP1_MP_C18','SPP1_MP_C2'),
                                    'IFI30_High','IFI30_Low')
plot_genes_violin(spp1.data[c('CD274'),],
                  group_cells_by = 'IFI30_Level',
                  normalize = T)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())

plot_genes_violin(spp1.data[c('CD274'),],
                  group_cells_by = 'IFI30_Level',
                  normalize = T)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())

cor.test(assay(spp1.data)['IFI30',],assay(spp1.data)['CD274',])

plot_genes_by_group(spp1.data,as.character(mLN.cluster.markers$gene_id),
                    group_cells_by="cluster.MP_SPP1",
                    ordering_type="cluster_row_col",
                    max.size=4)+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))+
  theme(axis.text.x = element_text(angle = 90))
saveRDS(spp1.data,file = './variables_v2/spp1.data.Rds')
spp1.data = readRDS(file = './variables_v2/spp1.data.Rds')

####4. IFI30 based difference####
spp1.sc = NormalizeData(spp1.sc)
spp1.sc$IFI30_level = ifelse(spp1.sc@assays$RNA@data['IFI30',] > median(spp1.sc@assays$RNA@data['IFI30',]),
                             'High_IFI30','Low_IFI30')
Idents(spp1.sc)=spp1.sc$IFI30_level
VlnPlot(spp1.sc[,spp1.sc@assays$RNA@data['CD274',]>0],features = c('CD274'),group.by = 'IFI30_level',pt.size = 0.1,adjust = 1)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())
#spp1.sc = NormalizeData(spp1.sc)
IFI30.diff = FindMarkers(spp1.sc,ident.1 = 'High_IFI30',logfc.threshold = -Inf)
write.csv(IFI30.diff,file = './tables_v2/IFI30.diff.csv')
spp1.data[['IFI30_Level']] = spp1.sc$IFI30_level

IFI30.diff = read.csv(file = './tables_v2/IFI30.diff.csv',row.names = 1)
IFI30.diff = FindMarkers(spp1.sc,ident.1 = 'High_IFI30')
IFI30.log2FC = IFI30.diff$avg_log2FC
names(IFI30.log2FC)=rownames(IFI30.diff)
# IFI30.log2FC = apply(log2(spp1.sc@assays$RNA@data[,spp1.sc$IFI30_level == 'High_IFI30']+1),1,median)-
#   apply(log2(spp1.sc@assays$RNA@data[,spp1.sc$IFI30_level != 'High_IFI30']+1),1,median)
IFI30.diff = IFI30.diff[order(IFI30.diff$avg_log2FC),]
IFI30.diff.top = IFI30.diff[c(1:10,6841:6850),]
IFI30.diff.top$label = factor(rownames(IFI30.diff.top),levels=rownames(IFI30.diff.top))
ggplot(IFI30.diff.top,aes(x=avg_log2FC,y=label,fill=avg_log2FC))+geom_bar(stat = 'identity',width = 0.6)+
  theme_cowplot(font_size = 8)+
  scale_fill_gradientn(colours = rev(hcl.colors(100,'Spectral')))


library(clusterProfiler)
go.items = read.gmt(file('./c5.all.v2024.1.Hs.symbols.gmt'))
hallmark.items = read.gmt(file('./MSigDB_Hallmark.gmt'))
IFI30.log2FC.s = sort(IFI30.log2FC,decreasing = T)

IFI30.gsea = GSEA(IFI30.log2FC.s,minGSSize = 5,TERM2GENE = hallmark.items)

IFI30.gsea.res = IFI30.gsea@result
ridgeplot(IFI30.gsea,showCategory = 10,fill = 'NES')+
  theme_cowplot(font_size = 7)+
  scale_fill_gradientn(colors  = hcl.colors(100))
write.csv(IFI30.gsea.res,file = './tables_v2/IFI30 gsea res.csv')
IFI30.gsea.go = GSEA(IFI30.log2FC.s,minGSSize = 5,TERM2GENE = go.items[grepl('GOBP',go.items$term),])

IFI30.gsea.go.res = IFI30.gsea.go@result
ridgeplot(IFI30.gsea.go,showCategory = 15,fill = 'NES')+
  theme_cowplot(font_size = 7)+
  scale_fill_gradientn(colors  = hcl.colors(100))


IFI30.gsea.M12 = GSEA(IFI30.log2FC.s,minGSSize = 5,TERM2GENE = m12.df[,c(2,1)],pvalueCutoff = 1.2)
gseaplot(IFI30.gsea.M12,geneSetID = 'M2')


#####Phagocytosis####

Phagocytosis.markers = c('CD74','FCGR2B','GSN','HSP90AA1','ITGAM','ITGB1',
                         'MSR1','CCL2','TGM2','THBS1','TLR2','CD93','CTSL',
                         'TUBA1C','CEBPB','CTSD','CLEC4W','NFKBIA','TREM1')
pha.diff = IFI30.diff[Phagocytosis.markers,]
spp1.sc = AddModuleScore(spp1.sc,features = list(cc=Phagocytosis.markers),name = 'phagoScore')

VlnPlot(spp1.sc,features = c('phagoScore1'),group.by = 'IFI30_level',pt.size = 0,adjust = 3)+
  ggpubr::stat_compare_means()

meta.data.show = spp1.sc@meta.data
meta.data.show$cluster_C12 = ifelse(meta.data.show$cluster.MP_SPP1 == 'SPP1_MP_C12','C12','Others')
ggplot(meta.data.show,aes(x=cluster_C12,y=phagoScore1,fill=cluster_C12))+
  geom_violin(draw_quantiles = 0.5,scale = 'width')+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))
ggplot(meta.data.show,aes(x=From,y=phagoScore1,fill=From))+
  geom_violin(draw_quantiles = 0.5,scale = 'width')+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))


ggplot(meta.data.show,aes(x=IFI30_level,y=M1Score1,fill=IFI30_level))+
  geom_violin(draw_quantiles = 0.5,scale = 'width')+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))
ggplot(meta.data.show,aes(x=IFI30_level,y=M2Score1,fill=IFI30_level))+
  geom_violin(draw_quantiles = 0.5,scale = 'width')+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))

# ####5. M1-M2 scores####
# library(scGSVA)
# hsMS <- buildMSIGDB(species = "human", keytype = "SYMBOL", anntype = "HALLMARK")
# 
M1.markers = read.delim(file = './M1 markers.txt',header = F)
M2.markers = read.delim(file = './M2 markers.txt',header = F)

M1.markers = unique(as.character(M1.markers$V1))
M2.markers = unique(as.character(M2.markers$V1))
spp1.sc = AddModuleScore(spp1.sc,features = list(cc=M1.markers),name = 'M1Score')
spp1.sc = AddModuleScore(spp1.sc,features = list(cc=M2.markers),name = 'M2Score')

# 
m12.df = data.frame(GeneID = c(M1.markers,M2.markers),
                    Term = c(rep('M1',length(M1.markers)),
                             rep('M2',length(M2.markers))),
                    Annot = c(rep('M1',length(M1.markers)),
                              rep('M2',length(M2.markers)))
)
# 
# mLN.up = rownames(spp1.sc.diff[spp1.sc.diff$p_val_adj<0.01 & spp1.sc.diff$avg_log2FC> 1.5,])
# 
# mLN.up.M1M2 = enricher(mLN.up,TERM2GENE = m12.df[,c(2,1)])@result
# 
# anno.m12 = hsMS
# anno.m12@annot = m12.df
# 
# 
# gsva.m12 = scgsva(as.matrix(assay(Myeloid.data)),anno.m12,method="ssgsea",cores = 16)
# 
# Myeloid.data[['M1_score']] = gsva.m12@gsva$M1
# Myeloid.data[['M2_score']] = gsva.m12@gsva$M2
# 
# mp.data = Myeloid.data[,Myeloid.data[['subCellType_level1']] == 'MP']
# 
# meta.data.mp.show = as.data.frame(colData(mp.data))
# 
# ggplot(meta.data.mp.show,aes(x=subCellType,y=M1_score-M2_score,fill=subCellType))+geom_boxplot(outlier.size = 0.1)+
#   theme_cowplot(font_size = 8)+
#   scale_fill_manual(values = cols.cell)
# 
# ggplot(meta.data.mp.show,aes(x=From,y=M2_score,fill=From))+geom_boxplot(outlier.size = 0.1)+
#   theme_cowplot(font_size = 8)+
#   scale_fill_manual(values = cols.From)+
#   theme(axis.text.x = element_text(angle = 90))
# 
# spp1.data[['IFI30_level']]= spp1.sc$IFI30_level
# 
# meta.data.mp.show = as.data.frame(colData(spp1.data))
# 
# ggplot(meta.data.mp.show,aes(x=IFI30_level,y=M2_score,fill=IFI30_level))+geom_boxplot(outlier.size = 0.1)+
#   theme_cowplot(font_size = 8)
# ggplot(meta.data.mp.show,aes(x=IFI30_level,y=M1_score,fill=IFI30_level))+geom_boxplot(outlier.size = 0.1)+
#   theme_cowplot(font_size = 8)

####6. spatial analysis#####
st.m = readRDS('./variables_v2/st.m_bayesSpace.Rds')
st.m.myeloid = st.m[,st.m$cluster.name %in% c('Myeloid','Plasma+Myeloid')]
SpatialFeaturePlot(st.m,
                   features =c('SPP1'),slot = 'data',
                   ncol = 8,
                   alpha = 0.7,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0)
SpatialFeaturePlot(st.m.myeloid,
                   features =c('SPP1'),slot = 'data',
                   ncol = 8,
                   alpha = 0.7,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0)

SpatialFeaturePlot(st.m.myeloid,
                   features =c('CCL21'),slot = 'data',
                   images = c('P6L','P6T'),
                   ncol = 2,
                   alpha = 0.7,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0.6)


SpatialDimPlot(st.m.myeloid,group.by = 'cluster.name',
               ncol = 4,
               alpha = 0.8,
               pt.size.factor = 1.5,
               stroke = 0,
               image.alpha = 1,
               cols = cols.cl)

save.image('myeloid_v2.RData')


#####Spatial clustering#####
st.f.myeloid = st.m.myeloid[,st.m.myeloid$Source == 'PLM']
cds.st.myeloid = new_cell_data_set(as(as.matrix(st.f.myeloid@assays$Spatial@counts[rownames(st.f.myeloid),]), "sparseMatrix"),
                                   cell_metadata = st.f.myeloid@meta.data,
                                   gene_metadata = data.frame(id = rownames(st.f.myeloid),gene_short_name=rownames(st.f.myeloid),
                                                              row.names = rownames(st.f.myeloid))
)
cds.st.myeloid <- cds.st.myeloid[,Matrix::colSums(exprs(cds.st.myeloid)) != 0]
dim(cds.st.myeloid) ###26262  5259
cds.st.myeloid <- estimate_size_factors(cds.st.myeloid)
cds.st.myeloid <- monocle3::preprocess_cds(cds.st.myeloid, num_dim = 50,method="PCA",norm_method="log") 
plot_pc_variance_explained(cds.st.myeloid)


cds.st.myeloid <- align_cds(cds.st.myeloid, num_dim = 20, alignment_group = "orig.ident")
cds.st.myeloid <- reduce_dimension(cds.st.myeloid,cores = 20)
cds.st.myeloid <- cluster_cells(cds.st.myeloid, resolution=1e-3)

cols.cl = colorRampPalette(jcolors())(length(unique(clusters(cds.st.myeloid))))

plot_cells(cds.st.myeloid,group_label_size = 4,cell_size = 1,cell_stroke = 0)+
  ggplot2::scale_color_manual(values = cols.cl)


marker_test_res.epi.st <- top_markers(cds.st.myeloid, group_cells_by="cluster", genes_to_test_per_group=100,
                                      reference_cells=1000, cores=24)
marker_test_res.epi.st$cell_group = paste0('SC',marker_test_res.epi.st$cell_group)
write.csv(marker_test_res.epi.st,file = './tables_v2/myeloid marker_test_res.epi.st.csv')

marker_test_res.epi.st.filter = filter(marker_test_res.epi.st,marker_test_p_value < 0.01) %>% group_by(cell_group) %>% top_n(20,specificity)

cds.st.myeloid[['cluster_ml']] = paste0('SC',clusters(cds.st.myeloid))
cds.st.myeloid[['cluster_ml']] =factor(cds.st.myeloid[['cluster_ml']] ,levels = paste0('SC',1:7))

st.f.myeloid$cluster_ml = paste0('SC',clusters(cds.st.myeloid))

cols.st =  cols.default[1:length(unique(st.f.myeloid$cluster_ml))]
names(cols.st) = unique(st.f.myeloid$cluster_ml)
SpatialDimPlot(st.f.myeloid,group.by = 'cluster_ml',ncol = 4,image.alpha = 0.4,
               images = c('P3L','P5L','P6L','P7L','P3T','P5T','P6T','P7T'), 
               pt.size.factor = 1.5,stroke = 0,cols =cols.st)

SpatialDimPlot(st.f.myeloid,group.by = 'cluster_ml',ncol =2,image.alpha = 0.4,
               images = c('P6L','P6T'), 
               pt.size.factor = 1.8,stroke = 0,cols =cols.st)


#####Plot#####
plot_genes_violin(cds.st.myeloid[ c('CCL21'),],group_cells_by = 'cluster_ml',ncol = 4,normalize = F)+
  theme(axis.text.x = element_text(angle = 90))+scale_fill_manual(values = cols.st)

plot_genes_violin(cds.st.myeloid[ c('CCL21','LIPA','HLA-DRB5','ACP5','GPNMB'),],group_cells_by = 'cluster_ml',ncol = 1,normalize = T)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.st)


SpatialFeaturePlot(st.f.myeloid,
                   features =c('CCL21'),
                   images = c('P6L','P6T'), 
                   alpha = 0.8,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0.5,
                   ncol = 4)

SpatialFeaturePlot(st.f.myeloid,
                   features =c('SPP1'),
                   alpha = 0.8,
                   pt.size.factor = 1.8,
                   stroke = 0,
                   image.alpha = 0.5,
                   ncol = 2)

SpatialFeaturePlot(st.m,
                   features =c('GPNMB'),
                   images = c('P6L','P6T'), 
                   alpha = 0.8,
                   pt.size.factor = 1.5,
                   stroke = 0,
                   image.alpha = 0.5,
                   ncol = 4)
SpatialFeaturePlot(st.m.myeloid,
                   features =c('CCL21'),slot = 'data',
                   images = c('P6L','P6T'),
                   ncol = 2,
                   alpha = 0.7,
                   pt.size.factor = 1.8,
                   stroke = 0,
                   image.alpha = 0.6)

VlnPlot(st.m,features = c('GPNMB'),group.by = 'From',pt.size = 0)


plot_genes_violin(cds.st.m.myeloid[c('FCGR2A','FCGR2B','FCGR3A','FCGR3B'),],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

plot_genes_violin(Myeloid.data[c('GPNMB'),Myeloid.data[['From']] %in% c('mLN','mLN+ PT','mLN- PT')],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))
plot_genes_violin(Myeloid.data[c('LIPA'),Myeloid.data[['From']] %in% c('mLN','mLN+ PT','mLN- PT')],
                  group_cells_by = 'From',
                  ncol = 4)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))

saveRDS(cds.st.myeloid,file = './variables_v2/cds.st.myeloid.Rds')
cds.st.myeloid = readRDS(file = './variables_v2/cds.st.myeloid.Rds')
write.csv(marker_test_res.epi.st,file = './tables_v2/spatial myeloid clusters markers.csv')




plot_genes_violin(Myeloid.data[c('CCL21','LIPA','HLA-DRB5','ACP5','GPNMB','CCL18','CPVL','FCGR2B'),
                               Myeloid.data[['From']] %in% c('mLN','mLN+ PT','mLN- PT')],
                  group_cells_by = 'From',ncol = 4)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.From)
plot_genes_violin(Myeloid.data[c('CCL21','LIPA','HLA-DRB5','ACP5','GPNMB','CCL18','CPVL','FCGR2B'),
                               Myeloid.data[['From']] %in% c('mLN','mLN+ PT','mLN- PT')],
                  group_cells_by = 'subCellType',ncol = 1)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.cell)

plot_genes_by_group(cds.st.myeloid,as.character(marker_test_res.epi.st.filter[marker_test_res.epi.st.filter$cell_group == 'SC7',]$gene_id),
                    group_cells_by="cluster_ml",
                    ordering_type="cluster_row_col",
                    max.size=3)+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))+
  theme(axis.text.x = element_text(angle = 90))

#####Spatial cluster number diff#####
meta.data.cd = as.data.frame(colData(cds.st.myeloid))

f.data=table(meta.data.cd$orig.ident,meta.data.cd$cluster_ml)
f.data.r = melt(f.data,meansure.vars = colnames(f.data))
f.data.r$From = substr(f.data.r$Var1,3,3)

ggplot(f.data.r,aes(x=From,fill=From,y=value))+geom_boxplot(outlier.size = 0.3)+
  ggpubr::stat_compare_means()+
  facet_wrap(~Var2,scales = 'free')+
  theme_cowplot(font_size = 8)

#####Distance based analysis#####

st.f = st.m[,st.m$Source == 'PLM']
samples = unique(names(st.f@orig.ident))

unique(st.f$orig.ident)

myeloid.spots = c()
myeloid.distance.min = c()
myeloid.distance.median = c()
orig.flags = c()

for(sample in samples){
  
  st.f.i.myeloid = colnames(st.f[,st.f$orig.ident == sample & st.f$cluster.name == 'Myeloid'])
  st.f.i.epi = colnames(st.f[,st.f$orig.ident == sample & st.f$cluster.name == 'Epithelial'])
  
  st.f.i.distance = epiDistance(st.f,sample,st.f.i.myeloid,st.f.i.epi)
  
  
  
  myeloid.spots = append(myeloid.spots,st.f.i.myeloid)
  myeloid.distance.min = append(myeloid.distance.min,
                                apply(st.f.i.distance,1,min))
  myeloid.distance.median = append(myeloid.distance.median,
                                   apply(st.f.i.distance, 1, median))
  orig.flags = append(orig.flags,rep(sample,length(st.f.i.myeloid)))
  print(sample)
}

myeloid.distance.res = data.frame(row.names = myeloid.spots,
                                  distance.min = myeloid.distance.min,
                                  distance.median = myeloid.distance.median,
                                  sample = orig.flags)
min(myeloid.distance.res$distance.min)
myeloid.distance.res$distance.tag = ifelse(myeloid.distance.res$distance.min> 3,'Far','Near')
saveRDS(myeloid.distance.res,file = './variables_v2/myeloid.distance.res-v2.Rds')
st.f.myeloid.2 = st.f[,st.f$cluster.name == 'Myeloid']
st.f.myeloid.2 = NormalizeData(st.f.myeloid.2)
st.f.myeloid.2$distance.tag = myeloid.distance.res[colnames(st.f.myeloid.2),'distance.tag']
VlnPlot(st.f.myeloid.2,features = 'SPP1',group.by = 'distance.tag')

Idents(st.f.myeloid.2) = st.f.myeloid.2$distance.tag

st.f.myeloid.diff.dis = FindMarkers(st.f.myeloid.2,ident.1 = 'Near')
cds.st.f.myeloid.2 = new_cell_data_set(as(as.matrix(st.f.myeloid.2@assays$Spatial@counts[rownames(st.f.myeloid.2),]), "sparseMatrix"),
                                     cell_metadata = st.f.myeloid.2@meta.data,
                                     gene_metadata = data.frame(id = rownames(st.f.myeloid.2),gene_short_name=rownames(st.f.myeloid.2),
                                                                row.names = rownames(st.f.myeloid.2))
)

cds.st.f.myeloid.2[['From']] = factor(cds.st.f.myeloid.2[['From']],levels = c('mLN','mLN+ PT','mLN- PT'))
plot_genes_violin(cds.st.f.myeloid.2[c('SPP1','CXCL9','APOE','GPNMB','IFI30','PGK1'),],
                  group_cells_by = 'distance.tag',ncol = 3)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())
saveRDS(st.f.myeloid.2,file = './variables_v2/st.f.myeloid.2.Rds')
st.f.myeloid.2 = readRDS(file = './variables_v2/st.f.myeloid.2.Rds')
