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

T.data = readRDS(file = './variables_v2/T.data.Rds')

marker.genes = c('CD3E','CD3D','CD4','CD8A','CD8B',
                 'FCGR3A','XCL1','KLRF1',
                 'IL7R','SELL','TCF7','CCR7','LEF1',
                 'RORA','CD69','CD40LG',
                 'LAG3','TIGIT','PDCD1','HAVCR2','CTLA4','RGS1',
                 'GZMA','GNLY','PRF1','GZMB','GZMK','IFNG','NKG7',
                 'IL2RA','FOXP3',
                 'S100A4','S100A6','IL32','ANXA1','CCL5',
                 'CXCL13','ARHGAP15','ALDOA',
                 'TRDC',
                 'BCL6','CXCR5')
plot_genes_by_group(T.data,marker.genes,
                    group_cells_by="subCellType",
                    ordering_type="cluster_row_col",norm_method = 'size_only',
                    max.size=4,scale_max = 10)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))
T.data[['subCellType']] = factor(T.data[['subCellType']],levels = c('ALDOA+ T',"ARHGAP15+ T",'CD4 Naive T','CD4 Resident Effector Memory',
                                                                    'CD4 Effector Memory', "CXCL13+ T", "Exhausted T",'Treg',
                                                                    "CD8 Cytotoxic T","CD8 Exhausted T",'NK',"¦Ã¦Ä T"))


cols.cell = cols.default[1:length(unique(colData(T.data)$subCellType))]
names(cols.cell)=levels(colData(T.data)$subCellType)
scater::plotReducedDim(T.data,dimred = 'UMAP',colour_by = 'subCellType',
                       point_size = 1)+
  ggplot2::scale_color_manual(values = cols.cell)+
  theme_void()

####2. statistics on the cell sub types####

cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')

meta.data.show = T.data@colData
meta.data.show = as.data.frame(meta.data.show)
meta.data.show$Disease2 = substr(meta.data.show$Disease,1,4)
meta.data.show$subCellType = as.character(meta.data.show$subCellType)


ggplot2::ggplot(meta.data.show,ggplot2::aes(x=From,fill=subCellType))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+facet_wrap(~Disease2,ncol = 1)+
  theme_cowplot(font_size = 8)+
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
  facet_wrap(~Var2,scales = 'free',ncol = 6)+
  theme_cowplot(font_size = 8)+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 0))+
  scale_fill_manual(values = cols.From)



####3. From difference####

marker_test_res.T = top_markers(T.data,group_cells_by = 'subCellType',genes_to_test_per_group = 100,reference_cells = 1000,cores = 32)
marker_test_res.T.filter = filter(marker_test_res.T,fraction_expressing > 0.5 & marker_test_p_value < 0.01) %>% group_by(cell_group)

cxcl13.markers = filter(marker_test_res.T,marker_test_p_value < 0.01 & 
                               substr(gene_id,1,2) %in% c('MT','RP') == F & 
                               cell_group == 'CXCL13+ T' &
                               gene_id != 'MALAT1') %>% 
  arrange(marker_test_p_value)
cxcl13.markers$label = ''
cxcl13.markers$label[1:20]=cxcl13.markers$gene_id[1:20]
ggplot(cxcl13.markers,aes(color=log2(mean_expression+1),y=-log10(marker_test_p_value),x=specificity))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label = label),size=2,max.overlaps = 1000)+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours = hcl.colors(100))






####4. exhaustion scores####
st.t = CreateSeuratObject(counts = assay(T.data))
st.t = NormalizeData(st.t)

ex.markers = c("CXCL13", "HAVCR2", "PDCD1", "TIGIT", "LAG3", "CTLA4", "LAYN", "RBPJ", "VCAM1", "GZMB", "TOX", "MYO7A")
ef.markers = c('PRF1', 'IFNG', 'CCL4', 'HLA-DQA1', 'GZMK', 'GZMA', 'GZMH', 'CD44', 'DUSP2', 'KLRB1', 'KLRD1' , 'CTSW')
st.t = AddModuleScore(st.t,features = list(cc=ex.markers),name = 'exScore')
T.data[['exScore']]= st.t$exScore1
#T.data[['efScore']]= apply(st.t@assays$RNA@data[ef.markers,], 2, mean)

meta.data.show = as.data.frame(colData(T.data))

saveRDS(meta.data.show,file='./tables_v2/T_meta.data.show.Rds')

ggplot2::ggplot(meta.data.show[!grepl('NK',meta.data.show$subCellType),],
                ggplot2::aes(x=reorder(subCellType,exScore,FUN=median),y=exScore,fill=subCellType))+
  ggplot2::geom_boxplot(outlier.size = 0.5)+
  theme_cowplot()+
  ggplot2::scale_fill_manual(values = cols.cell)+
  theme(axis.text  = element_text(size = 8),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))
st.t = AddMetaData(st.t,as.data.frame(colData(T.data)))
Idents(st.t) = st.t$subCellType

cxcl13.markers.2 = FindMarkers(st.t,ident.1 = 'CXCL13+ T')
write.csv(cxcl13.markers.2,file = './tables_v2/seurat cxcl13 t markers.csv')

cxcl13.markers.2.sig = cxcl13.markers.2[cxcl13.markers.2$p_val_adj<0.01,]
cxcl13.markers.2.sig = cxcl13.markers.2.sig[order(cxcl13.markers.2.sig$avg_log2FC),]

cxcl13.markers.top = cxcl13.markers.2.sig[c(1:10,680:689),]
cxcl13.markers.top$gene = factor(rownames(cxcl13.markers.top),levels = rownames(cxcl13.markers.top))
ggplot(cxcl13.markers.top,aes(x=avg_log2FC,y=gene,fill=avg_log2FC))+geom_bar(stat = 'identity',width = 0.6)+
   scale_fill_gradientn(colors = hcl.colors(100,'Spectral'))+
  theme_cowplot(font_size = 8)


#pathway
library(clusterProfiler)

go.items = read.gmt(file('./c5.all.v2024.1.Hs.symbols.gmt'))
hallmark.items = read.gmt(file('./MSigDB_Hallmark.gmt'))
c12.log2FC = cxcl13.markers.2$avg_log2FC
names(c12.log2FC)=rownames(cxcl13.markers.2)
c12.log2FC.s = sort(c12.log2FC,decreasing = T)
c12.gsea = GSEA(c12.log2FC.s,minGSSize = 5,TERM2GENE = go.items,eps = 0)
c12.gsea.res = c12.gsea@result
#ridgeplot(c12.gsea)
c12.gsea.res.gobp = c12.gsea.res[grepl('GOBP',c12.gsea.res$ID),] 
c12.gsea.res.gobp$Type = ifelse(c12.gsea.res.gobp$NES>0,'Up','Down')
c12.gsea.res.gobp.top = filter(c12.gsea.res.gobp,p.adjust<0.05) %>% group_by(Type) %>% top_n(20,abs(NES))
c12.gsea.res.gobp.top$ID = gsub('GOBP_','',c12.gsea.res.gobp.top$ID)
c12.gsea.res.gobp.top$ID =  factor(c12.gsea.res.gobp.top$ID,levels = c12.gsea.res.gobp.top$ID)
ggplot(c12.gsea.res.gobp.top,aes(x=NES,y=ID,fill=-log10(p.adjust)))+
  geom_bar(stat = 'identity',width = 0.65)+
  scale_fill_gradientn(colours = rev(hcl.colors(100)))+
  theme_cowplot(font_size = 8)

c12.gsea.hallmark = GSEA(c12.log2FC.s,minGSSize = 5,TERM2GENE = hallmark.items,eps = 0)
c12.gsea.res.h = c12.gsea.hallmark@result
gseaplot(c12.gsea.hallmark,geneSetID = 'TNFA_SIGNALING_VIA_NFKB')


####5. cytoTRACE analysis####

expr.T = as.matrix(T.data@assays@data$counts)

# T.data.LUSC = preprocess_cds(T.data.LUSC)
# T.data.LUSC = reduce_dimension(T.data.LUSC,cores = 20)
# 
selected_genes <- rownames(expr.T)[ Matrix::rowSums(T.data@assays@data$counts > 0)> 100]
cyto.T <- CytoTRACE(expr.T[selected_genes,],batch = as.character(colData(T.data)$Source),ncores = 12)
#plotCytoGenes(cyto.LUSC, numOfGenes = 10)

gc()

cyto.genes = cyto.T$cytoGenes


sub.cells = sample(colnames(T.data),500)
cyto.genes.temp = cyto.genes[substr(names(cyto.genes),1,2) %in% c('MT','RP')==F]
top.genes = c(cyto.genes.temp[1:25],rev(cyto.genes.temp)[1:25])
plot.mat = log2( as.data.frame(t(expr.T[names(top.genes),sub.cells]))+1)

plot.mat$CytoTRACE = as.numeric(cyto.T[['CytoTRACE']][sub.cells])
plot.mat = plot.mat[order(plot.mat$CytoTRACE),]

pheatmap::pheatmap(t(as.data.frame(lapply(plot.mat,minMax),check.names = F))[c('CytoTRACE',names(top.genes)),],color = rev(hcl.colors(100,'Spectral')),show_colnames = F,cluster_cols = F,
                   cluster_rows = F,gaps_row = c(1,26),fontsize = 8)


T.data[['CytoTRACE']] = cyto.T[['CytoTRACE']]
rm(expr.T)



meta.data.T = as.data.frame(colData(T.data))
ggplot(meta.data.T,aes(y=CytoTRACE,x=subCellType,fill=subCellType))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values =cols.cell)+
  theme(axis.text.x = element_text(angle=90))


save.image(file = 'T_image_v2.RData')

####6.Spatial analysis####

st.m = readRDS(file='./variables_v2/st.m_bayesSpace.Rds')
st.f.ls = st.m[,st.m$cluster.name == 'LS-like']
SpatialDimPlot(st.m,group.by = 'cluster.name',ncol =2,image.alpha = 0.6,
               alpha = 0.7,
               images = c('P6L','P6T'), 
               pt.size.factor = 1.8,stroke = 0)

VlnPlot(st.m,features = c('CXCL13','ICA1','IL6ST','NMB','FKBP5','PDCD1','TOX2','CTLA4'),
        group.by = 'cluster.name',pt.size = 0,ncol = 4,slot='count',adjust = 8)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())
st.m = AddModuleScore(st.m,features = list(cc=ex.markers),name = 'exScore')
meta.data.st = st.m@meta.data
ggplot(meta.data.st,aes(x=From,y=exScore1,fill=From))+geom_boxplot()+
  theme_cowplot(font_size = 8)+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(2,3),c(1,3)))

cds.st.m = new_cell_data_set(as(as.matrix(st.m@assays$Spatial@counts[rownames(st.m),]), "sparseMatrix"),
                             cell_metadata = st.m@meta.data,
                             gene_metadata = data.frame(id = rownames(st.m),gene_short_name=rownames(st.m),
                                                        row.names = rownames(st.m))
)
plot_genes_violin(cds.st.m[c('CXCL13','ICA1','IL6ST','NMB','FKBP5','PDCD1','TOX2','CTLA4'),],
                  group_cells_by = 'cluster.name',ncol = 4,normalize = F)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

cds.st.m[['cluster.2']] = ifelse(cds.st.m[['cluster.name']] == 'LS-like','LS-like','Others')
plot_genes_violin(cds.st.m[c('CXCL13','TRBC2','PDCD1','CTLA4'),],
                  group_cells_by = 'cluster.2',ncol = 4,normalize = F)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())
