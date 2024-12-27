.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM/')
library(BayesSpace)
library(scater)
#library(harmony)
library(dplyr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(reshape2)
library(Seurat)
####1. read and process by seurat####
st.m = readRDS(file = '/share/dichen/lungNodeM/st_E13530&Our.Rds')

cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559","#32709f",
                 '#b8a3c4','#60428d','#9f5633','grey','#c92e2a','#c79795','#d9c17c',"#d0896b",'#5d94b4',"black","blue",'cyan')



image.info = data.frame()
for(i in 1:length(unique(st.m$region))){
  
  image.i = st.m@images[[i]]
  image.i = image.i@coordinates
  
  image.info = rbind(image.info,image.i)
}

st.m = AddMetaData(st.m,image.info)

st.m = SCTransform(st.m,variable.features.n = 2000,assay = "Spatial")
st.m = RunPCA(st.m, verbose = T, npcs = 30)
st.m = harmony::RunHarmony(st.m,group.by.vars = 'region',assay.use = 'SCT')
st.m = RunUMAP(st.m,reduction = "harmony", dims = 1:30)
#saveRDS(st.m,file = './variables/st.merge_E13530.Rds')


####2. Change as one bayesSpace object####
colData <- st.m@meta.data

colData$spot <- rownames(colData)  
colData$in_tissue = colData$tissue
#colData <- colData[colData$in_tissue > 0, ]

counts <- Matrix::as.matrix(st.m@assays$Spatial@counts)

rowData <- data.frame(gene_id = rownames(counts),gene_name = rownames(counts))

sce <- SingleCellExperiment(assays = list(counts = counts), 
                            rowData = rowData, colData = colData)
metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "Visium"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

clusterPlot(sce, label="region", color = NA)  #make sure no overlap between samples

sce$row[sce$region == "P11_T1"] = 
  100 + sce$row[sce$region == "P11_T1"]
sce$row[sce$region == "P15_T2"] = 
  200 + sce$row[sce$region == "P15_T2"]
sce$row[sce$region == "P16_T1"] = 
  300 + sce$row[sce$region == "P16_T1"]

sce$row[sce$region == "P19_T1"] = 
  100 + sce$row[sce$region == "P19_T1"]
sce$row[sce$region == "P24_T1"] = 
  200 + sce$row[sce$region == "P24_T1"]
sce$row[sce$region == "P25_T2"] = 
  300 + sce$row[sce$region == "P25_T2"]

sce$row[sce$region == "P3T"] = 
  100 + sce$row[sce$region == "P3T"]
sce$row[sce$region == "P5L"] = 
  200 + sce$row[sce$region == "P5L"]
sce$row[sce$region == "P5T"] = 
  300 + sce$row[sce$region == "P5T"]

sce$row[sce$region == "P6T"] = 
  100 + sce$row[sce$region == "P6T"]
sce$row[sce$region == "P7L"] = 
  200 + sce$row[sce$region == "P7L"]
sce$row[sce$region == "P7T"] = 
  300 + sce$row[sce$region == "P7T"]



sce$col[sce$region == "P17_T2"] = 
  150 + sce$col[sce$region == "P17_T2"]
sce$col[sce$region == "P19_T1"] = 
  150 + sce$col[sce$region == "P19_T1"]
sce$col[sce$region == "P24_T1"] = 
  150 + sce$col[sce$region == "P24_T1"]
sce$col[sce$region == "P25_T2"] = 
  150 + sce$col[sce$region == "P25_T2"]

sce$col[sce$region == "P3L"] = 
  300 + sce$col[sce$region == "P3L"]
sce$col[sce$region == "P3T"] = 
  300 + sce$col[sce$region == "P3T"]
sce$col[sce$region == "P5L"] = 
  300 + sce$col[sce$region == "P5L"]
sce$col[sce$region == "P5T"] = 
  300 + sce$col[sce$region == "P5T"]


sce$col[sce$region == "P6L"] = 
  450 + sce$col[sce$region == "P6L"]
sce$col[sce$region == "P6T"] = 
  450 + sce$col[sce$region == "P6T"]
sce$col[sce$region == "P7L"] = 
  450 + sce$col[sce$region == "P7L"]
sce$col[sce$region == "P7T"] = 
  450 + sce$col[sce$region == "P7T"]



clusterPlot(sce, "region", color = NA)  #make sure no overlap between samples

####3. preprocess####
sce = spatialPreprocess(sce, n.PCs = 50) #lognormalize, PCA
sce = runUMAP(sce, dimred = "PCA")
colnames(reducedDim(sce, "UMAP")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce$region))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()

reducedDim(sce,'HARMONY') = st.m@reductions$harmony@cell.embeddings
reducedDim(sce,'UMAP.HARMONY') = st.m@reductions$umap@cell.embeddings

colnames(reducedDim(sce, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce$region))) +
  geom_point(size=0.1) +
  labs(color = "Sample") +
  theme_bw()
####4. clustering####
# sce <- qTune(sce, qs=seq(5, 25),use.dimred = "HARMONY", platform="Visium", d=20,
#              gamma=3,
#              nrep=1000, burn.in=100)
# 
# qPlot(sce)

sce = spatialCluster(sce, use.dimred = "HARMONY", q = 21, 
                     platform="Visium", d=20,
                     init.method="mclust", model="t", gamma=3,
                     nrep=1000, burn.in=100) #use HARMONY
print(Sys.time())#16:00 to ?"2024-12-19 16:45:53 +08"
clusterPlot(sce, color = NA,
            size=0.1) + #plot clusters
  labs(title = "BayesSpace joint clustering")
saveRDS(sce,'./variables_v2/sce_bayesSpace.Rds')
#saveRDS(sce,'./variables/sce_bayesSpace_withMPLC.Rds')
sce = readRDS('./variables_v2/sce_bayesSpace.Rds')

featurePlot(sce,feature = 'CCL19',color = NA)+scale_fill_gradientn(colors = hcl.colors(100,'Inferno'))

st.m$bs.cluster = paste0('SC',sce[['spatial.cluster']])

#SpatialDimPlot(st.m,group.by = 'bs.cluster')
######Sub-clustering#####
sce = readRDS('./variables_v2/sce_bayesSpace.Rds')

sce.sub = sce[,sce[['cluster.name']] == 'LS-like']
BayesSpace::clusterPlot(sce.sub, "region", color = NA)  #make sure no overlap between samples
featurePlot(sce.sub,feature = 'CXCL13',color=NA)
sce.sub = spatialCluster(sce.sub, use.dimred = "HARMONY", q = 2, 
                     platform="Visium", d=10,
                     init.method="mclust", model="t", gamma=3,
                     nrep=1000, burn.in=100) #use HARMONY

clusterPlot(sce.sub, color = NA,
            palette = ggsci::pal_npg()(2),
            size=0.1) + #plot clusters
  labs(title = "BayesSpace joint clustering")

####5. Find markers####

Idents(st.m)=st.m$bs.cluster

de_markers = FindAllMarkers(st.m, assay = 'SCT',logfc.threshold = 0.1, 
                            test.use = "wilcox", 
                            min.pct = 0.01, 
                            min.diff.pct = 0.1, 
                            only.pos = TRUE)
#de_markers$cluster = paste0('SC',de_markers$cluster)
top_markers = filter(de_markers,p_val_adj<0.01) %>% group_by(cluster) %>% top_n(6,avg_log2FC) %>% filter(pct.1>0.3)

DoHeatmap(st.m[,sample(colnames(st.m),2500)],features = unique(top_markers$gene))
DotPlot(st.m,assay='SCT',features =  unique(top_markers$gene))+coord_flip()+
  theme_cowplot(font_size = 6)+
  scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))+
  scale_size_continuous(range = c(0,2.5))+
  theme(axis.text.x = element_text(angle = 90))


write.csv(de_markers,file = './tables_v2/top markers of spatial clusters_bayesspace SC21.csv')

####6. Statistics on clusters####
cols.cl = cols.default[1:21]
names(cols.cl) = paste0(
  'SC',1:21
)


st.m$bs.cluster = factor(st.m$bs.cluster,levels = names(cols.cl))


ggplot(st.m@meta.data,aes(x=region,fill=bs.cluster))+geom_bar(position = 'fill',width = 0.75)+
  theme_cowplot()+
  scale_fill_manual(values = cols.cl)+
  theme(axis.text.x = element_text(angle = 90))


com2orig = table(st.m$region,st.m$bs.cluster)

f.data<-as.matrix(com2orig)
f.data = f.data/apply(f.data, 1, sum)
#f.data = t(f.data)
#rownames(ids.all)=ids.all$id
meta.data.u = st.m@meta.data
meta.data.u$orig.ident = meta.data.u$region
meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
rownames(meta.data.u)=meta.data.u$orig.ident

com2orig.df = melt(f.data,measure.vars = colnames(f.data))
#com2orig.df$patient = meta.data.u[as.character(com2orig.df$Var1),'PatientID']
com2orig.df$From = meta.data.u[as.character(com2orig.df$Var1),'From']

ggplot(com2orig.df,aes(x=From,y=value,fill=From,color=From))+geom_boxplot(outlier.size = 1,outlier.alpha = 0.7)+
  facet_wrap(~Var2,scales = 'free',ncol = 9)+
  ggpubr::stat_compare_means(aes(label = ..p.signif..),comparisons = list(c(1,2),c(1,3),c(2,3)))+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 0))+
  scale_fill_manual(values = alpha(brewer.pal(7,'Set1'),0.7))+
  scale_color_manual(values = brewer.pal(7,'Set1'))



####~~~~Plot some features####

featurePlot(sce, "CCL21",platform = 'Visium',size=0,color=NA)+scale_fill_gradientn(colours = hcl.colors(256))
featurePlot(sce, "CXCL13",platform = 'Visium',size=0,color=NA)+scale_fill_gradientn(colours = hcl.colors(256))
featurePlot(sce, "CST1",platform = 'Visium',size=0.5,color=NA)+scale_fill_gradientn(colours = hcl.colors(256))
featurePlot(sce, "KRT18",platform = 'Visium',size=0.5,color=NA)+scale_fill_gradientn(colours = hcl.colors(256))
featurePlot(sce, "MT-ND5",platform = 'Visium',size=0.5,color=NA)+scale_fill_gradientn(colours = hcl.colors(256))
featurePlot(sce, "RGS1",platform = 'Visium',size=0.5,color=NA)+scale_fill_gradientn(colours = hcl.colors(256))
featurePlot(sce, "IFI30",platform = 'Visium',size=0.5,color=NA)+scale_fill_gradientn(colours = hcl.colors(256))


st.m = PercentageFeatureSet(st.m, "^MT-", col.name = "percent_mito")
VlnPlot(st.m,features = "percent_mito",group.by = 'bs.cluster')
VlnPlot(st.m,features = c('FTH1','FTL'),group.by = 'From')


####7. Name the clusters####
st.m$cluster.name = dplyr::recode(st.m$bs.cluster,
                                  'SC1'='Epithelial',#
                                  'SC2'='Epithelial',#
                                  'SC3'='LS-like',
                                  'SC4'= 'LS-like',
                                  'SC5'='Fibroblast',#Fibroblast
                                  'SC6'='Fibroblast',
                                  'SC7'='Epithelial',
                                  'SC8'='Myeloid',#IFI30
                                  'SC9'='Epithelial',#
                                  'SC10'='Plasma',
                                  'SC11'='MT-high',
                                  'SC12'='Epithelial',
                                  'SC13'='Epithelial',
                                  'SC14'='Plasma',#Fibroblast
                                  'SC15'='Epithelial',
                                  'SC16'='Plasma',
                                  'SC17'='Fibroblast',
                                  'SC18'='Plasma+Myeloid',#IGHG,C1QB,CD68
                                  'SC19'='Plasma+Myeloid',
                                  'SC20'='Epithelial',
                                  'SC21'='Epithelial')#
sce[['cluster.name']] = st.m$cluster.name

cols.cl = c('lightgrey',cols.default[2:7])
names(cols.cl)=c('MT-high','LS-like','Plasma','Epithelial','Plasma+Myeloid','Myeloid','Fibroblast')
clusterPlot(sce, color = NA,label='cluster.name',
            palette = cols.cl,
            size=0.1) + #plot clusters
  labs(title = "BayesSpace joint clustering")


VlnPlot(st.m[,st.m$cluster.name == 'Epithelial'],features = c('FTH1','FTL'),group.by = 'From',
        pt.size = 0)

saveRDS(sce,'./variables_v2/sce_bayesSpace.Rds')
saveRDS(st.m,file='./variables_v2/st.m_bayesSpace.Rds')

st.m = readRDS(file='./variables_v2/st.m_bayesSpace.Rds')
SpatialDimPlot(st.m,cols = cols.cl,pt.size.factor = 1.2,ncol = 4,group.by = 'cluster.name',
               stroke = 0)
Idents(st.m) = st.m$cluster.name
cluster.markers = FindAllMarkers(st.m,slot = 'data')
write.csv(cluster.markers,file = './tables_v2/ST all spot cellType markers.csv')
cluster.markers = read.csv(file = './tables_v2/ST all spot cellType markers.csv',row.names = 1)
top_markers = filter(cluster.markers,p_val_adj<0.01) %>% group_by(cluster) %>% top_n(8,avg_log2FC)

DotPlot(st.m,assay='SCT',features =  unique(top_markers$gene))+coord_flip()+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours = rev(hcl.colors(100,'Spectral')))+
  scale_size_continuous(range = c(0,3))+
  theme(axis.text.x = element_text(angle = 90))


SpatialDimPlot(st.m,images=c('P3L12','P5L','P6L','P7L5'),image.alpha = 0,
               cols = cols.cl,pt.size.factor = 1.5,ncol = 4,group.by = 'cluster.name',
               stroke = 0)
SpatialDimPlot(st.m,images=c('P3T','P5T','P6T','P7T'),image.alpha = 0,
               cols = cols.cl,pt.size.factor = 1.5,ncol = 4,group.by = 'cluster.name',
               stroke = 0)

SpatialDimPlot(st.m,images=c('P10_T1','P11_T1','P15_T2','P16_T1',
                             'P17_T2','P19_T1','P24_T1','P25_T2'),image.alpha = 0,
               cols = cols.cl,pt.size.factor = 1.5,ncol = 4,group.by = 'cluster.name',
               stroke = 0)

SpatialFeaturePlot(st.m,images=c('P3T','P5T','P6T','P7T'),image.alpha = 0,
                   pt.size.factor = 1.5,ncol = 4,
                   feature = c('CCL19','MS4A1'),
                   stroke = 0)
SpatialFeaturePlot(st.m,images=c('P3L12','P5L','P6L','P7L5'),image.alpha = 0,
                   pt.size.factor = 1.5,ncol = 4,
                   feature = c('CCL19','MS4A1'),
                   stroke = 0)

####8. Statistics on cluster.name####



com2orig = table(st.m$region,st.m$cluster.name)

f.data<-as.matrix(com2orig)
f.data = f.data/apply(f.data, 1, sum)
#f.data = t(f.data)
#rownames(ids.all)=ids.all$id
meta.data.u = st.m@meta.data
meta.data.u$orig.ident = meta.data.u$region
meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
rownames(meta.data.u)=meta.data.u$orig.ident

com2orig.df = melt(f.data,measure.vars = colnames(f.data))
#com2orig.df$patient = meta.data.u[as.character(com2orig.df$Var1),'PatientID']
com2orig.df$From = meta.data.u[as.character(com2orig.df$Var1),'From']
cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')
com2orig.df$From = factor(com2orig.df$From,levels = c('mLN','mLN+ PT','mLN- PT'))
ggplot(com2orig.df,aes(x=From,y=value,fill=From))+
  geom_boxplot(outlier.size = 1,outlier.alpha = 0.7)+
  facet_wrap(~Var2,scales = 'free',ncol = 9)+
  ggpubr::stat_compare_means(comparisons = list(c(1,3),c(2,3)),method = 't.test',method.args = list(alternative = 'greater'))+
  theme_cowplot()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 0))+
  scale_fill_manual(values = cols.From)

####9.neighborhood analysis, update 2024-08-28#### 

nei.matrix = data.frame()
names(st.m@images) = unique(st.m$orig.ident)
orig.idents = unique(st.m$orig.ident)
for(orig.ident.i in orig.idents){
  print(orig.ident.i)
  st.i = st.m[,st.m$orig.ident == orig.ident.i]
  spots.info = st.m@images[[orig.ident.i]]@coordinates
  colnames(spots.info)=c('in_tissue','array_row','array_col','pxl_col','pxl_row')
  spots.info$pairxy=paste(spots.info$array_row,spots.info$array_col)
  
  n.spot = ncol(st.i)
  #spot.matrix = matrix(0,nrow = n.spot,ncol = n.spot,dimnames = list(colnames(st.i),colnames(st.i)))
  
  
  spot.matrix = Reduce(cbind,lapply(colnames(st.i),function(a){
    
    coord.i = spots.info[a,]
    row.i = as.integer(coord.i['array_row'])
    col.i = as.integer(coord.i['array_col'])
    candi.coord.i=data.frame(rowc=c(row.i,row.i,row.i-1,row.i-1,row.i+1,row.i+1),
                             colc=c(col.i-2,col.i+2,col.i-1,col.i+1,col.i-1,col.i+1))
    candi.coord.i$pairxy = paste(candi.coord.i$rowc,candi.coord.i$colc)
    
    neis.a = rownames(spots.info[spots.info$pairxy %in% as.character(candi.coord.i$pairxy),])
    
    one.row = data.frame(row.names = colnames(st.i),V1=rep(0,ncol(st.i)))
    one.row[neis.a,'V1'] = 1
    
    return(one.row)
    
  }))
  
  spot.matrix = as.matrix(spot.matrix)
  
  clusters = unique(st.m$cluster.name)
  clusters = clusters[!is.na(clusters)]
  meta.data = st.i@meta.data
  cluster.matrix = Reduce(cbind,lapply(clusters,function(a){
    
    cluster.spots = rownames(meta.data[meta.data$cluster.name == a,])
    cluster.spots = cluster.spots[cluster.spots != 'NA']
    one.row = data.frame(row.names = colnames(st.i),V1=rep(0,ncol(st.i)))
    one.row[cluster.spots,'V1'] = 1
    colnames(one.row) = a
    
    return(one.row)
    
  }))
  
  
  
  cluster.matrix = as.matrix(cluster.matrix)
  
  CC.nei.matrix = t(cluster.matrix) %*% spot.matrix %*% cluster.matrix
  
  diag(CC.nei.matrix) = 0.5*diag(CC.nei.matrix)
  
  
  nei.matrix.df = melt(CC.nei.matrix,measure.vars = colnames(CC.nei.matrix))
  nei.matrix.df$orig.ident = orig.ident.i
  
  nei.matrix = rbind(nei.matrix,nei.matrix.df)
  
}

nei.matrix$Patient = substr(nei.matrix$orig.ident,1,2)
nei.matrix$Source = ifelse(nei.matrix$Patient %in% c('P3','P5','P6','P7'),'PLM','Ana')
nei.matrix$From = ifelse(nei.matrix$Patient %in% c('P3','P5','P6','P7') == F,'mLN- PT',ifelse(grepl('L',nei.matrix$orig.ident),'mLN','mLN+ PT'))

nei.matrix$From = factor(nei.matrix$From,levels=c('mLN','mLN+ PT','mLN- PT'))
ggplot(nei.matrix,aes(x=Var1,y=Var2))+geom_raster(aes(fill=log2(value+1)))+facet_grid(~From)+
  scale_fill_gradientn(colours = rev(hcl.colors(100,'Spectral')))+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank(),
        legend.position = 'bottom')

write.csv(nei.matrix,file = './tables_v2/neighborhood_nei.matrix.csv')
nei.matrix = read.csv(file = './tables_v2/neighborhood_nei.matrix.csv',row.names = 1)
nei.matrix$From = factor(nei.matrix$From,levels=c('mLN','mLN+ PT','mLN- PT'))

diag.merge.res = filter(nei.matrix,Var1 == Var2)
cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')

ggplot(diag.merge.res[diag.merge.res$Var1 == 'LS-like',],aes(x=From,y=log2(value+1),fill=From))+geom_violin(draw_quantiles = 0.5,scale = 'width')+
  geom_jitter()+
  theme_cowplot()+
  ggpubr::stat_compare_means(comparisons = list(c(1,3),c(2,3)),method = 't.test',method.args = list(alternative = 'greater'))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.From)

plasma.epi = filter(nei.matrix,Var1 %in% c('Epithelial') )
plasma.sum = group_by(plasma.epi,orig.ident) %>% summarise(TN=sum(value))
plasma.sum = as.data.frame(plasma.sum)
rownames(plasma.sum) = plasma.sum$orig.ident
plasma.epi$prop = unlist(lapply(1:nrow(plasma.epi), function(i){
  cc = as.double(plasma.epi[i,'value'])/as.integer(plasma.sum[as.character(plasma.epi[i,'orig.ident']),'TN'])
  
}
       ))


ggplot(plasma.epi[plasma.epi$Var2 == 'Plasma',],aes(x=From,y=prop,fill=From))+geom_violin(draw_quantiles = 0.5,scale = 'width')+
  geom_jitter()+
  theme_cowplot()+
  ggpubr::stat_compare_means(comparisons = list(c(1,3),c(2,3)),method = 't.test',method.args = list(alternative = 'greater'))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.From)

fibro.epi = filter(nei.matrix,Var2 == 'Epithelial' & Var1 != 'Epithelial')
ggplot(fibro.epi,aes(x=reorder(Var1,value,median),y=value,fill=Var1))+geom_boxplot(outlier.size = 0.1)+
  facet_wrap(~From)+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.cl)

fibro.B = filter(nei.matrix,Var2 == 'Fibroblast' & Var1 != 'Fibroblast')
ggplot(fibro.B,aes(x=reorder(Var1,value,median),y=value,fill=Var1))+geom_boxplot(outlier.size = 0.1)+
  facet_wrap(~From)+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.cl)


myeloid.ls= filter(nei.matrix,Var2 == 'Myeloid' & Var1 != 'Myeloid')

ggplot(myeloid.ls[myeloid.ls$Var1 == 'LS-like',],aes(x=reorder(From,value,median),y=value,fill=From))+geom_boxplot(outlier.size = 0.1)+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.From)
ggplot(myeloid.ls[myeloid.ls$From == 'mLN',],aes(x=reorder(Var1,value,median),y=value,fill=Var1))+geom_boxplot(outlier.size = 0.1)+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.cl)

ggplot(myeloid.ls[myeloid.ls$From == 'mLN+ PT',],aes(x=reorder(Var1,value,median),y=value,fill=Var1))+geom_boxplot(outlier.size = 0.1)+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.cl)
ggplot(myeloid.ls[myeloid.ls$From == 'mLN- PT',],aes(x=reorder(Var1,value,median),y=value,fill=Var1))+geom_boxplot(outlier.size = 0.1)+
  coord_flip()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.cl)


####10. plot####

SpatialDimPlot(st.m,
               cols = cols.cl,
               group.by = 'cluster.name',
               pt.size.factor = 1.3,
               stroke = 0,
               alpha = 0,image.alpha = 1,
               label.box = F,
               ncol = 4,
               combine = T)

####11. difference in LS####

st.m.ls = st.m[,st.m$cluster.name == 'LS-like']
Idents(st.m.ls) = st.m.ls$From
ex.markers = c("CXCL13", "HAVCR2", "PDCD1", "TIGIT", "LAG3", "CTLA4", "LAYN", "RBPJ", "VCAM1", "GZMB", "TOX", "MYO7A")
ef.markers = c('PRF1', 'IFNG', 'CCL4', 'HLA-DQA1', 'GZMK', 'GZMA', 'GZMH', 'CD44', 'DUSP2', 'KLRB1', 'KLRD1' , 'CTSW')
st.m.ls = AddModuleScore(st.m.ls,features = list(cc=ex.markers),name = 'exh_Score',nbin=20)
VlnPlot(st.m.ls,features = c('exh_Score1'),pt.size = 0,group.by = 'From')
st.m.ls.meta = st.m.ls@meta.data

ggplot(st.m.ls.meta,aes(x=From,fill=From,y=exh_Score1))+geom_violin(draw_quantiles = c(0.5),scale = 'width')+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
  theme_cowplot()+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90))
  

st.m.ls.t = st.m.ls[,st.m.ls$From != 'mLN']
Idents(st.m.ls.t)=st.m.ls.t$From
st.m.ls.diff = FindMarkers(st.m.ls.t,ident.1 = 'mLN+ PT',
                           logfc.threshold = 0.1,min.pct = 0)
st.m.ls.diff$gene = rownames(st.m.ls.diff)
write.csv(st.m.ls.diff,file = './tables_v2/st.m.ls.diff.csv')


c12.log2FC = st.m.ls.diff$avg_log2FC
names(c12.log2FC)=rownames(st.m.ls.diff)
c12.log2FC.s = sort(c12.log2FC,decreasing = T)
c12.gsea = GSEA(c12.log2FC.s,minGSSize = 5,TERM2GENE = hallmark.items,eps = 0,pvalueCutoff = 2)
c12.gsea.res.hallmark = c12.gsea@result
write.csv(c12.gsea.res.hallmark,file = './tables_v2/LS diff hallmark gsea.csv')

GC.markers = c('BCL6','SUGCT','SSBP2','STAG3','RGS13',
               'LCK','MYO1E','IRAG2','LAT2','RRAS2',
               'MARCKSL1','SERPINA9','PRPSAP2',
               'LMO2','CD22','AICDA')

st.m.ls = AddModuleScore(st.m.ls,features = list(GC= intersect(rownames(st.m.ls),GC.markers)),name = 'GC_Score',nbin = 18)

VlnPlot(st.m.ls,features = 'GC_Score1',group.by = 'From')
VlnPlot(st.m.ls,features = 'HAVCR2',group.by = 'From',adjust = 100)

st.meta = st.m.ls@meta.data

ggplot(st.meta,aes(x=From,fill=From,y=GC_Score1))+geom_boxplot(outlier.size = 0.1)+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(2,3),c(1,3)))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = cols.From)

#####Neighbors####

st.m.lym = st.m.ls


lym.neis = colnames(st.m.lym)

samples = unique(names(st.m@images))

for(sample in samples){
  st.m.sample = st.m[,st.m$orig.ident == sample]
  
  lym.i = colnames(st.m.sample[,st.m.sample$cluster.name == 'LS-like'])
  
  i.coord = st.m@images[[sample]]@coordinates
  i.coord$pairxy = paste(i.coord$row,i.coord$col)
  
  for(spot.i in lym.i){
    row.i = as.integer(i.coord[spot.i,'row'])
    col.i = as.integer(i.coord[spot.i,'col'])
    candi.coord.i=data.frame(rowc=c(row.i,row.i,row.i-1,row.i-1,row.i+1,row.i+1),
                             colc=c(col.i-2,col.i+2,col.i-1,col.i+1,col.i-1,col.i+1))
    candi.coord.i$pairxy = paste(candi.coord.i$rowc,candi.coord.i$colc)
    spots.info.c.i = i.coord[i.coord$pairxy %in% as.character(candi.coord.i$pairxy),]
    
    lym.neis = union(lym.neis,rownames(spots.info.c.i))
    
  }
  
  print(sample)
}

st.m.lym.neis = st.m[,lym.neis]
#st.m.lym.Epi = st.m[,st.m$cluster.name %in% c('LS-like','Epithelial')]


for(sample in samples){
  st.m.sample = st.m.lym.neis[,st.m.lym.neis$orig.ident == sample]
  
  lym.i = colnames(st.m.sample)
  
  i.coord = st.m@images[[sample]]@coordinates
  i.coord$pairxy = paste(i.coord$row,i.coord$col)
  
  for(spot.i in lym.i){
    row.i = as.integer(i.coord[spot.i,'row'])
    col.i = as.integer(i.coord[spot.i,'col'])
    candi.coord.i=data.frame(rowc=c(row.i,row.i,row.i-1,row.i-1,row.i+1,row.i+1),
                             colc=c(col.i-2,col.i+2,col.i-1,col.i+1,col.i-1,col.i+1))
    candi.coord.i$pairxy = paste(candi.coord.i$rowc,candi.coord.i$colc)
    spots.info.c.i = i.coord[i.coord$pairxy %in% as.character(candi.coord.i$pairxy),]
    
    lym.neis = union(lym.neis,rownames(spots.info.c.i))
    
  }
  
  print(sample)
}
st.m.lym.neis.2 = st.m[,lym.neis]


SpatialDimPlot(st.m.lym.neis.2,
               cols = cols.cl,
               group.by = 'cluster.name',
               pt.size.factor = 1.5,
               stroke = 0,
               alpha = 1,image.alpha =1,
               label.box = F,
               ncol =8,
               combine = T)


SpatialDimPlot(st.m.lym.neis.2,
               cols = cols.cl,images = c('P7L','P5L'),
               group.by = 'cluster.name',
               pt.size.factor = 1.5,
               stroke = 0,
               alpha = 1,image.alpha =1,
               label.box = F,
               ncol =2,
               combine = T)
saveRDS(st.m.lym.neis.2,file = './variables_v2/st.m.LS.neis.twoStep.Rds')
st.nei.meta = st.m.lym.neis.2@meta.data

st.nei.meta.sum = dplyr::group_by(st.nei.meta,orig.ident,From,cluster.name) %>% dplyr::summarise(n=n())
colnames(st.nei.meta.sum)[3] = 'cellType'


st.nei.meta.sum.2 = st.nei.meta.sum
st.nei.meta.sum.2.sum = dplyr::group_by(st.nei.meta.sum.2,orig.ident) %>% dplyr::summarise(n=sum(n))
st.nei.meta.sum.2.sum = as.data.frame(st.nei.meta.sum.2.sum)
rownames(st.nei.meta.sum.2.sum)=st.nei.meta.sum.2.sum$orig.ident

st.nei.meta.sum.2$score = unlist(lapply(1:nrow(st.nei.meta.sum.2),function(i){
  
  a= as.double(st.nei.meta.sum[i,'n']/st.nei.meta.sum.2.sum[as.character(st.nei.meta.sum[i,'orig.ident']),'n'])
})
)
ggplot(st.nei.meta.sum,aes(y=cellType,fill=cellType,x=n))+geom_boxplot(outlier.size = 0.1)+
  facet_wrap(~From)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = cols.cl)


ggplot(st.nei.meta.sum.2,aes(y=cellType,fill=cellType,x=score))+geom_boxplot(outlier.size = 0.1)+
  facet_wrap(~From)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = cols.cl)

st.nei.meta.sum.2$From = factor(st.nei.meta.sum.2$From,levels = c('mLN','mLN+ PT','mLN- PT'))
ggplot(st.nei.meta.sum.2,aes(x=From,fill=From,y=score))+geom_boxplot(outlier.size = 0.1)+
  facet_wrap(~cellType)+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = cols.From)


ggplot(st.nei.meta.sum.2[st.nei.meta.sum.2$cellType == 'Epithelial',],aes(x=From,fill=From,y=score))+geom_boxplot(outlier.size = 0.1)+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = cols.From)

write.csv(st.nei.meta.sum.2,file = './tables_v2/st.nei.meta.sum.2_scores.csv')

Idents(st.m.lym.neis.2)=st.m.lym.neis.2$cluster.name

marker.neis.2 = FindAllMarkers(st.m.lym.neis.2)

st.only.neis = st.m.lym.neis.2[,st.m.lym.neis.2$cluster.name != 'LS-like']

Idents(st.only.neis)=st.only.neis$From

st.only.neis.diff = FindMarkers(st.only.neis,ident.1 = 'mLN+ PT',ident.2 = 'mLN- PT')

SpatialFeaturePlot(st.m.lym.neis.2,
                   features = c('C3','FTL','IGHG4','VIM','APOC1'),
              images = c('P5T','P11_T1'),
               pt.size.factor = 1.5,
               stroke = 0,
               alpha = 1,image.alpha =1,
               ncol =4,
               combine = T)

st.m.mye = st.m[,st.m$cluster.name == 'Myeloid']
st.m.mye$neiType = ifelse(colnames(st.m.mye) %in% colnames(st.m.lym.neis.2),'Nei_LS','Others')

Idents(st.m.mye)=st.m.mye$neiType
diff.mye.nei = FindMarkers(st.m.mye,ident.1 = 'Nei_LS')
