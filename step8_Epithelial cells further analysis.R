####Further analysis about epithelial cells, 2024-09-23####
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
library(clusterProfiler)
library(GSVA)
library(survival)
library(survminer)
generate_subColors = function(cols.cell,sub.labels){
  
  cols.generate = c()
  major.labels = names(cols.cell)
  sub.labels.major =  as.character(unlist(lapply(as.character(sub.labels),function(a){
    if(grepl('_',a)){
      return(strsplit(a,'_')[[1]][1])
    }else{
      return(a)
    }
  })))
  
  cols = character(length = length(sub.labels))
  
  for(label in major.labels){
    
    sub.i.ls = sub.labels[sub.labels.major == label]
    
    cols.i = colorRampPalette(c('white',cols.cell[label]))(1+length(sub.i.ls))
    
    cols[sub.labels.major == label]=cols.i[-1]
    
  }
  
  return(cols)
}
cols.default = c("#93bcd0","#ebad71","#bb8f7c","#639f95","#a0ca8a","#479544","#e58b8e","#d65559","#32709f",
                 '#b8a3c4','#60428d','#9f5633','grey','#c92e2a','#c79795','#d9c17c',"#d0896b",'#5d94b4',"black","blue",'cyan')

epi.data = readRDS(file = './variables_v2/epi.data.Rds')
epi.data = epi.data[,epi.data[['subCellType']] != 'Remove']

####8.1 plot the basic information####
dim(epi.data)#14787 181421
cols.cell = cols.default[-2][1:length(unique(SummarizedExperiment::colData(epi.data)$subCellType))]
names(cols.cell)=unique(SummarizedExperiment::colData(epi.data)$subCellType)[order(unique(SummarizedExperiment::colData(epi.data)$subCellType))]
# plot_cells(epi.data,group_cells_by = 'subCellType',color_cells_by = 'subCellType',
#            cell_stroke = 0,cell_size = 0.4,label_cell_groups = F,
#            rasterize = T,
#            group_label_size = 3)+scale_color_manual(values = cols.cell )+theme_void()
scater::plotReducedDim(epi.data,dimred = 'UMAP',colour_by = 'subCellType',point_size = 1)+
  ggplot2::scale_color_manual(values = cols.cell)+
  theme_void()

cols.clusters = generate_subColors(cols.cell ,unique(epi.data[['assigned_cell_type']]))
names(cols.clusters)= as.character(unique(epi.data[['assigned_cell_type']]))

scater::plotReducedDim(epi.data,dimred = 'UMAP',colour_by = 'assigned_cell_type',point_size = 1)+
  ggplot2::scale_color_manual(values = cols.clusters)+
  theme_void()

other.marker = c("AGER","RTKN2","CLIC5",#AT1
                 "SFTPC","SFTPA1",'SFTPA2','WFDC2','SFTPB',#AT2
                 'GSTA1','DTHD1','DCDC2B','PIFO','FOXJ1','TPPP3','CAPS',#Ciliated
                 'PIP','SUCNR1',#PIP
                 'SCGB1A1','SCGB3A2', 'SCGB3A1',#Club
                 'KRT5','TP63','LY6D','KRT6A','KRT15'#Basal
                 
)
plot_genes_by_group(epi.data,other.marker,
                    group_cells_by="subCellType",
                    ordering_type="none",
                    max.size=4)+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colours =rev(hcl.colors(100,'Spectral')))+
  theme(axis.text.x = element_text(angle = 90))

epi.PIP = epi.data[,epi.data[['subCellType']] == "PIP+ Epi"]
scater::plotReducedDim(epi.PIP,dimred = 'UMAP',colour_by = 'PatientID',point_size = 1)+
  ggplot2::scale_color_manual(values = c(as.character(jcolors('pal2')),cols.default))+
  theme_void()


epi.s = epi.data[,epi.data[['subCellType']] == "SUCNR1+ Epi"]
scater::plotReducedDim(epi.s,dimred = 'UMAP',colour_by = 'PatientID',point_size = 1)+
  ggplot2::scale_color_manual(values = c(as.character(jcolors('pal2')),cols.default))+
  theme_void()

####scater::plot no point bounds 

cols.diseases = jcolors()[2:5]
names(cols.diseases)=unique(epi.data[['Disease']])
plot_cells(epi.data, color_cells_by="Disease", label_cell_groups=FALSE,cell_size = 0.1,cell_stroke = 0)+
  facet_wrap(~Disease,ncol = 4)+
  scale_color_manual(values = cols.diseases)+
  theme(legend.position = 'none')+
  theme_void()
scater::plotReducedDim(epi.data,dimred = 'UMAP',colour_by = 'Disease',point_size = 0.1)+
  ggplot2::scale_color_manual(values = cols.diseases)+
  theme_void()

cols.From = as.character(pal_npg()(6))
names(cols.From)= c('mLN+ N','mLN+ PT','mLN','mLN- N','mLN- PT','nLN')



table(epi.data[['sorting']],epi.data[['cluster.epi']])

table(epi.data[['orig.ident']],epi.data[['cluster.epi']])

table(epi.data[['From']],epi.data[['cluster.epi']])

epi.data[['From.2']] = ifelse(epi.data[['From']] %in% c('mLN','mLN+ PT'), 'mLN+','mLN-')
marker_test_res.From.2 = top_markers(epi.data, group_cells_by="From.2", genes_to_test_per_group=200,
                                     reference_cells=1000, cores=32)

####8.2 statistics on cell types####
meta.data.show = epi.data@colData
meta.data.show = as.data.frame(meta.data.show)
meta.data.show$Disease2 = substr(meta.data.show$Disease,1,4)

ggplot2::ggplot(meta.data.show[meta.data.show$subCellType %in% c("SUCNR1+ Epi","PIP+ Epi") == F,],ggplot2::aes(x=From,fill=subCellType))+
  ggplot2::geom_bar(position = 'fill',width = 0.6)+facet_wrap(~Disease2,ncol = 1)+
  theme_cowplot(font_size = 8)+
  ggplot2::scale_fill_manual(values =  cols.cell)+
  theme(axis.text  = element_text(size = 8),
        axis.text.x = element_text(angle = 90,hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))+
  ylab('Proportion')
table(epi.data[['orig.ident']],epi.data[['cluster.epi']])[,'Epi_C12']

com2orig = table(epi.data[['From']],epi.data[['assigned_cell_type']])

f.data<-as.matrix(com2orig)
#f.data = f.data/apply(f.data, 1, sum)
#f.data = t(f.data)
#rownames(ids.all)=ids.all$id
meta.data.u = meta.data.show
meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
rownames(meta.data.u)=meta.data.u$orig.ident

com2orig.df = melt(f.data,measure.vars = colnames(f.data))
com2orig.df = com2orig.df[com2orig.df$Var2 != 'Unassigned',]##zero values

com2orig.df$Var2 = as.character(com2orig.df$Var2)
com2orig.df$Var2 = factor(com2orig.df$Var2,levels = c('AT1',paste0('AT2_C',1:8),paste0('Basal_C',1:4),
                                                      'Ciliated','Club','PIP+ Epi','SUCNR1+ Epi' ))



ggplot(com2orig.df,aes(x='',y=value,fill=Var1,color=Var1))+
  geom_bar(stat="identity", width=1, position = 'fill')+coord_polar("y", start=0)+
  facet_wrap(~Var2,ncol = 6)+
  theme_void()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 0))+
  scale_fill_manual(values = alpha(cols.From,0.7))+
  scale_color_manual(values = alpha(cols.From,0.7))

####8.3 Difference between mLN+ and mLN- epithelial cells####

epi.data[['From.2']] = ifelse(epi.data[['From']] %in% c('mLN','mLN+ PT'), 'mLN+','mLN-')

epi.data.AT2 = epi.data[,epi.data[['subCellType']] == 'AT2']
marker_test_res.AT2.From2 = top_markers(epi.data.AT2, group_cells_by="From.2", genes_to_test_per_group=100,
                                     reference_cells=1000, cores=32)
marker_diff.AT2.mLN = marker_test_res.AT2.From2[marker_test_res.AT2.From2$cell_group == 'mLN+',]
rownames(marker_diff.AT2.mLN)=marker_diff.AT2.mLN$gene_id
marker_diff.AT2.mLN = marker_diff.AT2.mLN[marker_diff.AT2.mLN$marker_test_p_value<0.01,]
marker_diff.AT2.mLN = marker_diff.AT2.mLN[marker_diff.AT2.mLN$gene_id != 'MALAT1' & substr(marker_diff.AT2.mLN$gene_id,1,2) %in% c('MT','RP') == F,]
write.csv(marker_diff.AT2.mLN,file = './tables_v2//AT2 mLN high markers.csv')


epi.data.Basal = epi.data[,epi.data[['subCellType']] == 'Basal']
marker_test_res.Basal.From2 = top_markers(epi.data.Basal, group_cells_by="From.2", genes_to_test_per_group=100,
                                        reference_cells=1000, cores=32)

marker_diff.Basal.mLN = marker_test_res.Basal.From2[marker_test_res.Basal.From2$cell_group == 'mLN+',]
marker_diff.Basal.mLN = marker_diff.Basal.mLN[marker_diff.Basal.mLN$marker_test_p_value <0.01,]
marker_diff.Basal.mLN = marker_diff.Basal.mLN[marker_diff.Basal.mLN$gene_id != 'MALAT1' & substr(marker_diff.Basal.mLN$gene_id,1,2) %in% c('MT','RP') == F,]
rownames(marker_diff.Basal.mLN)=marker_diff.Basal.mLN$gene_id

write.csv(marker_diff.Basal.mLN,file = './tables_v2//Basal mLN high markers.csv')

top.20.union = union(top_n(marker_diff.AT2.mLN,20,mean_expression)$gene_id,top_n(marker_diff.Basal.mLN,20,mean_expression)$gene_id)

top.20.AT2 = top_n(marker_diff.AT2.mLN,20,mean_expression)
top.20.AT2 = arrange(top.20.AT2,mean_expression)
top.20.AT2$gene_id = factor(top.20.AT2$gene_id,levels = top.20.AT2$gene_id)
top.20.AT2$CellType = 'AT2'

top.20.Basal = top_n(marker_diff.Basal.mLN,20,mean_expression)
top.20.Basal = arrange(top.20.Basal,mean_expression)
top.20.Basal$gene_id = factor(top.20.Basal$gene_id,levels = top.20.Basal$gene_id)
top.20.Basal$CellType = 'Basal'

top.res = rbind(top.20.AT2[,colnames(top.20.Basal)],top.20.Basal)
write.csv(top.res,file = './tables_v2//AT2 Basal PLM specific genes_top.csv')
top.res = read.csv(file = './tables_v2//AT2 Basal PLM specific genes_top.csv',row.names = 1)
top.res = top.res[order(top.res$mean_expression),]
top.res$gene_id = factor(top.res$gene_id,levels = top.res$gene_id[!duplicated(top.res$gene_id)])
ggplot(top.res,aes(x=mean_expression,y=gene_id,fill=CellType))+geom_bar(stat = 'identity',width = 0.7)+facet_grid(~CellType)+
  theme_cowplot(font_size = 8)+scale_fill_manual(values = cols.cell[c('AT2','Basal')])+
  theme(text = element_text(size = 8),
        axis.text = element_text(size = 8),
        strip.background = element_blank())

plot_genes_violin(epi.data.Basal[c('IGHG1','IGHG4','IGHA1','IGLC2'),] ,
           group_cells_by = 'From')+ scale_fill_manual(values = cols.From)

plot_genes_violin(epi.data[c('FTH1'),] ,normalize = F,
                  group_cells_by = 'From')+ scale_fill_manual(values = cols.From)

rm(epi.data.AT2)
rm(epi.data.Basal)
gc()

plot_genes_violin(epi.data[c('IGHG1','IGHG4','IGHA1','IGLC2'),] ,
                  group_cells_by = 'From')+ 
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())
####8.4 CytoTRACE analysis ####

#####LUSC#####

library(CytoTRACE)



minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}


epi.data.LUSC = epi.data[,epi.data[['Disease']] %in% c('LUSC','LUSC_Normal')]
expr.LUSC = as.matrix(epi.data.LUSC@assays@data$counts)

# epi.data.LUSC = preprocess_cds(epi.data.LUSC)
# epi.data.LUSC = reduce_dimension(epi.data.LUSC,cores = 20)
# 
selected_genes <- rownames(expr.LUSC)[ Matrix::rowSums(epi.data.LUSC@assays@data$counts > 0)> 200]
cyto.LUSC <- CytoTRACE(expr.LUSC[selected_genes,],batch = as.character(colData(epi.data.LUSC)$Source),ncores = 12)
#plotCytoGenes(cyto.LUSC, numOfGenes = 10)

gc()

cyto.genes = cyto.LUSC$cytoGenes


sub.cells = sample(colnames(epi.data.LUSC),500)
cyto.genes.temp = cyto.genes[substr(names(cyto.genes),1,2) %in% c('MT','RP')==F]
top.genes = c(cyto.genes.temp[1:25],rev(cyto.genes.temp)[1:25])
plot.mat = log2( as.data.frame(t(expr.LUSC[names(top.genes),sub.cells]))+1)

plot.mat$CytoTRACE = as.numeric(cyto.LUSC[['CytoTRACE']][sub.cells])
plot.mat = plot.mat[order(plot.mat$CytoTRACE),]

pheatmap::pheatmap(t(as.data.frame(lapply(plot.mat,minMax)))[c('CytoTRACE',names(top.genes)),],color = rev(hcl.colors(100,'Spectral')),show_colnames = F,cluster_cols = F,
                   cluster_rows = F,gaps_row = c(1,26),fontsize = 8)


epi.data.LUSC[['CytoTRACE']] = cyto.LUSC[['CytoTRACE']]
epi.data.LUSC[['CytoTRACE_Level']] = ifelse(epi.data.LUSC[['CytoTRACE']] > median(epi.data.LUSC[['CytoTRACE']]),'High_CytoTRACE','Low_CytoTRACE')



meta.data.epi.LUSC = as.data.frame(colData(epi.data.LUSC))
meta.data.epi.LUSC$CD44 = as.numeric(assay(epi.data.LUSC)['CD44',])

meta.data.epi.LUSC$From = factor(meta.data.epi.LUSC$From,
                                     levels=c('mLN+ N','mLN+ PT','mLN','mLN- PT','mLN- N'))
ggplot(meta.data.epi.LUSC,aes(y=CytoTRACE,x=From,fill=From,color=From))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+scale_fill_manual(values = alpha(cols.From,0.7))+
  scale_color_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle=90))

ggplot(meta.data.epi.LUSC,aes(y=log2(CD44+1),x=CytoTRACE_Level,fill=CytoTRACE_Level))+geom_boxplot(outlier.size = 0.2)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle=90))

saveRDS(cyto.LUSC,file = './variables_v2/cytoTRACE_LUSC.Rds')
#saveRDS(epi.data.LUSC,file = './variables_v2/epi.data.LUSC.Rds')
cyto.LUSC = readRDS(file = './variables_v2/cytoTRACE_LUSC.Rds')
rm(expr.LUSC)
gc()
#####LUAD#####

epi.data.LUAD = epi.data[,epi.data[['Disease']] %in% c('LUAD','LUAD_Normal')]
expr.LUAD = as.matrix(epi.data.LUAD@assays@data$counts)
selected_genes <- rownames(expr.LUAD)[ Matrix::rowSums(epi.data.LUAD@assays@data$counts > 0)> 200]
set.seed(11223)
cyto.LUAD <- CytoTRACE(expr.LUAD[selected_genes,],ncores = 12,batch = as.character(colData(epi.data.LUAD)$Source))

epi.data.LUAD[['CytoTRACE']] = cyto.LUAD[['CytoTRACE']]
epi.data.LUAD[['CytoTRACE_Level']] = ifelse(epi.data.LUAD[['CytoTRACE']] > median(epi.data.LUAD[['CytoTRACE']]),'High_CytoTRACE','Low_CytoTRACE')

meta.data.epi.LUAD = as.data.frame(colData(epi.data.LUAD))

meta.data.epi.LUAD$From = factor(meta.data.epi.LUAD$From,
                                     levels=c('mLN+ N','mLN+ PT','mLN','mLN- PT','mLN- N'))
meta.data.epi.LUAD$CD44 = as.numeric(assay(epi.data.LUAD)['CD44',])
ggplot(meta.data.epi.LUAD,aes(y=CytoTRACE,x=From,fill=From,color=From))+geom_boxplot(outlier.size = 0.2)+
  theme_cowplot(font_size = 8)+scale_fill_manual(values = alpha(cols.From,0.7))+
  scale_color_manual(values = cols.From)+
  theme(axis.text.x = element_text(angle=90))

ggplot(meta.data.epi.LUAD,aes(y=log2(CD44+1),x=CytoTRACE_Level,fill=CytoTRACE_Level))+geom_boxplot(outlier.size = 0.2)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle=90))

cyto.genes = cyto.LUAD$cytoGenes
sub.cells = sample(colnames(expr.LUAD),500)
cyto.genes.temp = cyto.genes[substr(names(cyto.genes),1,2) %in% c('MT','RP')==F]
top.genes = c(cyto.genes.temp[1:25],rev(cyto.genes.temp)[1:25])
plot.mat = log2( as.data.frame(t(expr.LUAD[names(top.genes),sub.cells]))+1)

plot.mat$CytoTRACE = as.numeric(cyto.LUAD[['CytoTRACE']][sub.cells])
plot.mat = plot.mat[order(plot.mat$CytoTRACE),]
plot.mat.new = t(as.data.frame(lapply(plot.mat,minMax)))
#rownames(plot.mat.new)[24] = "C1QTNF3-AMACR"
pheatmap::pheatmap(plot.mat.new[c('CytoTRACE',names(top.genes)),],color = rev(hcl.colors(100,'Spectral')),show_colnames = F,cluster_cols = F,
                   cluster_rows = F,gaps_row = c(1,26),fontsize = 8)

saveRDS(cyto.LUAD,file = './variables/cytoTRACE_LUAD.Rds')
cyto.LUAD = readRDS(file = './variables/cytoTRACE_LUAD.Rds')
#####Pathway analysis#####
load(file = './path.info.RData')
path2gene.all = read.gmt(gmtfile = './MSigDB_Hallmark.gmt')

rank.genes = data.frame(gene = names(cyto.genes),scores = as.numeric(cyto.genes))
gsea.res = GSEA(cyto.genes,minGSSize = 3,seed=1234,TERM2GENE = path2gene.all)
gsea.res.LUSC = gsea.res@result
gsea.res.LUAD = gsea.res@result
ridgeplot(gsea.res,showCategory = 30)

gsea.res.LUSC$Disease = 'LUSC'
gsea.res.LUAD$Disease = 'LUAD'

path.res = rbind(gsea.res.LUSC,
                 gsea.res.LUAD)
res.filter = filter(path.res,pvalue<0.01) %>% group_by(Disease) 
res.filter$Description = factor(res.filter$Description,levels = rev(as.character(union(gsea.res.LUAD$Description,gsea.res.LUSC$Description))))
ggplot(res.filter,aes(x=-log10(p.adjust),y=Description,fill=NES))+geom_bar(stat = 'identity',width = 0.65)+
  facet_grid(~Disease,scales = 'free',space = 'free')+geom_vline(xintercept = -log10(0.05),lty=3)+
  geom_vline(xintercept = 0,lty=1)+
  theme_cowplot(font_size = 8)+
  scale_fill_gradientn(colours = rev(hcl.colors(100,'Spectral')))+
  theme(axis.text = element_text(size=8),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 8))
write.csv(path.res,file = './tables_v2/Epi cytoTRACE gsea res.csv')
path.res = read.csv(file = './tables_v2/Epi cytoTRACE gsea res.csv',row.names = 1)

cyto.genes.luad = names(cyto.LUAD$cytoGenes[abs(cyto.LUAD$cytoGenes)>0.3])
cyto.genes.lusc = names(cyto.LUSC$cytoGenes[abs(cyto.LUSC$cytoGenes)>0.3])


both.cyto.genes = intersect(cyto.genes.luad,cyto.genes.lusc)

both.cyto.path = enricher(both.cyto.genes,pvalueCutoff = 1.2,TERM2GENE = path2gene.all)@result

#####shared genes#####

cyto.genes.luad = names(cyto.LUAD$cytoGenes[abs(cyto.LUAD$cytoGenes)>0.5])
cyto.genes.lusc = names(cyto.LUSC$cytoGenes[abs(cyto.LUSC$cytoGenes)>0.5])


both.cyto.genes = intersect(cyto.genes.luad,cyto.genes.lusc)

plot.both.data = data.frame(gene = both.cyto.genes,
                            LUAD_Corr = as.numeric(cyto.LUAD$cytoGenes[both.cyto.genes]),
                            LUSC_Corr = as.numeric(cyto.LUSC$cytoGenes[both.cyto.genes]))


plot.both.data$sum = plot.both.data$LUAD_Corr*plot.both.data$LUSC_Corr

top.both = union(names(cyto.LUAD$cytoGenes[cyto.LUAD$cytoGenes>0.6]),
                 names(cyto.LUSC$cytoGenes[cyto.LUSC$cytoGenes>0.6]))
plot.both.data$label = ''
plot.both.data$label[plot.both.data$gene %in% top.both] = plot.both.data$gene[plot.both.data$gene %in% top.both]

ggplot(plot.both.data,aes(x=LUSC_Corr,y=LUAD_Corr,color = sum))+geom_point()+
  ggrepel::geom_text_repel(aes(label = label),size=3,max.overlaps = 5000)+
  theme_cowplot(font_size = 8)+
  scale_color_gradientn(colors=hcl.colors(100))
write.csv(plot.both.data,file = './table results/LUAD and LUSC cytoTRACE corr both.csv')

#####Plot cytoTRACE scores in umap#####

meta.data.epi = as.data.frame(colData(epi.data))
meta.data.epi$cytoTRACE = NA
meta.data.epi[names(cyto.LUAD$CytoTRACE),'cytoTRACE'] = cyto.LUAD$CytoTRACE
meta.data.epi[names(cyto.LUSC$CytoTRACE),'cytoTRACE'] = cyto.LUSC$CytoTRACE
epi.data[['CytoTRACE']] = meta.data.epi$cytoTRACE


scater::plotReducedDim(epi.data,dimred = 'UMAP',colour_by = 'cytoTRACE',point_size = 0.2)+
  ggplot2::scale_color_manual(values = cols.cell)+
  scale_color_gradientn(colours=rev(hcl.colors(100,'Spectral')))+
  theme_void()

scater::plotReducedDim(epi.data.LUSC,dimred = 'UMAP',colour_by = 'CytoTRACE',point_size = 0.2)+
  ggplot2::scale_color_manual(values = cols.cell)+
  scale_color_gradientn(colours=rev(hcl.colors(100,'Spectral')))+
  theme_void()+
  theme(legend.position = 'none')
scater::plotReducedDim(epi.data.LUAD,dimred = 'UMAP',colour_by = 'CytoTRACE',point_size = 0.2)+
  ggplot2::scale_color_manual(values = cols.cell)+
  scale_color_gradientn(colours=rev(hcl.colors(100,'Spectral')))+
  theme_void()+
  theme(legend.position = 'none')

#####NicheNet#####
  
####8.5 survival analysis####

tumor='BRCA'
mat<-read.delim(file=paste0("./data/",tumor,'ExpressionNormalized.txt'),sep='\t',check.names = FALSE)
clin<-read.delim(file=paste0('./data/',tumor,'_clin.txt'),sep='\t')
notDead <- is.na(clin$days_to_death)
if (length(notDead) > 0) {
  clin[notDead,]$days_to_death <-
    clin[notDead,]$days_to_last_follow_up
}
rownames(clin)=clin$bcr_patient_barcode
samplesTP<-colnames(mat)[substr(colnames(mat),14,15)=='01']

#library(limma)

mat.tp = log2(mat[,samplesTP]+1)

colnames(mat.tp) = substr(colnames(mat.tp),1,12)

cyto.genes.luad = names(cyto.LUAD$cytoGenes[abs(cyto.LUAD$cytoGenes)>0.5])
cyto.genes.lusc = names(cyto.LUSC$cytoGenes[abs(cyto.LUSC$cytoGenes)>0.5])


both.cyto.genes = intersect(cyto.genes.luad,cyto.genes.lusc)#46 genes
AA_Genes = intersect(both.cyto.genes,rownames(mat.tp))
clin.cls = intersect(clin$bcr_patient_barcode,colnames(mat.tp))

clin.av = clin[clin.cls,]
clin.av$AA_PS = apply(mat.tp[AA_Genes,clin.cls],2,mean)
clin.av$AA_PS_Level = ifelse(clin.av$AA_PS > median(clin.av$AA_PS),'H','L')

# AA_gsva = gsva(as.matrix(mat.tp[,clin.cls]),list('CSC_Score'=AA_Genes),parallel.sz=12)
# clin.av$AA_GSVA = as.double(AA_gsva[,clin.cls])


#clin.av$stage = gsub('a|b','',clin.av$tumor_stage)
# ggplot(clin.av,aes(x=stage,y=AA_PS,fill=stage))+geom_violin(draw_quantiles = 0.5,scale = 'width')+
#   stat_compare_means(comparisons = list(c(2,3),c(2,4),c(2,5)))+
#   theme_cowplot()+
#   theme(axis.text.x = element_text(angle = 90,hjust = 1))
# ggplot(clin.av,aes(x=stage,y=AA_GSVA,fill=stage))+geom_violin(draw_quantiles = 0.5,scale = 'width')+
#   stat_compare_means(comparisons = list(c(2,3),c(2,4),c(2,5)))+
#   theme_cowplot()+
#   theme(axis.text.x = element_text(angle = 90,hjust = 1))
#clin.av$AA_GSVA_Level = ifelse(clin.av$AA_GSVA > median(clin.av$AA_GSVA),'H','L')

clin.av$status <- grepl("dead", clin.av$vital_status, ignore.case = TRUE)

clin.av$days_to_death = as.double(clin.av$days_to_death)

coxph(Surv(days_to_death,status) ~ AA_PS,clin.av)

clin.av$AA_PS_Level = ifelse(clin.av$AA_PS > quantile(clin.av$AA_PS,0.5),'H','L')
ggsurvplot(survfit(Surv(days_to_death,status) ~ AA_PS_Level,clin.av),data = clin.av,pval = T,
           conf.int = T,
           risk.table = T,
           palette  = brewer.pal(3,'Set1')[c(1,2)]
)

write.csv(clin.av,file='./tables_v2/CSC clin kirc.csv')

clin.av = read.csv(file='./tables_v2/CSC clin luad.csv',row.names = 1)

write.csv(meta.data.epi,file = './tables_v2/epi.data metadata.csv')


####8.6 IgG relevant difference####

IgG.plot1 = plot_genes_violin(epi.data[c('IGHG1'),],group_cells_by = 'From')
IgG.plot2 = plot_genes_violin(epi.data[c('IGLC2'),],group_cells_by = 'From')
IgG.plot3 = plot_genes_violin(epi.data[c('IGHA1'),],group_cells_by = 'From')
IgG.plot4 = plot_genes_violin(epi.data[c('IGHG4'),],group_cells_by = 'From')


epi.data.meta = as.data.frame(colData(epi.data))

epi.data.meta.2 = epi.data.meta
epi.data.meta.2$Disease2 = substr(epi.data.meta.2$Disease,1,4)
#epi.data.meta.2$CytoTRACE = epi.data[['cytoTRACE']]
#epi.data.meta.2$IGHG4 = IGHG4.plot$data$expression
epi.data.meta.2$IGHG4 = IgG.plot4$data$expression
epi.data.meta.2$IGHG1 = IgG.plot1$data$expression
epi.data.meta.2$IGLC2 = IgG.plot2$data$expression
epi.data.meta.2$IGHA1 = IgG.plot3$data$expression
epi.data.meta.2$IgG.score = apply(epi.data.meta.2[,c('IGHG4','IGHG1','IGLC2','IGHA1')],1,median)
epi.data.meta.2$IgG_Level = ifelse(epi.data.meta.2$IgG.score  > median(epi.data.meta.2$IgG.score ),'High_Ig','Low_Ig')
pathway.p=c()
pathway.fc = c()

pathways = colnames(cds.gsva.2@gsva)

saveRDS(epi.data.meta.2,file='./variables_v2/epi.data.meta.2.Rds')

for(pathway in pathways){
  
  test.pathway = kruskal.test(as.formula(paste0(pathway,'~From')),epi.data.meta.2)
  pathway.p = append(pathway.p,test.pathway$p.value)
  
}
res.pathway = data.frame(pathways,pathway.p)

ggplot(epi.data.meta.2,aes(x=From,y=log2(IgG.score+1),fill=From))+
  geom_violin(draw_quantiles = 0.5,scale = 'width')+
  facet_wrap(~Disease2)+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.From)
ggplot(epi.data.meta.2,aes(x=Stage,y=log2(IgG.score+1),fill=Stage))+
  geom_violin(draw_quantiles = 0.5,scale = 'width')+
  facet_wrap(~Source)+
  stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
  theme_cowplot(font_size = 8)

ggplot(epi.data.meta.2,aes(x=cluster.epi,y=APOPTOSIS,fill=cluster.epi))+geom_violin(draw_quantiles = 0.5,scale = 'width')+
  theme_cowplot(font_size = 8)

ggplot(epi.data.meta.2,aes(x=cluster.epi,fill=From))+geom_bar(position = 'fill')+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.From)

ggplot(epi.data.meta.2,aes(x=From,y=EPITHELIAL_MESENCHYMAL_TRANSITION,fill=From))+geom_violin(draw_quantiles = 0.5,scale = 'width')+
  theme_cowplot(font_size = 8)+
  scale_fill_manual(values = cols.From)

ggplot(epi.data.meta.2,aes(x=IGHG4_Level,y=APOPTOSIS,fill=IGHG4_Level))+
  stat_compare_means(comparisons = list(c(1,2)))+
  geom_violin(draw_quantiles = 0.5,scale = 'width')+
  facet_wrap(From~Disease2)+
  theme_cowplot(font_size = 8)

ggplot(epi.data.meta.2,aes(x=From,fill=IGHG4_Level))+
  geom_bar(position = 'fill')+
  facet_wrap(~Disease2)+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

IGHG4.corrs = unlist(lapply(pathways, function(a){
  cor(epi.data.meta.2[epi.data.meta.2$From %in% c('mLN','mLN+ PT'),]$IGHG4, epi.data.meta.2[epi.data.meta.2$From %in% c('mLN','mLN+ PT'),a],method = 'spearman')
}))

IGHG4.corr.res = data.frame(pathways,IGHG4.corrs)

#####Malignant or not####
annotations.row = read.csv(file = './tables_v2/all inferCNV annotations of cells.csv',row.names = 1)
malignant.cells.PLM = rownames(annotations.row[annotations.row$CNV.Type == 'Malignant',])#95716

epi.data.PLM = epi.data[,epi.data[['Source']] == 'PLM']
epi.data.PLM[['Malignant']] = ifelse(colnames(epi.data.PLM) %in% malignant.cells.PLM,'Malignant','Non-malignant')
epi.data.PLM[['Disease2']] = substr(epi.data.PLM[['Disease']],1,4)
plot_genes_violin(epi.data.PLM[c('IGHA1','IGHG1','IGHG4','IGLC2'),],group_cells_by = 'Malignant',ncol=2)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_blank())

# epi.data.PLM[['IgG.score']] = epi.data.meta.2[colnames(epi.data.PLM),'IgG.score']
# epi.data.PLM[['IgG_Level']] = ifelse(log2(epi.data.PLM[['IgG.score']]+1) > median(log2(epi.data.PLM[['IgG.score']]+1)),'H','L') 
# epi.data.PLM[['Disease2']] = substr(epi.data.PLM[['Disease']],1,4)
# Ig.diff = top_markers(epi.data.PLM,group_cells_by = 'IgG_Level',genes_to_test_per_group = 200,
#                          reference_cells=1000, cores=32)
# 
# Ig.diff = read.csv(file = './tables_v2/IgGroup diff.csv')
# write.csv(Ig.diff,file = './tables_v2/IgGroup diff.csv')


ggplot(as.data.frame(colData(epi.data.PLM)),aes(x=IgG_Level,y=cytoTRACE,fill=IgG_Level))+
  ggpubr::stat_compare_means(comparisons = list(c(1,2)))+
  geom_violin(draw_quantiles = 0.5,scale = 'width')+
  facet_wrap(~Disease2)+
  theme_cowplot(font_size = 8)+
  theme(strip.background = element_blank())

epi.data.PLM.sc = CreateSeuratObject(counts = assay(epi.data.PLM))
epi.data.PLM.sc = AddMetaData(epi.data.PLM.sc,as.data.frame(colData(epi.data.PLM)))
epi.data.PLM.sc = AddModuleScore(epi.data.PLM.sc,features = list('IgG'=c('IGHG4','IGHG1','IGLC2','IGHA1')),name = 'IgG_Score')
epi.data.PLM.sc$IgG_Level = ifelse(epi.data.PLM.sc$IgG_Score1 > median(epi.data.PLM.sc$IgG_Score1),'H','L')
Idents(epi.data.PLM.sc)=epi.data.PLM.sc$IgG_Level
epi.data.PLM[['IgG_Level']] = epi.data.PLM.sc$IgG_Level
plot_genes_violin(epi.data.PLM[c('HLA-DRB1','CST3','HLA-DRA','SLPI'),],group_cells_by = 'IgG_Level',ncol = 4)
plot_genes_violin(epi.data.PLM[c('IFNGR1','BCL10','HMOX1','TAP1'),],group_cells_by = 'IgG_Level',ncol = 4,normalize = F)

VlnPlot(epi.data.PLM.sc,features = c('IFNGR1','BCL10','HMOX1','TAP1'),group.by = 'IgG_Level',pt.size = 0,adjust = 10,ncol = 4)

saveRDS(epi.data.PLM.sc,file = './variables_v2/epi.data.PLM.sc.Rds')

Ig.diff.2 = FindMarkers(epi.data.PLM.sc,ident.1 = 'H',logfc.threshold = 0,min.pct = 0)
write.csv(Ig.diff.2,file = './tables_v2/IgGroup diff_seurat version.csv')

library(clusterProfiler)
go.items = read.gmt(file('./c5.all.v2024.1.Hs.symbols.gmt'))
hallmark.items = read.gmt(file('./MSigDB_Hallmark.gmt'))
c12.log2FC = Ig.diff.2$avg_log2FC
names(c12.log2FC)=rownames(Ig.diff.2)
c12.log2FC.s = sort(c12.log2FC,decreasing = T)
c12.gsea = GSEA(c12.log2FC.s,minGSSize = 5,TERM2GENE = hallmark.items,eps = 0,pvalueCutoff = 2)
c12.gsea.res.hallmark = c12.gsea@result
gseaplot(c12.gsea,geneSetID = 'APOPTOSIS')
ridgeplot(c12.gsea)
write.csv(c12.gsea.res.hallmark,file = './tables_v2/GSEA_hallmark_IgGroup diff_seurat version.csv')
c12.gsea.res.hallmark = read.csv(file = './tables_v2/GSEA_hallmark_IgGroup diff_seurat version.csv',row.names = 1)
c12.gsea.res.hallmark$ID = factor(c12.gsea.res.hallmark$ID,levels = rev(c12.gsea.res.hallmark$ID))
ggplot(c12.gsea.res.hallmark[1:8,],aes(y=ID,x=NES,fill=-log10(p.adjust)))+geom_bar(stat = 'identity',width = 0.7)+
  theme_cowplot(font_size = 8)+
  scale_fill_gradientn(colours = hcl.colors(100))

c12.gsea.2 = GSEA(c12.log2FC.s,minGSSize = 5,TERM2GENE = go.items,eps = 0,pvalueCutoff = 2)
c12.gsea.res.go = c12.gsea.2@result
c12.gsea.res.gobp = c12.gsea.res.go[grepl('GOBP',c12.gsea.res.go$ID),] 
c12.gsea.res.gobp$Type = ifelse(c12.gsea.res.gobp$NES>0,'Up','Down')
c12.gsea.res.gobp.top = filter(c12.gsea.res.gobp,p.adjust<0.05) %>% group_by(Type) %>% top_n(10,abs(NES))
c12.gsea.res.gobp.top$ID = gsub('GOBP_','',c12.gsea.res.gobp.top$ID)
c12.gsea.res.gobp.top$ID =  factor(c12.gsea.res.gobp.top$ID,levels = c12.gsea.res.gobp.top$ID)

write.csv(c12.gsea.res.go,file = './tables_v2/GSEA_go_IgGroup diff_seurat version.csv')
path2cate = read.csv(file='/share/dichen/lungNodeM/KEGG hsa pathway list.csv')
load('/share/dichen/lungNodeM/path.info.RData')

c12.gsea.3 = GSEA(c12.log2FC.s,minGSSize = 5,TERM2GENE = path.info,eps = 0,pvalueCutoff = 2)
c12.gsea.res.kegg = c12.gsea.3@result
write.csv(c12.gsea.res.kegg,file = './tables_v2/GSEA_kegg_IgGroup diff_seurat version.csv')


# stat.value = c()
# pathway.p=c()
# 
# for(pathway in pathways){
#   
#   test.pathway = kruskal.test(as.formula(paste0(pathway,'~IgG_Level')),epi.data.meta.2[colnames(epi.data.PLM),])
#   pathway.p = append(pathway.p,test.pathway$p.value)
#   stat.value = append(stat.value,test.pathway$statistic)
# }
# 
# res.pathway = data.frame(pathways,pathway.p,stat.value)
# res.pathway$FDR = p.adjust(res.pathway$pathway.p)
# ggplot(epi.data.meta.2,aes(x=IgG_Level,y=ANGIOGENESIS,fill=IgG_Level))+
#   geom_violin(draw_quantiles = 0.5,scale = 'width')+
#   theme_cowplot(font_size = 8)
# 
# 
# PLM.infor = read.csv(file = './documents/Patient Info.csv')
# PLM.infor$Stage2 = gsub('A|B','',PLM.infor$Stage)
# rownames(PLM.infor)=PLM.infor$Patient_ID
# colData(epi.data.PLM)$Stage = PLM.infor[colData(epi.data.PLM)$PatientID,'Stage2']
# 
# 
# ggplot(as.data.frame(colData(epi.data.PLM)),aes(x=Stage,y=log2(IgG.score+1),fill=Stage))+
#   geom_violin(draw_quantiles = 0.5,scale = 'width')+
#   stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))+
#   theme_cowplot(font_size = 8)

marker_test_res.epi = read.csv(file = './tables_v2/Epi assigned cell type marker.csv')

#####IgG high clusters####
com2orig = table(as.character(epi.data[['From']]),epi.data[['cluster.epi']])

f.data<-as.matrix(com2orig)
#f.data = f.data/apply(f.data, 1, sum)
#f.data = t(f.data)
#rownames(ids.all)=ids.all$id
meta.data.u = meta.data.show
meta.data.u = meta.data.u[!duplicated(meta.data.u$orig.ident),]
rownames(meta.data.u)=meta.data.u$orig.ident

com2orig.df = melt(f.data,measure.vars = colnames(f.data))
#com2orig.df = com2orig.df[com2orig.df$Var2 != 'Unassigned',]##zero values

com2orig.df$Var2 = as.character(com2orig.df$Var2)


ggplot(com2orig.df[com2orig.df$Var2 %in% c('Epi_C5','Epi_C9','Epi_C11'),],aes(x='',y=value,fill=Var1,color=Var1))+
  geom_bar(stat="identity", width=1, position = 'fill')+coord_polar("y", start=0)+
  facet_wrap(~Var2,ncol = 8)+
  theme_void()+
  theme(axis.text = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 0))+
  scale_fill_manual(values = cols.From,0.7)+
  scale_color_manual(values = cols.From,0.7)

plot_genes_violin(epi.data[c('IGHG1','IGHG4','IGHA1','IGLC2'),] ,
                  group_cells_by = 'cluster.epi')+ scale_fill_manual(values = cols.default)+
  theme(axis.text.x = element_text(angle = 90))


epi.sub = epi.data[,epi.data[['cluster.epi']] %in% c('Epi_C5','Epi_C9','Epi_C11')]

plot_genes_violin(epi.sub[c('IGHG1','IGHG4','IGHA1','IGLC2'),] ,
                  group_cells_by = 'From')+ scale_fill_manual(values = cols.From)+
  facet_wrap(~cluster.epi)+
  theme(axis.text.x = element_text(angle = 90))


#####IGHG in ST####

st.m.epi = st.m[,st.m$cluster.name == 'Epithelial']
SpatialFeaturePlot(st.m.epi,features = 'IGHG4',pt.size.factor = 1.2,ncol=8)
VlnPlot(st.m.epi,features = 'IGHG4',pt.size = 0,group.by = 'From')


####Spatial analysis####

#####differences of Ig among three types #####
st.m = readRDS('./variables_v2/st.m_bayesSpace.Rds')

st.m.epi = st.m[,st.m$cluster.name == 'Epithelial']
cds.st.epi = new_cell_data_set(as(as.matrix(st.m.epi@assays$Spatial@counts), "sparseMatrix"),
                               cell_metadata = st.m.epi@meta.data,
                               data.frame(id = rownames(st.m.epi@assays$Spatial@counts),gene_short_name=rownames(st.m.epi@assays$Spatial@counts),
                                          row.names = rownames(st.m.epi@assays$Spatial@counts)))
cds.st.epi[['From']] = factor(cds.st.epi[['From']],levels = c('mLN','mLN+ PT','mLN- PT'))
plot_genes_violin(cds.st.epi[c('IGHA1','IGHG1','IGHG4','IGLC2'),],group_cells_by = 'From',
                  ncol = 2)+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(2,3),c(1,3)))+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())+
  scale_fill_manual(values = cols.From)

epiDistance = function(st.f,sample,input.cells,epi.cells){
  
  i.coord = st.f@images[[sample]]@coordinates
  
  distance.mat = matrix(,nrow = length(input.cells),ncol = length(epi.cells),
                        dimnames = list(input.cells,
                                        epi.cells))
  
  for(i in 1:nrow(distance.mat)){
    row.i.coord = i.coord[input.cells[i],]
    for(j in 1:ncol(distance.mat)){
      
      col.j.coord = i.coord[epi.cells[j],]
      
      distance.mat[i,j] = sqrt((row.i.coord[1,'row']-col.j.coord[1,'row'])**2+(row.i.coord[1,'col']-col.j.coord[1,'col'])**2)
      
      
    }
  }
  
  return(distance.mat)
  
}

st.f = st.m[,st.m$Source == 'PLM']

samples = unique(st.f$orig.ident)

unique(st.f$orig.ident)

epi.spots = c()
orig.flags = c()
epi.distance.min = c()

for(sample in samples){
  
  st.f.i.epi = colnames(st.f[,st.f$orig.ident == sample & st.f$cluster.name == 'Epithelial'])
  st.f.i.epi.not = colnames(st.f[,st.f$orig.ident == sample & st.f$cluster.name != 'Epithelial'])
  
  st.f.i.distance = epiDistance(st.f,sample,st.f.i.epi,st.f.i.epi.not)
  epi.spots = append(epi.spots,st.f.i.epi)
  orig.flags = append(orig.flags,rep(sample,length(st.f.i.epi)))
  epi.distance.min = append(epi.distance.min,
                            apply(st.f.i.distance,1,min))
  print(sample)
}

epi.distance.res = data.frame(row.names = epi.spots,
                              distance.min = epi.distance.min,
                              sample = orig.flags)
epi.distance.res$patient = substr(epi.distance.res$sample,1,2)

epi.distance.res$distance.tag = ''


ggplot(epi.distance.res,aes(x=sample,y=distance.min,fill=sample))+geom_boxplot()


dplyr::group_by(epi.distance.res,patient) %>% dplyr::summarise(mean(distance.min))
# 1 P3                      1.85
# 2 P5                      2.98
# 3 P6                      2.82
# 4 P7                      1.54
epi.distance.res[epi.distance.res$patient %in% c('P5','P6'),'distance.tag'] = ifelse(epi.distance.res[epi.distance.res$patient %in% c('P5','P6'),'distance.min'] < 3,'Edge','Inside')
epi.distance.res[epi.distance.res$patient %in% c('P3','P7'),'distance.tag'] = ifelse(epi.distance.res[epi.distance.res$patient %in% c('P3','P7'),'distance.min'] < 2,'Edge','Inside')

saveRDS(epi.distance.res,file = './variables_v2/epi.distance.res.Rds')



st.f.epi = st.f[,st.f$cluster.name == 'Epithelial']
st.f.epi$distance.tag = epi.distance.res[colnames(st.f.epi),'distance.tag']

st.f.epi = NormalizeData(st.f.epi)
SpatialDimPlot(st.f.epi,cols = c( 'Edge'='red','Inside'='blue'),
               group.by =c('distance.tag'),
               images = c('P3L','P5L','P6L','P7L','P3T','P5T','P6T','P7T'), 
               pt.size.factor = 1.5,
               stroke = 0,
               image.alpha = 0.5,
               ncol=4)

SpatialDimPlot(st.f.epi,
               group.by =c('distance.tag'),
               images = c('P6T','P6L'), 
               pt.size.factor = 1.5,
               stroke = 0,
               image.alpha = 0.5,
               ncol=2)

Idents(st.f.epi)=st.f.epi$distance.tag
markers.distance.epi.2 = FindAllMarkers(st.f.epi,only.pos = T,min.pct = 0.01,
                                        logfc.threshold = 0.1)
#IGLC2,IGLC3,IGHM,IGHG4,IGHG1,IGHG3
cds.st.epi = new_cell_data_set(as(as.matrix(st.f.epi@assays$Spatial@counts), "sparseMatrix"),
                               cell_metadata = st.f.epi@meta.data,
                               data.frame(id = rownames(st.f.epi@assays$Spatial@counts),gene_short_name=rownames(st.f.epi@assays$Spatial@counts),
                                          row.names = rownames(st.f.epi@assays$Spatial@counts)))

plot_genes_violin(cds.st.epi[c('IGHA1','IGHG1','IGHG4','IGLC2'),],group_cells_by = 'distance.tag',
                  ncol = 4)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())
plot_genes_violin(cds.st.epi[c('COL1A1','COL3A1','COL1A2','FN1'),],group_cells_by = 'distance.tag',
                  ncol = 4)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())

plot_genes_violin(cds.st.epi[c('B2M','HLA-B','HLA-DRA','HLA-A'),],group_cells_by = 'distance.tag',
                  ncol = 4)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())



saveRDS(st.f.epi,file = './variables_v2/st.f.epi.Rds')

st.f.epi = readRDS(file = './variables_v2/st.f.epi.Rds')

write.csv(markers.distance.epi.2,file='./tables_v2/ST markers for Epi Edge vs Inside.csv')
markers.distance.epi.2 = read.csv(file='./tables_v2/ST markers for Epi Edge vs Inside.csv',row.names = 1)

SpatialFeaturePlot(st.f.epi,features = c('IGHG4'),pt.size.factor = 1.5,stroke = 0,images = 'P5T')

st.f.epi.Plasma = st.f[,st.f$cluster.name %in% c('Epithelial','Plasma','Plasma+Myeloid')]

p2=SpatialDimPlot(st.f.epi.Plasma,group.by = 'cluster.name',
               pt.size.factor = 1.5,images = c('P6T','P6L'),cols = cols.cl)
p1=SpatialFeaturePlot(st.f.epi,features = c('IGHG4'),
                   pt.size.factor = 1.5,stroke = 0,images = c('P6T','P6L'))
SpatialFeaturePlot(st.f.epi,features = c('COL1A1'),
                   pt.size.factor = 1.5,stroke = 0,images = c('P6T','P6L'))

#####Edge and Inside pathway####
markers.distance.epi.edge = filter(markers.distance.epi.2,cluster == 'Edge' &
                                     p_val_adj<0.01)
markers.distance.epi.inside = filter(markers.distance.epi.2,cluster == 'Inside' &
                                       p_val_adj<0.01)

edge.path = enricher(markers.distance.epi.edge$gene,TERM2GENE = path2gene.all,
                     pvalueCutoff = 1.1)@result

inside.path = enricher(markers.distance.epi.inside$gene,TERM2GENE = path2gene.all,
                       pvalueCutoff = 1.1)@result


edge.path$From = 'Edge'
inside.path$From = 'Inside'

path.res = rbind(edge.path,inside.path)
res.filter = filter(path.res,pvalue<0.01 ) %>% group_by(From) %>% top_n(-5,pvalue) %>% arrange(From,pvalue) 

res.filter$Description = factor(res.filter$Description,levels = res.filter$Description)

ggplot(res.filter,aes(x=-log10(p.adjust),y=Description,fill=From))+geom_bar(stat = 'identity',width = 0.65)+
  facet_grid(~From,scales = 'free_y',space = 'free')+geom_vline(xintercept = 2,lty=3)+
  geom_vline(xintercept = 0,lty=1)+
  theme_cowplot(font_size = 8)+
  theme(axis.text = element_text(size=8),
        axis.line.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.title = element_text(size=8),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size = 8))
saveRDS(path.res,file = './variables_v2//Epi inside edge pathways_dis4.Rds')
path.res = readRDS(file = './variables_v2//Epi inside edge pathways_dis4.Rds')

#####Spatial CytoTRACE####
expr.st = as.matrix(st.f.epi@assays$Spatial@counts)

# epi.data.LUSC = preprocess_cds(epi.data.LUSC)
# epi.data.LUSC = reduce_dimension(epi.data.LUSC,cores = 20)
# 
selected_genes <- rownames(expr.st)[ Matrix::rowSums(expr.st > 0)> 50]
cyto.st <- CytoTRACE(expr.st[selected_genes,],batch = as.character(substr(st.f.epi$orig.ident,1,2)),ncores = 12)

st.f.epi$CytoTRACE = cyto.st$CytoTRACE

SpatialFeaturePlot(st.f.epi,features = c('CytoTRACE'),
                   pt.size.factor = 1.5,stroke = 0,images = c('P6T','P6L'))
saveRDS(cyto.st,file = './variables_v2/spatial cytoTRACE.Rds')

VlnPlot(st.f.epi,features = 'CytoTRACE',group.by = 'distance.tag',pt.size = 0)

meta.data = st.f.epi@meta.data
ggplot(meta.data,aes(x=distance.tag,y=CytoTRACE,fill=distance.tag))+geom_boxplot()+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)
#####NicheNet spatial interaction#####
library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
sc.t = st.m[,st.m$Source == 'PLM']
Idents(sc.t)=sc.t$cluster.name

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(file = './data/ligand_target_matrix_nsga2r_final.rds')
weighted_networks <- readRDS(file = "./data/weighted_networks_nsga2r_final.rds")

receiver = "Epithelial"
expressed_genes_receiver <- get_expressed_genes(receiver, sc.t, pct = 0.1)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()


sender_celltypes <- c('Plasma','Plasma+Myeloid','Fibroblast','Myeloid')

list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, sc.t, 0.25)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
## [1] 8702
length(potential_ligands)
## [1] 597
length(potential_ligands_focused)
## [1]365

geneset_oi <- unique(rownames(markers.distance.epi.2[markers.distance.epi.2$avg_log2FC>0.5 & 
                                                       markers.distance.epi.2$p_val_adj < 0.05 &
                                                       markers.distance.epi.2$cluster=='Edge',]))


background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
## [1] 6299
length(geneset_oi)
## [1] 1824

# geneset_oi = geneset_oi[grepl('IG',geneset_oi)]
# geneset_oi = geneset_oi[1:13]
# geneset_oi = geneset_oi[-4]
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = order(aupr_corrected,decreasing = T))
ligand_activities


best_upstream_ligands <- ligand_activities %>% top_n(50, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.1) 


order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
vis_ligand_target.row.mean = apply(vis_ligand_target, 1, mean)
vis_ligand_target.row.mean = sort(vis_ligand_target.row.mean,decreasing = T)
vis_ligand_target.high = vis_ligand_target[names(vis_ligand_target.row.mean)[20:1],]
vis_ligand_target.high = vis_ligand_target.high[,apply(vis_ligand_target.high, 2, sum) >0]
make_heatmap_ggplot(vis_ligand_target.high, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  theme(axis.text =  element_text(size = 8))+
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")
saveRDS(vis_ligand_target,file = './variables_v2/nicheNet_epiEdge.Rds')


VlnPlot(sc.t,features = c('TGFB1','TGFB2','BMP4','IL13','INHBA','POSTN','CCN2','DKK1'),group.by = 'cluster.name',pt.size = 0,adjust = 5)
cds.st = new_cell_data_set(as(as.matrix(sc.t@assays$Spatial@counts), "sparseMatrix"),
                               cell_metadata = sc.t@meta.data,
                               data.frame(id = rownames(sc.t@assays$Spatial@counts),gene_short_name=rownames(sc.t@assays$Spatial@counts),
                                          row.names = rownames(sc.t@assays$Spatial@counts)))

plot_genes_violin(cds.st[ c('IL13','SLPI','LPA'),],group_cells_by = 'cluster.name',
                  ncol = 4)+
  ggpubr::stat_compare_means()+
  theme_cowplot(font_size = 8)+
  theme(axis.text.x = element_text(angle = 90),
        strip.background = element_blank())
