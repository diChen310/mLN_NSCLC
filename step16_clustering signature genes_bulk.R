.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')


library(MuSiC)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(cowplot)
library(ggplot2)
cds.all = readRDS(file = './variables_v2/inter_cds_refine.Rds')

marker_test_res = top_markers(cds.all, group_cells_by="subCellType", genes_to_test_per_group=50,
                              reference_cells=1000, cores=24)
write.csv(marker_test_res,file = './tables_v2/all_refine_markers.csv')
marker_test_res = read.csv(file = './tables_v2/all_refine_markers.csv')

#################Based on the genes##############################
#int.genes = intersect(rownames(dataBulk),marker_test_res$gene_id)
tumor = 'LUAD'
clin.luad = read.delim(file = './data/LUAD_clin.txt')

notDead <- is.na(clin.luad$days_to_death)
if (length(notDead) > 0) {
  clin.luad[notDead,]$days_to_death <-
    clin.luad[notDead,]$days_to_last_follow_up
}



mat.1<-read.delim(file=paste0("./data/",tumor,'ExpressionNormalized.txt'),sep='\t',check.names = FALSE)
marker_test_res.sig = marker_test_res[marker_test_res$marker_test_p_value<0.001,]
both.genes = intersect(union(unique(marker_test_res$gene_id),both.cyto.genes),rownames(mat))

dataBulk = mat.1[both.genes,]

samplesTP<-colnames(dataBulk)[substr(colnames(dataBulk),14,15)=='01']
samplesTP.2<-colnames(mat.2)[substr(colnames(mat.2),14,15)=='01']

# 
dataBulk.tp = dataBulk[,samplesTP]
colnames(dataBulk.tp) = substr(colnames(dataBulk.tp),1,12)

rownames(clin.luad)=clin.luad$bcr_patient_barcode
clin.cls = intersect(clin.luad$bcr_patient_barcode,colnames(dataBulk.tp))
rownames(clin.luad)=clin.luad$submitter_id
clin.av = clin.luad[clin.cls,]
clin.av$stage = gsub('a|b','',clin.av$tumor_stage)
#clin.av = cbind(clin.av,t(dataBulk.tp[,clin.cls]) )
clin.av$status <- grepl("dead", clin.av$vital_status, ignore.case = TRUE)

clin.av$days_to_death = as.double(clin.av$days_to_death)



exp.sel = log2(dataBulk.tp+1)
cluster.res = ConsensusClusterPlus(d=scale(exp.sel),maxK = 4,clusterAlg = 'pam',seed = 33355 )
saveRDS(cluster.res,file = './variables_v2/LUAD_genes_clusters.Rds')
cluster.res.2 = cluster.res[[4]]$consensusClass

clin.av$SCC = paste('SCC',cluster.res.2[rownames(clin.av)])
survdiff(Surv(days_to_death,status) ~ SCC,clin.av)

ggsurvplot(survfit(Surv(days_to_death,status) ~ SCC,clin.av),data = clin.av,pval = T,
           conf.int = F)
ggsurvplot(survfit(Surv(days_to_death,status) ~ SCC=='SCC 3',clin.av),data = clin.av,pval = T,
           conf.int = F)
survdiff(Surv(days_to_death,status) ~ SCC=='SCC 3',clin.av)
clin.av$Group = ifelse(clin.av$SCC == 'SCC 3','Poor','Others')

ggsurvplot(survfit(Surv(days_to_death,status) ~ Group,clin.av),
           pval = 5e-06, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_cowplot(font_size = 8), # Change ggplot2 theme
           palette = c("lightblue","#e58b8e")
)
scc4.samples = rownames(clin.av[clin.av$SCC == 'SCC 3',])
other.samples = rownames(clin.av[clin.av$SCC != 'SCC 3',])

scc4.p = unlist(lapply(1:nrow(dataBulk.tp), function(i){
  wilcox.test(as.double(exp.sel[i,scc4.samples]),as.double(exp.sel[i,other.samples]))$p.value
}))
scc4.FC = unlist(lapply(1:nrow(dataBulk.tp), function(i){
 mean(as.double(log2(dataBulk.tp[i,scc4.samples]+1)))-
    mean(log2(as.double(dataBulk.tp[i,other.samples])+1))
}))
res.scc4 = data.frame(gene=rownames(dataBulk.tp),pvalue = scc4.p,log2FC=scc4.FC)
res.scc4=res.scc4[order(res.scc4$pvalue),]
#res.scc4.sig = res.scc4[res.scc4$pvalue<0.05,]
#res.scc4.sig = res.scc4.sig[order(res.scc4.sig$log2FC),]
boxplot(as.double(exp.sel['FTH1',scc4.samples]),as.double(exp.sel['FTH1',other.samples]))

wilcox.test(as.double(exp.sel['FTH1',scc4.samples]),as.double(exp.sel['FTH1',other.samples]))

write.csv(res.scc4,file = './tables_v2/dataBulk_LUAD_clustering_diff.csv')


load('/share/dichen/lungNodeM/path.info.RData')

exp.sel = log2(mat.1[,samplesTP]+1)
colnames(exp.sel)=substr(colnames(exp.sel),1,12)
scc4.p = unlist(lapply(1:nrow(exp.sel), function(i){
  wilcox.test(as.double(exp.sel[i,scc4.samples]),as.double(exp.sel[i,other.samples]))$p.value
}))
scc4.FC = unlist(lapply(1:nrow(exp.sel), function(i){
  mean(as.double(exp.sel[i,scc4.samples]))-
    mean(as.double(exp.sel[i,other.samples]))
}))
res.scc4.2 = data.frame(gene=rownames(exp.sel),pvalue = scc4.p,log2FC=scc4.FC)
names(scc4.FC)=rownames(exp.sel)
luad.mut = readRDS(file = '/share/ubuntu/dichen/KRAS/luad.mut.Rds')
mutsig.luad <- read.delim(file='/share/ubuntu/dichen/KRAS/LUAD-TCGA.sig_genes.txt')
luad.top20 <- as.character(mutsig.luad$gene)[1:20]
luad.sig.genes <- as.character(mutsig.luad[mutsig.luad$q < 0.1,'gene'])

luad.mut.scc3= subsetMaf(luad.mut,tsb = scc4.samples)
luad.mut.other= subsetMaf(luad.mut,tsb = other.samples)

oncoplot(luad.mut.scc3,genes = luad.sig.genes)
oncoplot(luad.mut.other,genes = luad.sig.genes)
diff.mut = mafCompare(m1=luad.mut.scc3,m2=luad.mut.other,m1Name = 'BC3',m2Name = 'Others')
diff.mut.res = diff.mut$results
forestPlot(mafCompareRes = diff.mut, pVal = 0.005, color = c('royalblue', 'maroon'), geneFontSize = 0.8)


saveRDS(diff.mut,file = './variables_v2/TCGA-LUAD BC3 diff mutations.Rds')
library(clusterProfiler)
gsea.scc4 = GSEA(scc4.FC[order(scc4.FC,decreasing = T)],pvalueCutoff = 1.2,TERM2GENE = path.info)
gsea.res = gsea.scc4@result
library(pheatmap)

cls = names(cluster.res.2)[cluster.res[[4]]$consensusTree$order]
# exp.sel = mat[,samplesTP]
# colnames(exp.sel) = substr(colnames(exp.sel),1,12)
# exp.sel = log2(exp.sel+1)
mat.input = t(scale(t(exp.sel[c('MUC21','NAPSA','C16orf89','SFTPC','PGC',
                                    'SFTPA1','HLA-DQB2','SLC34A2','CD74','HLA-DPB1',
                                    'HLA-DRB1','HLA-DRA','CD37','CXCR3','CD1E',
                                    'FTH1','PBK','AKR1C2','AKR1C1','GPX2',
                                    'MAGEA3','KRT6A','ALDH3A1','AKR1C3','KRT17',
                                    'TOP2A','MKI67','AURKB','BIRC5','TPX2'),cls])))
mat.input[mat.input<=-3]=-3
pheatmap(mat.input,cluster_cols = F,
         color = hcl.colors(256),
         annotation_col = data.frame(SSC = clin.av[cls,'SCC'],row.names = cls,
                                     Gender = clin.av[cls,'gender'],
                                     Race = clin.av[cls,'race'] ,
                                     stage = clin.av[cls,'stage']),
         show_colnames = F,
         fontsize = 8,
         annotation_colors = list(SSC = c('SCC 1'='red','SCC 2'='orange','SCC 3'='blue','SCC 4'='cyan'))
)


save.image('deconvBulkLuad.RData')




####GEO validation
library(GEOquery)
gseid<-"GSE30219"
Data <- getGEO(gseid,destdir="./",parseCharacteristics = T)
#Data = getGEO(filename='E:/project/metaboliteProteinInteraction/GSE63898_series_matrix.txt.gz',destdir="./")
myMatrix <- Data[[1]]@assayData$exprs
(nrow(myMatrix))
(ncol(myMatrix))

#Phenotype infomation
myPDfile <- phenoData(Data[[1]])
myPDfile = myPDfile@data
(nrow(myPDfile))
(ncol(myPDfile))
#myPDfile$`disease state:ch1`= ifelse(grepl('Non tumor NT',as.character(myPDfile$title)),'N','T')
##Change gene id into gene symbos
genes.raw.names <- rownames(myMatrix)
#gpl <- getGEO(filename='GPL3921.soft')

gene.anno = Data[[1]]@featureData@data
#gene.anno$Gene_symbol = unlist(lapply(gene.anno$gene_assignment,function(a){defineGeneName(a,gene.info)}))
#gene.anno = Table(gpl)

symbol <- gene.anno[c('ID','Symbol','GB_ACC')]
colnames(symbol)<-c('ID','Gene Symbol','Definition')
head(symbol)
symbol <- symbol[!is.na(symbol$ID),]
#symbol <- symbol[!is.na(symbol$`Gene Symbol`),]

rownames(symbol) <- symbol$ID
rows.symbols = symbol[as.character(rownames(myMatrix)),'Gene Symbol']
dup.rows = unique(rows.symbols[duplicated(rows.symbols)])
uni.rows = rows.symbols[!duplicated(rows.symbols)]
uni.rows = uni.rows[uni.rows %in% dup.rows == F]
uni.ids = symbol[symbol$`Gene Symbol` %in% uni.rows,'ID']
uni.mat <- myMatrix[uni.ids,]
rownames(uni.mat) <- symbol[uni.ids,'Gene Symbol']

##Merge duplicated gene ids by the median value
dup.mat <- matrix(0,nrow=length(dup.rows),ncol=ncol(myMatrix),dimnames = list(dup.rows,colnames(myMatrix)))
for( i in 1:length(dup.rows)){
  dup = dup.rows[i]
  ids.dup = symbol[symbol$`Gene Symbol` == dup,'ID']
  ids.dup = intersect(ids.dup,rownames(myMatrix))
  dup.mat[i,] = apply(myMatrix[ids.dup,],2,median)
}
matrix.retain = rbind(uni.mat,dup.mat)

save(matrix.retain,file = './variables/GSE10143_retainMatrix.RData')

save(myPDfile,file = './variables/GSE10143_pdFile.RData')

#myPDfile = read.csv(file = './data/GEO/GSE30219_clin.csv',row.names = 1)
#####GSE30219#####
matrix.retain <- read.csv(file = './data/GEO/GSE30219_matrix_retain.csv',row.names = 1)
matrix.retain <- log2(matrix.retain+1)
mat.new = scale(t(matrix.retain))
control <- trainControl(method='CV',number=10)
train.x = t(scale(exp.sel))[,intersect(rownames(exp.sel),colnames(mat.new))]
ene.train.2 <- train(x= train.x,
                   y= ifelse(cluster.res[[4]]$consensusClass == '3','Poor','Other'), 
                   method = 'glmnet',
                   trControl = control)
ene.train.2
#Accuracy was used to select the optimal model using the largest value.
#The final values used for the model were alpha = 0.1 and lambda = 0.03411048.
enet <- glmnet(x=train.x,y=ifelse(cluster.res[[4]]$consensusClass == '3','Poor','Other'),family = 'binomial',
               alpha = 0.1,lambda =0.1078668)



gse.pred = predict(enet,newx = mat.new[,rownames(enet$beta)],type = 'class')
clin_gse <- read.csv(file = './data/GEO/GSE30219_clin.csv',row.names = 1)
clin_gse = clin_gse[,c('id','status','days_to_death')]
clin_gse = clin_gse[clin_gse$id %in% rownames(mat.new),]
rownames(clin_gse) = clin_gse$id
clin_gse$Group = as.character(gse.pred[rownames(clin_gse),'s0'])
clin_gse$status = grepl("dead", clin_gse$status, ignore.case = TRUE)
km.by.cluster = survfit(Surv(days_to_death, status) ~ Group, data = clin_gse)
survdiff(Surv(days_to_death, status) ~ Group, data = clin_gse,rho = 1)
ggsurvplot(km.by.cluster,
           pval =0.009, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_cowplot(font_size = 8), # Change ggplot2 theme
           palette = c("lightblue","#e58b8e")
)
saveRDS(clin_gse,file = './variables_v2/GSE30219_clin.Rds')
rownames(myPDfile)=myPDfile$geo_accession

myPDfile$Group = clin_gse[myPDfile$geo_accession,'Group']

ggplot(myPDfile[!is.na(myPDfile$Group),],aes(x=Group,fill=`pn stage:ch1`))+geom_bar(position = 'fill')+
  theme_cowplot()

chisq.test(table(myPDfile[!is.na(myPDfile$Group),]$Group,myPDfile[!is.na(myPDfile$Group),]$`pn stage:ch1` == 'N0')) #p-value = 0.0002672

clin_gse = clin_gse[order(clin_gse$days_to_death),]
clin_gse = clin_gse[order(clin_gse$status),]

clin_gse = clin_gse[order(clin_gse$Group),]
cls = rownames(clin_gse)
mat.input =mat.new[cls,intersect(c('CXCR3','CD37','HLA-DPB1','CD74','HLA-DRA','HLA-DQB2','CD1E',
                                   'PGC','SFTPC','C16orf89','NAPSA','SLC34A2','PBK','AURKB',
                                   'BIRC5','TPX2','TOP2A','MKI67','KRT6A','AKR1C3','ALDH3A1','AKR1C1',
                                   'GPX2','FTH1'),colnames(mat.new))]
mat.input[mat.input<=-3]=-3
mat.input[mat.input>=3.5]=3.5
pheatmap(t(mat.input),cluster_cols = F,clustering_distance_rows = 'correlation',
         color = hcl.colors(256),
         annotation_col = data.frame(Group = clin_gse[cls,'Group'],row.names = cls),
         show_colnames = F,
         fontsize = 8,cluster_rows = F)


####library(GEOquery)
gseid<-"GSE68465"
Data <- getGEO(gseid,destdir="./",parseCharacteristics = T)
#Data = getGEO(filename='E:/project/metaboliteProteinInteraction/GSE63898_series_matrix.txt.gz',destdir="./")
myMatrix <- Data[[1]]@assayData$exprs
(nrow(myMatrix))
(ncol(myMatrix))

#Phenotype infomation
myPDfile <- phenoData(Data[[1]])
myPDfile = myPDfile@data
(nrow(myPDfile))
(ncol(myPDfile))
#myPDfile$`disease state:ch1`= ifelse(grepl('Non tumor NT',as.character(myPDfile$title)),'N','T')
##Change gene id into gene symbos
genes.raw.names <- rownames(myMatrix)
#gpl <- getGEO(filename='GPL3921.soft')

gene.anno = Data[[1]]@featureData@data
#gene.anno$Gene_symbol = unlist(lapply(gene.anno$gene_assignment,function(a){defineGeneName(a,gene.info)}))
#gene.anno = Table(gpl)

symbol <- gene.anno[c('ID','Gene Symbol','GB_ACC')]
colnames(symbol)<-c('ID','Gene Symbol','Definition')
head(symbol)
symbol <- symbol[!is.na(symbol$ID),]
#symbol <- symbol[!is.na(symbol$`Gene Symbol`),]

rownames(symbol) <- symbol$ID
rows.symbols = symbol[as.character(rownames(myMatrix)),'Gene Symbol']
dup.rows = unique(rows.symbols[duplicated(rows.symbols)])
uni.rows = rows.symbols[!duplicated(rows.symbols)]
uni.rows = uni.rows[uni.rows %in% dup.rows == F]
uni.ids = symbol[symbol$`Gene Symbol` %in% uni.rows,'ID']
uni.mat <- myMatrix[uni.ids,]
rownames(uni.mat) <- symbol[uni.ids,'Gene Symbol']

##Merge duplicated gene ids by the median value
dup.mat <- matrix(0,nrow=length(dup.rows),ncol=ncol(myMatrix),dimnames = list(dup.rows,colnames(myMatrix)))
for( i in 1:length(dup.rows)){
  dup = dup.rows[i]
  ids.dup = symbol[symbol$`Gene Symbol` == dup,'ID']
  ids.dup = intersect(ids.dup,rownames(myMatrix))
  dup.mat[i,] = apply(myMatrix[ids.dup,],2,median)
}
matrix.retain = rbind(uni.mat,dup.mat)

write.csv(matrix.retain,file = './data/GEO/GSE68465_retainMatrix.csv')
myPDfile$status = myPDfile$characteristics_ch1.4
myPDfile$days_to_death = gsub('months_to_last_contact_or_death: ','',myPDfile$characteristics_ch1.11)
myPDfile$days_to_death = as.numeric(myPDfile$days_to_death)

write.csv(myPDfile,file = './data/GEO/GSE68465_clin.csv')
myPDfile = read.csv(file = './data/GEO/GSE68465_clin.csv',row.names=1)
#####GSE68465#####
matrix.retain <- read.csv(file = './data/GEO/GSE68465_retainMatrix.csv',row.names = 1)
matrix.retain <- log2(matrix.retain+1)
mat.new = scale(t(matrix.retain))
control <- trainControl(method='CV',number=10)
train.x = t(scale(exp.sel))[,intersect(rownames(exp.sel),colnames(mat.new))]
ene.train.2 <- train(x= train.x,
                     y= ifelse(cluster.res[[4]]$consensusClass == '3','Poor','Other'), 
                     method = 'glmnet',
                     trControl = control)
ene.train.2
# Accuracy was used to select the optimal model using the largest value.
# The final values used for the model were alpha = 0.1 and lambda = 0.1099218.

enet <- glmnet(x=train.x,y=ifelse(cluster.res[[4]]$consensusClass == '3','Poor','Other'),family = 'binomial',
               alpha = 0.1,lambda =0.1099218)



gse.pred = predict(enet,newx = mat.new[,rownames(enet$beta)],type = 'class')
clin_gse <- myPDfile
myPDfile$id = rownames(myPDfile)
clin_gse = myPDfile[,c('id','status','days_to_death')]
clin_gse = clin_gse[clin_gse$id %in% rownames(mat.new),]
clin_gse = clin_gse[!is.na(clin_gse$days_to_death),]
rownames(clin_gse) = clin_gse$id
clin_gse$Group = as.character(gse.pred[rownames(clin_gse),'s0'])
clin_gse$status = grepl("Dead", clin_gse$status, ignore.case = TRUE)
km.by.cluster = survfit(Surv(days_to_death, status) ~ Group, data = clin_gse)
survdiff(Surv(days_to_death, status) ~ Group, data = clin_gse,rho = 1)
ggsurvplot(km.by.cluster,
           pval =T, conf.int = FALSE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_cowplot(font_size = 8), # Change ggplot2 theme
           palette = c("blue","red")
)
write.csv(clin_gse[,c('days_to_death','status','Group')],file = './tables_v3/fig s8c.csv')
rownames(myPDfile)=myPDfile$geo_accession
myPDfile =myPDfile[rownames(clin_gse),]

myPDfile$Group = clin_gse[myPDfile$geo_accession,'Group']
myPDfile$NStage = substr(myPDfile$characteristics_ch1.7,17,18)

ggplot(myPDfile[myPDfile$NStage != 'p',],aes(x=Group,fill=NStage))+geom_bar(position = 'fill')+
  theme_cowplot()

chisq.test(table(myPDfile[myPDfile$NStage != 'p',]$Group,myPDfile[myPDfile$NStage != 'p',]$NStage == 'N0')) #p-value = 0.0.0207

clin_gse = clin_gse[order(clin_gse$days_to_death),]
clin_gse = clin_gse[order(clin_gse$status),]

clin_gse = clin_gse[order(clin_gse$Group),]
cls = rownames(clin_gse)
mat.input =mat.new[cls,intersect(c('CXCR3','CD37','HLA-DPB1','CD74','HLA-DRA','HLA-DQB2','CD1E',
                                   'PGC','SFTPC','C16orf89','NAPSA','SLC34A2','PBK','AURKB',
                                   'BIRC5','TPX2','TOP2A','MKI67','KRT6A','AKR1C3','ALDH3A1','AKR1C1',
                                   'GPX2','FTH1'),colnames(mat.new))]
mat.input[mat.input<=-3]=-3
mat.input[mat.input>=3.5]=3.5
pheatmap(t(mat.input),cluster_cols = F,clustering_distance_rows = 'correlation',
         color = hcl.colors(256),
         annotation_col = data.frame(Group = clin_gse[cls,'Group'],row.names = cls),
         show_colnames = F,
         fontsize = 8)


