library(pheatmap)

cds.all = readRDS(file = './variables_v2/inter_cds_refine.Rds')
cds = cds.all[,cds.all[['Source']] == 'PLM']#14787 334956
meta.data.cds = read.csv(file = './tables_v2/meta.data.cds with malignant or not Epi.csv',row.names = 1)##########################re-analysis need!!!!!!!!!!!!!!!!!!!!!!!
cds[['subCellType_v2']] = meta.data.cds$subCellType_v2##########################re-analysis need!!!!!!!!!!!!!!!!!!!!!!!
epi.data.mal = readRDS(file = './epi.data.mal.Rds')
epi.data.mal[['IgG_Group']] = ifelse(epi.data.mal[['cluster.epi.mal']] %in%
                                       c('MC21','MC 2','MC 3','MC 9','MC13',"MC 15",'MC 27','MC16'),'IgG+','IgG-')

meta.data.cds = as.data.frame(colData(cds))
meta.data.cds$subCellType_v3 = meta.data.cds$subCellType_v2
meta.data.cds[colnames(epi.data.mal),'subCellType_v3'] = epi.data.mal[['IgG_Group']]

meta.data.cds$subCellType_v4 = meta.data.cds$subCellType
meta.data.cds[colnames(epi.data.mal),'subCellType_v4'] = paste(meta.data.cds[colnames(epi.data.mal),'subCellType'],epi.data.mal[['IgG_Group']])

write.csv(meta.data.cds[,c('orig.ident','nCount_RNA','nFeature_RNA','PatientID','scS','Disease',
                           'percent_mito','Source','N0','From','assigned_cell_type','subCellType',
                           "majorCellType2","subCellType_v2", "subCellType_v3",'majorCellType', "subCellType_v4")],file = './tables_v2/meta.data.cds main info.csv')
PLM.infor = read.csv(file = './documents/Patient Info.csv')
PLM.infor$Stage2 = gsub('A|B','',PLM.infor$Stage)
rownames(PLM.infor)=PLM.infor$Patient_ID
com2orig.2 = table(meta.data.cds$orig.ident,meta.data.cds$subCellType_v2)

f.data<-as.matrix(com2orig.2)
rownames(f.data)[row.names(f.data) == 'P0.3A']='P0.L3'
f.data = f.data/apply(f.data, 1, sum)
f.data = scale(f.data,scale=T,center =T)
t.rows = unlist(lapply(rownames(f.data[grepl('T',rownames(f.data)),]),
                       function(a){strsplit(a,'\\.')[[1]][1]}))
pheatmap(f.data[grepl('T',rownames(f.data)),],
         clustering_distance_rows='correlation',
                   annotation_row = data.frame(row.names = rownames(f.data[grepl('T',rownames(f.data)),]),
                                               Stage = PLM.infor[t.rows,'Stage2'],
                                               N_Stage=PLM.infor[t.rows,'N_stage'],
                                               M_Stage = PLM.infor[t.rows,'M_stage'],
                                               Smoking = PLM.infor[t.rows,'Smoking_history'],
                                               Disease = PLM.infor[t.rows,'Disease']))
l.rows = unlist(lapply(rownames(f.data[grepl('L',rownames(f.data)),]),
                       function(a){strsplit(a,'\\.')[[1]][1]}))

pheatmap(f.data[grepl('L',rownames(f.data)),],
                   clustering_distance_rows='correlation',
                   annotation_row = data.frame(row.names = rownames(f.data[grepl('L',rownames(f.data)),]),
                                               Stage = PLM.infor[l.rows,'Stage2'],
                                               N_Stage=PLM.infor[l.rows,'N_stage'],
                                               M_Stage = PLM.infor[l.rows,'M_stage'],
                                               Smoking = PLM.infor[l.rows,'Smoking_history'],
                                               Disease = PLM.infor[l.rows,'Disease']))
rows = unlist(lapply(rownames(f.data),
                       function(a){strsplit(a,'\\.')[[1]][1]}))

pheatmap(f.data,
         clustering_distance_rows='correlation',
         annotation_row = data.frame(row.names = rownames(f.data),
                                     Stage = PLM.infor[rows,'Stage2'],
                                     N_Stage=PLM.infor[rows,'N_stage'],
                                     M_Stage = PLM.infor[rows,'M_stage'],
                                     Smoking = PLM.infor[rows,'Smoking_history'],
                                     Disease = PLM.infor[rows,'Disease'],
                                     From = unlist(lapply(rownames(f.data),
                                                          function(a){substr(strsplit(a,'\\.')[[1]][2],1,1)}))))


corrs.cell = cor(f.data,method = 'spearman')

pheatmap(corrs.cell,color = colorRampPalette(rev(hcl.colors(100,'Spectral')))(100),
         fontsize = 8)

library(reshape2)
com2orig.2 = table(meta.data.cds$orig.ident,meta.data.cds$subCellType_v2)

f.data<-as.matrix(com2orig.2)
rownames(f.data)[row.names(f.data) == 'P0.3A']='P0.L3'
f.data = f.data/apply(f.data, 1, sum)

f.data.change = melt(f.data,measure.vars = colnames(f.data))
f.data.change$PatientID = unlist(lapply(as.character(f.data.change$Var1),
                                        function(a){strsplit(a,'\\.')[[1]][1]}))
f.data.change$From = unlist(lapply(as.character(f.data.change$Var1),
                                          function(a){substr(strsplit(a,'\\.')[[1]][2],1,1)}))
f.data.change$Disease = PLM.infor[f.data.change$PatientID,'Disease']
f.data.change$N_Stage = substr(PLM.infor[f.data.change$PatientID,'N_stage'],1,2)
f.data.change$M_Stage = PLM.infor[f.data.change$PatientID,'M_stage']
f.data.change$Smoking = PLM.infor[f.data.change$PatientID,'Smoking_history']
f.data.change$Stage = PLM.infor[f.data.change$PatientID,'Stage2']

write.csv(f.data.change,file = './tables_v2/plm cds f.data.change.csv')
ggplot(f.data.change,aes(x=N_Stage,fill=N_Stage,y=value))+geom_boxplot()+
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),size=2)+
  facet_wrap(~Var2,nrow=4)+
  theme_cowplot(font_size = 8)


ggplot(f.data.change,aes(x=Stage,fill=Stage,y=value))+geom_boxplot()+
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),comparisons = list(c(1,2),c(2,3)),size=2)+
  facet_wrap(~Var2,nrow=4,scales = 'free')+
  theme_cowplot(font_size = 8)


ggplot(f.data.change[f.data.change$From == 'T',],aes(x=Smoking,fill=Smoking,y=value))+geom_boxplot(outlier.size = 0.25)+
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),size=2,method = 't.test')+
  facet_wrap(~Var2,nrow=4,scales = 'free')+
  theme_cowplot(font_size = 8)+
  theme(strip.background = element_blank())


ggplot(f.data.change[f.data.change$From == 'L',],aes(x=Smoking,fill=Smoking,y=value))+geom_boxplot(outlier.size = 0.25)+
  ggpubr::stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),size=2,method = 't.test')+
  facet_wrap(~Var2,nrow=4,scales = 'free')+
  theme_cowplot(font_size = 8)+
  theme(strip.background = element_blank())
