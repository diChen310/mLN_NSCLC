####inferCNV one by one
.libPaths(c( "/home/ubuntu/R/x86_64-pc-linux-gnu-library/4.2","/home/dichen/R/x86_64-pc-linux-gnu-library/4.2"
))
setwd('/share/dichen/lungNodeM')
library(infercnv)
library(dplyr)
library(pheatmap)
cds = readRDS(file = './variables_v2/inter_cds_refine.Rds')
cds.PLM = cds[,colData(cds)$Source == 'PLM']
meta.data.PLM = colData(cds.PLM)
####1. four samples with normal samples as reference####

for(pi in c('P3','P5','P6','P7')){ ##With normal samples
  
  cds.pi = cds.PLM[,cds.PLM[['PatientID']] == pi]
  cds.pi = cds.pi[,cds.pi[['majorCellType']]%in% c('Epi')]
  cells.T = colnames(cds.pi[,cds.pi[['From']]%in% c('mLN+ N')])
  
  # cds.pi = cds.pi[,cds.pi[['majorCellType']]%in% c('Epi','T')]
  # cells.T = colnames(cds.pi[,cds.pi[['majorCellType']]%in% c('T')])
  
  exp.rawdata = as.matrix(cds.pi@assays@data$counts)
  annotation = colData(cds.pi)
  annotation$CellTypeInput = paste0(annotation$PatientID,'_',annotation$From,'_','Epi')
  annotation$CellName = rownames(annotation)
  write.table(annotation[,c('CellName','CellTypeInput')],
              file = paste0('./documents/inferCNV_v2/each/',pi,'/anno.txt'),row.names = F,sep='\t',quote = F,col.names = F)
  ref.type = unique(as.character(annotation$CellTypeInput))[grepl(' N_Epi',unique(as.character(annotation$CellTypeInput)))]
  print(unique(annotation$CellTypeInput))
  
  infercnv_ti = CreateInfercnvObject(raw_counts_matrix = exp.rawdata,
                                     gene_order_file = './documents/inferCNV/genes location-v2.txt',
                                     annotations_file = paste0('./documents/inferCNV_v2/each/',pi,'/anno.txt'),
                                     ref_group_names = ref.type)
  
  gc()
  infercnv_obj = infercnv::run(infercnv_ti,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=paste0('./documents/inferCNV_v2/each/',pi),  
                               cluster_by_groups=T, denoise=TRUE,HMM=FALSE,
                               num_threads = 24,
                               output_format = 'pdf') 
  
  print(pi)
}

####2. the other samples with major immune cells as reference####

for(pi in unique(cds.PLM[['PatientID']])){
  if(pi %in% c('P3','P5','P6','P7') == F){
    cds.pi = cds.PLM[,cds.PLM[['PatientID']] == pi]
    cds.pi = cds.pi[,cds.pi[['majorCellType']]%in% c('Epi','Myeloid','T','B','Mast')]
    
    exp.rawdata = as.matrix(cds.pi@assays@data$counts)
    annotation = colData(cds.pi)
    annotation$CellTypeInput = paste0(annotation$PatientID,'_',annotation$From,'_',annotation$majorCellType)
    annotation$CellName = rownames(annotation)
    write.table(annotation[,c('CellName','CellTypeInput')],
                file = paste0('./documents/inferCNV_v2/each2/',pi,'/anno.txt'),row.names = F,sep='\t',quote = F,col.names = F)
    ref.type = unique(as.character(annotation$CellTypeInput))[!grepl('_Epi',unique(as.character(annotation$CellTypeInput)))]
    print(unique(as.character(annotation$CellTypeInput)))
    
    infercnv_ti = CreateInfercnvObject(raw_counts_matrix = exp.rawdata,
                                       gene_order_file = './documents/inferCNV/genes location-v2.txt',
                                       annotations_file = paste0('./documents/inferCNV_v2/each2/',pi,'/anno.txt'),
                                       ref_group_names = ref.type)
    
    gc()
    infercnv_obj = infercnv::run(infercnv_ti,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir=paste0('./documents/inferCNV_v2/each2/',pi),  
                                 cluster_by_groups=T, denoise=TRUE,HMM=FALSE,
                                 num_threads = 36,
                                 output_format = 'pdf') 
  }
}

####3.Merge the results from all patients####
annotations.row = data.frame()
# cl.cols = cols.default[1:10]
# names(cl.cols) = paste('CC',1:10)
cluster.scores = data.frame()
gene.sum = data.frame()
for(pi in paste0('P',0:18)){
  
  if(pi %in% c('P3','P5','P6','P7')){
    
    infercnv_obj = readRDS(file = paste0('./documents/inferCNV_v2/each/',pi,'/run.final.infercnv_obj'))
  }else{
    infercnv_obj = readRDS(file = paste0('./documents/inferCNV_v2/each2/',pi,'/run.final.infercnv_obj'))
    
  }
  expr = infercnv_obj@expr.data
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  
  names_test = c( names(test_loc)[grepl('PT',names(test_loc))],
                  names(test_loc)[grepl('mLN_',names(test_loc))])
  locs_test = Reduce(union,lapply(names_test, function(a){
    cl.pi.1 = infercnv_obj@tumor_subclusters[["hc"]][[a]]
    if(length(cl.pi.1$labels)>2){
      cl.pi.res.1 =  test_loc[[a]][infercnv_obj@tumor_subclusters[["hc"]][[a]][["order"]]]
      
    }else{
      cl.pi.res.1 = test_loc[[a]]
    }
    return(cl.pi.res.1)
  }))
  
  row.expr.sum = apply(abs( expr[,locs_test]-1),2,sum) 
  
  cnv.mat.p0=expr[,locs_test]
  
  cluster.res = Reduce(append,lapply(names_test, function(a){
    cl.pi.1 = infercnv_obj@tumor_subclusters[["hc"]][[a]]
    if(length(cl.pi.1$labels)>2){
      cl.pi.res.1 = cutree(cl.pi.1,k=min(c(10,length(cl.pi.1$labels))))
      cl.pi.res.1 = cl.pi.res.1[infercnv_obj@tumor_subclusters[["hc"]][[a]][["order"]]]
    }else{
      cl.pi.res.1 = rep(1,length(test_loc[[a]]))
    }
    cl.pi.res.1 = as.character(cl.pi.res.1)
    return(cl.pi.res.1)
  }))
  
  tissue.type = Reduce(append,lapply(names_test,function(a){
    rep(strsplit(a,'_')[[1]][2],
        length(test_loc[[a]]))
  }))
  annotation_row.i = data.frame(From = tissue.type,
                                Clusters = paste(tissue.type,'CC',
                                                 cluster.res),
                                Patient = pi,
                                row.names = colnames(cnv.mat.p0))
  
  cnv.mat.p0= as.data.frame(cnv.mat.p0)
  
  
  annotations.row=rbind(annotations.row,annotation_row.i)
  
  cluster.stat = data.frame(Clusters = paste(tissue.type,'CC',
                                             cluster.res),
                            score = row.expr.sum)
  
  cluster.stat.sum = as.data.frame(cluster.stat %>% group_by(Clusters) %>% summarise(mean=mean(score)))
  cluster.stat.sum$Patient = pi
  cluster.scores = rbind(cluster.scores,cluster.stat.sum)
  
  gene.sum.i = apply(abs( expr[,locs_test]),1,sum) 
  gene.sum.i = data.frame(gene=rownames(cnv.mat.p0),sum=gene.sum.i)
  gene.sum.i$Patient = pi
  
  gene.sum = rbind(gene.sum,gene.sum.i)
  
  # tiff(paste0("./documents/inferCNV/all/rebulid_cnv_all_",pi,".tiff"),width = 480,height = 540,units = 'px')
  # mat=cnv.mat.p0
  # mat[mat>=1.15]=1
  # mat[mat<0.85]=0.9
  # ht = pheatmap(t(mat),show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F,
  #             color = colorRampPalette(rev(hcl.colors(100,'RdBu')))(100),
  #             annotation_row = annotation_row.i,
  #             annotation_colors = list(Clusters=cl.cols,
  #                                      From=c('T'='red','L'='blue')))
  # grid::grid.newpage()
  # grid::grid.draw(ht$gtable)
  # dev.off()
  # print(paste('Finished for',pi))
  
}

write.csv(cluster.scores,file = './tables_v2/each inferCNV cluster scores.csv')
write.csv(gene.sum,file = './tables_v2/gene.sum.csv')
gene.sum = read.csv(file = './tables_v2/gene.sum.csv')
####Score of the genes####
gene.info = read.delim('/share/dichen/lungNodeM/documents/inferCNV/genes location-v2.txt',header = F)
rownames(gene.info)=gene.info$V1

gene.sum.top = group_by(gene.sum,Patient) %>% top_n(500,sum)

####LUSC
gene.sum.LUSC = gene.sum[gene.sum$Patient %in% c('P0','P1','P5','P9','P15'),]
gene.sum.top.LUSC = as.data.frame(gene.sum.top[gene.sum.top$Patient %in% c('P0','P1','P5','P9','P15'),])

genes.all.LUSC = unique(gene.sum.LUSC$gene)
gene.cnv.sum.LUSC = unlist(lapply(genes.all.LUSC,function(a){
  sum(gene.sum.LUSC[gene.sum.LUSC$gene == a,'sum'])
}))

gene.cnv.topN.LUSC = unlist(lapply(genes.all.LUSC,function(a){
  nrow(gene.sum.top.LUSC[gene.sum.top.LUSC$gene == a,])
}))


CNV.gene.res.LUSC = data.frame(gene=genes.all.LUSC,
                               sum.Score = gene.cnv.sum.LUSC,
                               sum.N = gene.cnv.topN.LUSC)
CNV.gene.res.LUSC$chr = gene.info[as.character(CNV.gene.res.LUSC$gene),'V2']

write.csv(CNV.gene.res.LUSC,file = './tables_v2/LUSC gene CNV score.csv')

####LUAD

gene.sum.LUAD = gene.sum[gene.sum$Patient %in% c('P0','P1','P5','P9','P15') == F,]
gene.sum.top.LUAD = as.data.frame(gene.sum.top[gene.sum.top$Patient %in% c('P0','P1','P5','P9','P15') == F,])

genes.all.LUAD = unique(gene.sum.LUAD$gene)
gene.cnv.sum.LUAD = unlist(lapply(genes.all.LUAD,function(a){
  sum(gene.sum.LUAD[gene.sum.LUAD$gene == a,'sum'])
}))

gene.cnv.topN.LUAD = unlist(lapply(genes.all.LUAD,function(a){
  nrow(gene.sum.top.LUAD[gene.sum.top.LUAD$gene == a,])
}))


CNV.gene.res.LUAD = data.frame(gene=genes.all.LUAD,
                               sum.Score = gene.cnv.sum.LUAD,
                               sum.N = gene.cnv.topN.LUAD)
CNV.gene.res.LUAD$chr = gene.info[as.character(CNV.gene.res.LUAD$gene),'V2']

write.csv(CNV.gene.res.LUAD,file = './tables_v2/LUAD gene CNV score.csv')


#####Check malignant or not####
annotations.row = data.frame()
cl.cols = cols.default[1:20]
#cl.cols = colorRampPalette(cl.cols)(40)
#names(cl.cols) = c(paste('mLN+ PT CC',1:20),paste('mLN CC',1:20))

names(cl.cols) = c(paste('mLN+ PT CC',1:10),paste('mLN CC',1:10))
cluster.scores = data.frame()
gene.sum = data.frame()

mat.list = list()

pi = 'P18' #From P5 20240429


if(pi %in% c('P3','P5','P6','P7')){
  
  infercnv_obj = readRDS(file = paste0('./documents/inferCNV_v2/each/',pi,'/run.final.infercnv_obj'))
}else{
  infercnv_obj = readRDS(file = paste0('./documents/inferCNV_v2/each2/',pi,'/run.final.infercnv_obj'))
  
}
expr = infercnv_obj@expr.data
test_loc <- infercnv_obj@observation_grouped_cell_indices

names_test = c( names(test_loc)[grepl('PT',names(test_loc))],
                names(test_loc)[grepl('mLN_',names(test_loc))])
locs_test = Reduce(union,lapply(names_test, function(a){
  cl.pi.1 = infercnv_obj@tumor_subclusters[["hc"]][[a]]
  if(length(cl.pi.1$labels)>2){
    cl.pi.res.1 =  test_loc[[a]][infercnv_obj@tumor_subclusters[["hc"]][[a]][["order"]]]
    
  }else{
    cl.pi.res.1 = test_loc[[a]]
  }
  return(cl.pi.res.1)
}))

row.expr.sum = apply(abs( expr[,locs_test]-1),2,sum) 

cnv.mat.p0=expr[,locs_test]

cluster.res = Reduce(append,lapply(names_test, function(a){
  cl.pi.1 = infercnv_obj@tumor_subclusters[["hc"]][[a]]
  if(length(cl.pi.1$labels)>2){
    cl.pi.res.1 = cutree(cl.pi.1,k=min(c(10,length(cl.pi.1$labels))))
    cl.pi.res.1 = cl.pi.res.1[infercnv_obj@tumor_subclusters[["hc"]][[a]][["order"]]]
  }else{
    cl.pi.res.1 = rep(1,length(test_loc[[a]]))
  }
  cl.pi.res.1 = as.character(cl.pi.res.1)
  return(cl.pi.res.1)
}))

tissue.type = Reduce(append,lapply(names_test,function(a){
  rep(strsplit(a,'_')[[1]][2],
      length(test_loc[[a]]))
}))
annotation_row.i = data.frame(From = tissue.type,
                              Clusters = paste(tissue.type,'CC',
                                               cluster.res),
                              Patient = pi,
                              row.names = colnames(cnv.mat.p0),
                              cellType = meta.data.PLM[colnames(cnv.mat.p0),'assigned_cell_type'])

cnv.mat.p0= as.data.frame(cnv.mat.p0)




cluster.stat = data.frame(Clusters = paste(tissue.type,'CC',
                                           cluster.res),
                          score = row.expr.sum)

cluster.stat.sum = as.data.frame(cluster.stat %>% group_by(Clusters) %>% summarise(mean=mean(score)))
cluster.stat.sum$Patient = pi


mat=cnv.mat.p0
mat[mat>=1.15]=1.15
mat[mat<0.85]=0.85
mat.list[[pi]]=mat
# cl.cols = cols.default[1:30]
# names(cl.cols) = c(paste('T CC',1:15),paste('L CC',1:15))

table(annotation_row.i$Clusters,annotation_row.i$cellType)
pheatmap(t(mat),show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(rev(hcl.colors(100,'RdBu')))(100),
         annotation_row = annotation_row.i[,-4],
         annotation_colors = list(Clusters=cl.cols,
                                  From=c('mLN+ PT'='red','mLN'='blue')))
gc()
normal.cls = c('mLN+ PT CC 7')
cluster.stat.sum$CNV.Type = ifelse(cluster.stat.sum$Clusters %in% normal.cls,'Normal','Malignant')
annotation_row.i$CNV.Type = ifelse(annotation_row.i$Clusters  %in% normal.cls,'Normal','Malignant')

cluster.scores = rbind(cluster.scores,cluster.stat.sum)
annotations.row=rbind(annotations.row,annotation_row.i)


write.csv(cluster.scores,file = './each inferCNV cluster scores.csv')
write.csv(annotations.row,file = './all inferCNV annotations of cells.csv')
saveRDS(mat.list,file='./inferCNV_v2_allMat.Rds')
table(annotations.row$Patient,annotations.row$CNV.Type)
annotations.row = read.csv(file = './all inferCNV annotations of cells.csv',row.names = 1)
####Merge all mat as one####
mat.list = readRDS('inferCNV_v2_allMat.Rds')
mat.1 = mat.list[[1]]
for(i in 2:length(mat.list)){
  mat.i = mat.list[[i]]
  genes.overlap = rownames(mat.1)[rownames(mat.1) %in% rownames(mat.i)]
  
  mat.1 = cbind(mat.1[genes.overlap,],mat.i[genes.overlap,])
}

dim(mat.1)

cols = colnames(mat.1)[sample(1:ncol(mat.1),15000)]
mat.1.sub = mat.1[,colnames(mat.1) %in% cols]

mat.1.sub[mat.1.sub>=1.05]=1.05
mat.1.sub[mat.1.sub<0.95]=0.95
cl.patients = cols.default[1:19]
names(cl.patients) = paste0('P',0:18)
cl.chrs = rep(c('grey','black'),11)
names(cl.chrs) = unique(gene.info[rownames(mat.1.sub),'V2'])
pheatmap(t(mat.1.sub),show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F,
         color = colorRampPalette(rev(hcl.colors(100,'RdBu')))(100),
         annotation_row = annotations.row[colnames(mat.1.sub),c('From','Patient','CNV.Type')],
         annotation_col = data.frame(chr = gene.info[rownames(mat.1.sub),'V2'],row.names = rownames(mat.1.sub)),
         annotation_colors = list(Patient=cl.patients,
                                  chr = cl.chrs,
                                  From=c('mLN+ PT'=as.character(pal_npg()(6))[2],'mLN'=as.character(pal_npg()(6))[3]),
                                  CNV.Type = c('Malignant'='red','Normal'='blue')))
