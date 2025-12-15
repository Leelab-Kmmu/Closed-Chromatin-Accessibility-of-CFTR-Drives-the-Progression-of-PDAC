###############################################################################################################################
###########################################Figure 5 ###########################################################################
###############################################################################################################################
library(Seurat)
#library(SeuratObject, lib.loc = "/home/wangmeiheng/R/library")
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(RColorBrewer)
library(clustree)
options(stringsAsFactors=F)
setwd('/400T/wangmeiheng/spatial_transcripts/A')
load('spatial_A_umap.Rdata')
setwd('/400T/wangmeiheng/spatial_transcripts/B')
load('spatial_B_umap.Rdata')
markerGenes  <- c(
  'KRT19','AMBP',"TM4SF1","S100A4",'MKI67',"MUC1",'EPCAM',
  'CFTR',
  'DCN','LUM','COL1A1',
  'CHGB',
  'PRSS1',
  'CD3D',
  'FCGR3A',
  'CD79A',
  'AIF1','CD1C',
  'CD68')
colors <- brewer.pal(n = 16, name = "Set3")
Idents(spatial_A) <- spatial_A$celltype
Idents(spatial_B) <- spatial_B$celltype
##################################################################################
###########################Figure S5 A Plot Starting##############################
##################################################################################
colors <- c("#2ECC71", "#E67E22","#1B4F72" , "#9B59B6", "#C0392B")
p1<-SpatialPlot(spatial_A, group.by = "celltype",image.alpha = 0)+
  scale_fill_manual(values = colors)#,alpha = 1,image.alpha = 0.2

p1
p2<-SpatialPlot(spatial_B, group.by = "celltype",image.alpha = 0)+
  scale_fill_manual(values = colors)#,alpha = 1,image.alpha = 0.2

p2
##################################################################################
###########################Figure S5 B Plot Starting##############################
##################################################################################
DotPlot(spatial_A, features = c(markerGenes),group.by = "celltype") + scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+RotatedAxis()
DotPlot(spatial_B, features = c(markerGenes),group.by = "celltype") + scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+RotatedAxis()
##################################################################################
###########################Figure S5 C Plot Starting##############################
##################################################################################
DoHeatmap(spatial_A, features = c(markerGenes), size = 3.5,group.colors = colors) + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill") + theme(axis.text.y = element_text(size = 9))
DoHeatmap(spatial_B, features = c(markerGenes), size = 3.5,group.colors = colors) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill") + theme(axis.text.y = element_text(size = 9))
##################################################################################
###########################Figure S4 A Plot Starting##############################
##################################################################################
colors <- c("#2166AC", "lightgrey", "#B2182B")
for (i in c("CFTR","KRT19")) {
  p<- SpatialFeaturePlot(spatial_A, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}
for (i in c("CFTR","KRT19")) {
  p<- SpatialFeaturePlot(spatial_B, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}
##################################################################################
###########################Figure S4 B Plot Starting##############################
##################################################################################
markerGenes<-c("CFTR","KRT19")
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_A")
files=list.files(pattern ='^GSM')
load(paste0(files[1],'/',files[1],"_seurat_ST.RData"))
A_GSM3036911<-seurat_ST
load(paste0(files[2],'/',files[2],"_seurat_ST.RData"))
A_GSM4100721<-seurat_ST
load(paste0(files[3],'/',files[3],"_seurat_ST.RData"))
A_GSM4100722<-seurat_ST
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B") 
files=list.files(pattern ='^GSM')
load(paste0(files[1],'/',files[1],"_seurat_ST.RData"))
B_GSM3405534<-seurat_ST
load(paste0(files[2],'/',files[2],"_seurat_ST.RData"))
B_GSM4100723<-seurat_ST
load(paste0(files[3],'/',files[3],"_seurat_ST.RData"))
B_GSM4100724<-seurat_ST
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_D")
files=list.files(pattern ='^GSM')
load(paste0(files[1],'/',files[1],"_seurat_ST.RData"))
D_GSM4100725<-seurat_ST
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_E")
files=list.files(pattern ='^GSM')
load(paste0(files[1],'/',files[1],"_seurat_ST.RData"))
E_GSM4100726<-seurat_ST
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_G")
files=list.files(pattern ='^GSM')
load(paste0(files[1],'/',files[1],"_seurat_ST.RData"))
G_GSM4100728<-seurat_ST
library(RColorBrewer)
library(Seurat)
library(SeuratObject)
library(ggplot2)
colors <- c("#2166AC", "lightgrey", "#B2182B")
for (i in markerGenes) {
  p<- SpatialFeaturePlot(A_GSM3036911, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}

for (i in markerGenes) {
  p<- SpatialFeaturePlot(A_GSM4100721, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}
for (i in markerGenes) {
  p<- SpatialFeaturePlot(A_GSM4100722, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}
for (i in markerGenes) {
  p<- SpatialFeaturePlot(B_GSM3405534, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}
for (i in markerGenes) {
  p<- SpatialFeaturePlot(B_GSM4100723, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}
for (i in markerGenes) {
  p<- SpatialFeaturePlot(B_GSM4100724, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}

for (i in markerGenes) {
  p<- SpatialFeaturePlot(D_GSM4100725, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}
for (i in markerGenes) {
  p<- SpatialFeaturePlot(E_GSM4100726, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}
for (i in markerGenes) {
  p<- SpatialFeaturePlot(G_GSM4100728, features = i, slot = "data",image.alpha = 0) + 
    scale_fill_gradientn(colors = colors) 
  print(p)
}


##################################################################################
###########################Figure S5 D Plot Starting##############################
##################################################################################
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_A") 
load('A_GSM3036911_umap.Rdata')
colors <- c("#1B4F72","#2ECC71", "#E67E22" , "#9B59B6")
Idents(A_GSM3036911) <- A_GSM3036911$celltype
p1<-SpatialPlot(A_GSM3036911, group.by = "celltype",image.alpha = 0)+
  scale_fill_manual(values = colors)#,alpha = 1,image.alpha = 0.2

p1
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B") 
load('B_GSM3405534_umap.Rdata')
Idents(B_GSM3405534) <- B_GSM3405534$celltype
p2<-SpatialPlot(B_GSM3405534, group.by = "celltype",image.alpha = 0)+
  scale_fill_manual(values = colors)#,alpha = 1,image.alpha = 0.2

p2

##################################################################################
###########################Figure 5 A Plot Starting###############################
##################################################################################
#######################data.merge
setwd("/400T/wangmeiheng/spatial_transcripts")
load("/400T/wangmeiheng/spatial_transcripts/data_merge_new.RData")
library(dplyr)
library(ggplot2)
# Retrieve expression data for CFTR and KRT19
cftr_expr <- FetchData(data.merge, vars = "CFTR")
krt19_expr <- FetchData(data.merge, vars = "KRT19")
threshold <- 0

data.merge$expression_status <- ifelse(cftr_expr > threshold & krt19_expr <= threshold, "CFTR exp, KRT19 no_exp", 
                                       ifelse(cftr_expr <= threshold & krt19_expr > threshold, "KRT19 exp, CFTR no_exp", 
                                              ifelse(cftr_expr > threshold & krt19_expr > threshold, "Both exp", 
                                                     "Both no_exp")))

data.merge$expression_status <- factor(data.merge$expression_status,levels = c("Both exp","Both no_exp","CFTR exp, KRT19 no_exp","KRT19 exp, CFTR no_exp"))
bar.df<-data.frame(data.merge$expression_status,names(data.merge$orig.ident),data.merge$orig.ident)

colnames(bar.df)<-c('expression_status','cell_ID','slice')


bar.df <- mutate(bar.df,name=factor(bar.df$slice))
bar.df$group <- factor(bar.df$slice, 
                       levels = c("spatial_A", "spatial_B"),
                       labels = c("spatial_A","spatial_B"))


text.df <- as.data.frame(table(bar.df$name)) 
color_cluster = c('Both exp'="#25ABAE",'Both no_exp'= "#E3837B",'CFTR exp, KRT19 no_exp'="#79B85B",'KRT19 exp, CFTR no_exp'='#C386B8')
SpatialDimPlot(data.merge, group.by = "expression_status",image.alpha = 0,
               pt.size.factor = 1.2,images = "image",cols=color_cluster,stroke = 0)
SpatialDimPlot(data.merge, group.by = "expression_status",image.alpha = 0,
               pt.size.factor = 1.2,images = "image.spatial_B",cols=color_cluster,stroke = 0)
ggplot(bar.df,aes(x=group))+
  geom_bar(aes(fill=expression_status),position = "fill",width = .7)+
  scale_x_discrete("")+
  scale_y_continuous("Total spot proportion",expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual("expression_status",values = color_cluster)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+ 
  scale_fill_manual("expression_status", 
                    values = color_cluster) +  
  coord_flip() 
save(data.merge,file = 'data_merge_new.RData')
library(dplyr)


##################################################################################
###########################Figure 5 B Plot Starting###############################
##################################################################################
##################################################################################
fs<-list.files('/400T/wangmeiheng/spatial_transcripts/',recursive = TRUE,full.names = TRUE,pattern = "GSM.*seurat_ST\\.RData$")
seurat_fs<-list()
setwd("/400T/wangmeiheng/spatial_transcripts/")
for (i in 1:length(fs)) {
  load(fs[i])
  seurat_ST<-RenameCells(seurat_ST,new.names=paste(Idents(seurat_ST),Cells(seurat_ST),sep = '_'))
  if ("KRT19" %in% rownames(seurat_ST))
  {seurat_fs[i]<-seurat_ST}
}
seurat_fs<-seurat_fs[-9]
names(seurat_fs)<-gsub(".*/(GSM\\d+).*", "\\1", fs[-9])
seurat_merge<-merge(seurat_fs[[1]],seurat_fs[-1])
load( 'GEO_ST_seurat_merge.RData')
#seurat_merge <- Reduce(function(x, y) merge(x, y), seurat_fs)
######################merge
library(dplyr)
library(ggplot2)
#Retrieve expression data for CFTR and KRT19
cftr_expr <- FetchData(seurat_merge, vars = "CFTR")
krt19_expr <- FetchData(seurat_merge, vars = "KRT19")
threshold <- 0

seurat_merge$expression_status <- ifelse(cftr_expr > threshold & krt19_expr <= threshold, "CFTR exp, KRT19 no_exp", 
                                         ifelse(cftr_expr <= threshold & krt19_expr > threshold, "KRT19 exp, CFTR no_exp", 
                                                ifelse(cftr_expr > threshold & krt19_expr > threshold, "Both exp", 
                                                       "Both no_exp")))

seurat_merge$expression_status <- factor(seurat_merge$expression_status,levels = c("Both exp","Both no_exp","CFTR exp, KRT19 no_exp","KRT19 exp, CFTR no_exp"))
bar.df<-data.frame(seurat_merge$expression_status,names(seurat_merge$orig.ident),seurat_merge$orig.ident)

colnames(bar.df)<-c('expression_status','cell_ID','slice')


bar.df <- mutate(bar.df,name=factor(bar.df$slice))
bar.df$group <- factor(bar.df$slice)


text.df <- as.data.frame(table(bar.df$name)) 
color_cluster = c('Both exp'="#25ABAE",'Both no_exp'= "#E3837B",'CFTR exp, KRT19 no_exp'="#79B85B",'KRT19 exp, CFTR no_exp'='#C386B8')

ggplot(bar.df,aes(x=group))+
  geom_bar(aes(fill=expression_status),position = "fill",width = .7)+
  scale_x_discrete("")+
  scale_y_continuous("Total spot proportion",expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual("expression_status",values = color_cluster)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+ 
  scale_fill_manual("expression_status", 
                    values = color_cluster) + 
  coord_flip()
save(seurat_merge,file = 'GEO_ST_seurat_merge.RData')
images<-names(seurat_merge@images)
for (i in 1:9) {
  p<-SpatialDimPlot(seurat_merge, group.by = "expression_status",image.alpha = 0,
                 pt.size.factor =3,images = images[i],cols=color_cluster)
  print(p)
}


##################################################################################
###########################Figure 5 C Plot Starting###############################
##################################################################################
#CFTR is expressed spatially_A in the DUC region(spatial A)
load("/400T/wangmeiheng/spatial_transcripts/A/spatial_A_umap.Rdata")
cftr_expr <- FetchData(spatial_A, vars = "CFTR")
threshold <- 0
spatial_A$CFTR_exp_region <- ifelse(cftr_expr > threshold &spatial_A$celltype!='DUC', "CFTR exp in no Duc region", 
                                    ifelse(cftr_expr <= threshold & spatial_A$celltype=='DUC', "CFTR no_exp in Duc region", 
                                           ifelse(cftr_expr > threshold & spatial_A$celltype=='DUC', "CFTR exp in Duc region", 
                                                  "CFTR no_exp in no Duc region")))

spatial_A$CFTR_exp_region <- factor(spatial_A$CFTR_exp_region,levels = c("CFTR exp in Duc region","CFTR no_exp in no Duc region","CFTR exp in no Duc region","CFTR no_exp in Duc region"))
table(spatial_A$CFTR_exp_region)
N <- sum(table(spatial_A$CFTR_exp_region))            
K <- length(which(spatial_A$CFTR_exp_region=='CFTR exp in Duc region')) + length(which(spatial_A$CFTR_exp_region=='CFTR exp in no Duc region'))         
M <- length(which(spatial_A$CFTR_exp_region=='CFTR exp in Duc region'))+ length(which(spatial_A$CFTR_exp_region=='CFTR no_exp in Duc region'))         
k <- length(which(spatial_A$CFTR_exp_region=='CFTR exp in Duc region'))              
#Hypergeometric test
p_value1 <- phyper(k - 1, M, N - M, K, lower.tail = F)
p_value_text <- paste0("p = ", formatC(p_value1, format = "e", digits = 2))
temp<-spatial_A$CFTR_exp_region
library(VennDiagram)
venn.plot <- draw.pairwise.venn(
  area1 = length(which(temp=='CFTR exp in Duc region')) + length(which(temp=='CFTR exp in Duc region')),   
  area2 = length(which(temp=='CFTR exp in Duc region')) + length(which(temp=='CFTR no_exp in Duc region')), 
  cross.area = length(which(temp=='CFTR exp in Duc region')),    
  category = c("CFTR exp", "Duc region"),
  fill = c('#E38AAE','#2B65EC'),  
  alpha = 0.6,                        
  lty = "blank",                      
  cex = 2,                          
  cat.cex = 2,                        
  cat.pos = c(-20, 20)                
)
grid.newpage()
grid.rect(gp = gpar(fill = "#FBE7A1", alpha = 0.8,col = NA))  
grid.draw(venn.plot)
grid.text(p_value_text, x = 0.5, y = 0.3, gp = gpar(fontsize = 16, col = "black"))
grid.text(N, x = 0.5, y = 0.9, gp = gpar(fontsize = 16, col = "black"))


#CFTR is expressed spatial_B in the DUC region(spatial B)
load("/400T/wangmeiheng/spatial_transcripts/B/spatial_B_umap.Rdata")
cftr_expr <- FetchData(spatial_B, vars = "CFTR")
threshold <- 0
spatial_B$CFTR_exp_region <- ifelse(cftr_expr > threshold &spatial_B$celltype!='DUC', "CFTR exp in no Duc region", 
                                    ifelse(cftr_expr <= threshold & spatial_B$celltype=='DUC', "CFTR no_exp in Duc region", 
                                           ifelse(cftr_expr > threshold & spatial_B$celltype=='DUC', "CFTR exp in Duc region", 
                                                  "CFTR no_exp in no Duc region")))

spatial_B$CFTR_exp_region <- factor(spatial_B$CFTR_exp_region,levels = c("CFTR exp in Duc region","CFTR no_exp in no Duc region","CFTR exp in no Duc region","CFTR no_exp in Duc region"))
table(spatial_B$CFTR_exp_region)
N <- sum(table(spatial_B$CFTR_exp_region))            
K <- length(which(spatial_B$CFTR_exp_region=='CFTR exp in Duc region')) + length(which(spatial_B$CFTR_exp_region=='CFTR exp in no Duc region'))         
M <- length(which(spatial_B$CFTR_exp_region=='CFTR exp in Duc region'))+ length(which(spatial_B$CFTR_exp_region=='CFTR no_exp in Duc region'))         
k <- length(which(spatial_B$CFTR_exp_region=='CFTR exp in Duc region'))              
#Hypergeometric test
p_value1 <- phyper(k - 1, M, N - M, K, lower.tail = F)
p_value_text <- paste0("p = ", formatC(p_value1, format = "e", digits = 2))
temp<-spatial_B$CFTR_exp_region
library(VennDiagram)
venn.plot <- draw.pairwise.venn(
  area1 = length(which(temp=='CFTR exp in Duc region')) + length(which(temp=='CFTR exp in Duc region')),   
  area2 = length(which(temp=='CFTR exp in Duc region')) + length(which(temp=='CFTR no_exp in Duc region')), 
  cross.area = length(which(temp=='CFTR exp in Duc region')),    
  category = c("CFTR exp", "Duc region"),
  fill = c('#E38AAE','#2B65EC'),  
  alpha = 0.6,                        
  lty = "blank",                      
  cex = 2,                          
  cat.cex = 2,                        
  cat.pos = c(-20, 20)                
)
grid.newpage()
grid.rect(gp = gpar(fill = "#FBE7A1", alpha = 0.8,col = NA))  
grid.draw(venn.plot)
grid.text(p_value_text, x = 0.5, y = 0.3, gp = gpar(fontsize = 16, col = "black"))
grid.text(N, x = 0.5, y = 0.9, gp = gpar(fontsize = 16, col = "black"))



###########################################################################################
#CFTR is expressed in the DUC region (A_GSM3036911)
load("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_A/A_GSM3036911_umap.Rdata")
#Obtain CFTR expression data
cftr_expr <- FetchData(A_GSM3036911, vars = "CFTR")
threshold <- 0

A_GSM3036911$CFTR_exp_region <- ifelse(cftr_expr > threshold &A_GSM3036911$celltype!='DUC', "CFTR exp in no Duc region", 
                                       ifelse(cftr_expr <= threshold & A_GSM3036911$celltype=='DUC', "CFTR no_exp in Duc region", 
                                              ifelse(cftr_expr > threshold & A_GSM3036911$celltype=='DUC', "CFTR exp in Duc region", 
                                                     "CFTR no_exp in no Duc region")))

A_GSM3036911$CFTR_exp_region <- factor(A_GSM3036911$CFTR_exp_region,levels = c("CFTR exp in Duc region","CFTR no_exp in no Duc region","CFTR exp in no Duc region","CFTR no_exp in Duc region"))
table(A_GSM3036911$CFTR_exp_region)
N <- sum(table(A_GSM3036911$CFTR_exp_region))            
K <- length(which(A_GSM3036911$CFTR_exp_region=='CFTR exp in Duc region')) + length(which(A_GSM3036911$CFTR_exp_region=='CFTR exp in no Duc region'))         
M <- length(which(A_GSM3036911$CFTR_exp_region=='CFTR exp in Duc region'))+ length(which(A_GSM3036911$CFTR_exp_region=='CFTR no_exp in Duc region'))         
k <- length(which(A_GSM3036911$CFTR_exp_region=='CFTR exp in Duc region'))                
# Hypergeometric test
p_value1 <- phyper(k - 1, M, N - M, K, lower.tail = F)
p_value_text <- paste0("p = ", formatC(p_value1, format = "e", digits = 2))
temp<-A_GSM3036911$CFTR_exp_region
library(VennDiagram)
venn.plot <- draw.pairwise.venn(
  area1 = length(which(temp=='CFTR exp in Duc region')) + length(which(temp=='CFTR exp in Duc region')),    
  area2 = length(which(temp=='CFTR exp in Duc region')) + length(which(temp=='CFTR no_exp in Duc region')),  
  cross.area = length(which(temp=='CFTR exp in Duc region')),     
  category = c("CFTR exp", "Duc region"),
  fill = c('#E38AAE','#2B65EC'), 
  alpha = 0.6,                       
  lty = "blank",                      
  cex = 2,                           
  cat.cex = 2,                         
  cat.pos = c(-20, 20)                
)
grid.newpage()
grid.rect(gp = gpar(fill = "#FBE7A1", alpha = 0.8,col = NA))  
grid.draw(venn.plot)
grid.text(p_value_text, x = 0.5, y = 0.3, gp = gpar(fontsize = 16, col = "black"))
grid.text(N, x = 0.5, y = 0.9, gp = gpar(fontsize = 16, col = "black"))
#CFTR is expressed in the DUC region (B_GSM3405534)
load("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B/B_GSM3405534_umap.Rdata")
cftr_expr <- FetchData(B_GSM3405534, vars = "CFTR")
threshold <- 0
B_GSM3405534$CFTR_exp_region <- ifelse(cftr_expr > threshold &B_GSM3405534$celltype!='DUC', "CFTR exp in no Duc region", 
                                       ifelse(cftr_expr <= threshold & B_GSM3405534$celltype=='DUC', "CFTR no_exp in Duc region", 
                                              ifelse(cftr_expr > threshold & B_GSM3405534$celltype=='DUC', "CFTR exp in Duc region", 
                                                     "CFTR no_exp in no Duc region")))

B_GSM3405534$CFTR_exp_region <- factor(B_GSM3405534$CFTR_exp_region,levels = c("CFTR exp in Duc region","CFTR no_exp in no Duc region","CFTR exp in no Duc region","CFTR no_exp in Duc region"))
table(B_GSM3405534$CFTR_exp_region)
N <- sum(table(B_GSM3405534$CFTR_exp_region))            
K <- length(which(B_GSM3405534$CFTR_exp_region=='CFTR exp in Duc region')) + length(which(B_GSM3405534$CFTR_exp_region=='CFTR exp in no Duc region'))         
M <- length(which(B_GSM3405534$CFTR_exp_region=='CFTR exp in Duc region'))+ length(which(B_GSM3405534$CFTR_exp_region=='CFTR no_exp in Duc region'))        
k <- length(which(B_GSM3405534$CFTR_exp_region=='CFTR exp in Duc region'))               
# Hypergeometric test
p_value1 <- phyper(k - 1, M, N - M, K, lower.tail = F)
p_value_text <- paste0("p = ", formatC(p_value1, format = "e", digits = 2))
temp<-B_GSM3405534$CFTR_exp_region
library(VennDiagram)
venn.plot <- draw.pairwise.venn(
  area1 = length(which(temp=='CFTR exp in Duc region')) + length(which(temp=='CFTR exp in Duc region')),   
  area2 = length(which(temp=='CFTR exp in Duc region')) + length(which(temp=='CFTR no_exp in Duc region')),   
  cross.area = length(which(temp=='CFTR exp in Duc region')),    
  category = c("CFTR exp", "Duc region"),
  fill = c('#E38AAE','#2B65EC'), 
  alpha = 0.6,                       
  lty = "blank",                     
  cex = 2,                             
  cat.cex = 2,                         
  cat.pos = c(-20, 20)                
)
grid.newpage()
grid.rect(gp = gpar(fill = "#FBE7A1", alpha = 0.8,col = NA))  
grid.draw(venn.plot)
grid.text(p_value_text, x = 0.5, y = 0.3, gp = gpar(fontsize = 16, col = "black"))
grid.text(N, x = 0.5, y = 0.9, gp = gpar(fontsize = 16, col = "black"))

##################################################################################
###########################Figure 5 D Plot Starting###############################
##################################################################################
#CFTR is not expressed in the malignant region (spatial_A)
load("/400T/wangmeiheng/spatial_transcripts/A/spatial_A_umap.Rdata")
threshold <- 0
cftr_expr <- FetchData(spatial_A, vars = "CFTR")
spatial_A$CFTR_exp_cancer <- ifelse(cftr_expr > threshold & spatial_A$celltype!="Malignant DUC", "CFTR exp in no Cancer region", 
                                    ifelse(cftr_expr <= threshold & spatial_A$celltype=="Malignant DUC", "CFTR no_exp in Cancer region", 
                                           ifelse(cftr_expr > threshold & spatial_A$celltype=='Malignant DUC', "CFTR exp in Cancer region", 
                                                  "CFTR no_exp in no Cancer region")))

spatial_A$CFTR_exp_cancer <- factor(spatial_A$CFTR_exp_cancer,levels = c("CFTR exp in Cancer region","CFTR no_exp in no Cancer region","CFTR exp in no Cancer region","CFTR no_exp in Cancer region"))
table(spatial_A$CFTR_exp_cancer)
N <- sum(table(spatial_A$CFTR_exp_cancer))            
K <- length(which(spatial_A$CFTR_exp_cancer=='CFTR no_exp in Cancer region')) + length(which(spatial_A$CFTR_exp_cancer=='CFTR no_exp in no Cancer region'))        
M <- length(which(spatial_A$CFTR_exp_cancer=='CFTR exp in Cancer region'))+ length(which(spatial_A$CFTR_exp_cancer=='CFTR no_exp in Cancer region'))        
k <- length(which(spatial_A$CFTR_exp_cancer=='CFTR no_exp in Cancer region'))               
# Hypergeometric test
p_value2 <- phyper(k - 1, M, N - M, K, lower.tail = F)
temp<-spatial_A$CFTR_exp_cancer
library(VennDiagram)
p_value_text <- paste0("p = ", formatC(p_value2, format = "e", digits = 2))

venn.plot <- draw.pairwise.venn(
  area1 = length(which(temp=='CFTR no_exp in Cancer region')) + length(which(temp=='CFTR no_exp in no Cancer region')),    
  area2 = length(which(temp=='CFTR no_exp in Cancer region')) + length(which(temp=='CFTR exp in Cancer region')),  
  cross.area = length(which(temp=='CFTR no_exp in Cancer region')), 
  category = c("CFTR no_exp", "cancer region"),
  fill = c('#E38AAE','#2B65EC'),
  alpha = 0.6,                     
  lty = "blank",                       
  cex = 2,                            
  cat.cex = 2,                         
  cat.pos = c(-20, 20)                 
)
grid.newpage()
grid.rect(gp = gpar(fill = "#FBE7A1", alpha = 0.8,col = NA))  
grid.draw(venn.plot)
grid.text(p_value_text, x = 0.5, y = 0.3, gp = gpar(fontsize = 16, col = "black"))
grid.text(N, x = 0.5, y = 0.9, gp = gpar(fontsize = 16, col = "black"))


#CFTR is not expressed in the malignant region (spatial_B)
load("/400T/wangmeiheng/spatial_transcripts/B/spatial_B_umap.Rdata")
threshold <- 0
cftr_expr <- FetchData(spatial_B, vars = "CFTR")
spatial_B$CFTR_exp_cancer <- ifelse(cftr_expr > threshold & spatial_B$celltype!="Malignant DUC", "CFTR exp in no Cancer region", 
                                    ifelse(cftr_expr <= threshold & spatial_B$celltype=="Malignant DUC", "CFTR no_exp in Cancer region", 
                                           ifelse(cftr_expr > threshold & spatial_B$celltype=='Malignant DUC', "CFTR exp in Cancer region", 
                                                  "CFTR no_exp in no Cancer region")))

spatial_B$CFTR_exp_cancer <- factor(spatial_B$CFTR_exp_cancer,levels = c("CFTR exp in Cancer region","CFTR no_exp in no Cancer region","CFTR exp in no Cancer region","CFTR no_exp in Cancer region"))
table(spatial_B$CFTR_exp_cancer)
N <- sum(table(spatial_B$CFTR_exp_cancer))            
K <- length(which(spatial_B$CFTR_exp_cancer=='CFTR no_exp in Cancer region')) + length(which(spatial_B$CFTR_exp_cancer=='CFTR no_exp in no Cancer region'))        
M <- length(which(spatial_B$CFTR_exp_cancer=='CFTR exp in Cancer region'))+ length(which(spatial_B$CFTR_exp_cancer=='CFTR no_exp in Cancer region'))        
k <- length(which(spatial_B$CFTR_exp_cancer=='CFTR no_exp in Cancer region'))               
# Hypergeometric test
p_value2 <- phyper(k - 1, M, N - M, K, lower.tail = F)
temp<-spatial_B$CFTR_exp_cancer
library(VennDiagram)
p_value_text <- paste0("p = ", formatC(p_value2, format = "e", digits = 2))

venn.plot <- draw.pairwise.venn(
  area1 = length(which(temp=='CFTR no_exp in Cancer region')) + length(which(temp=='CFTR no_exp in no Cancer region')),    
  area2 = length(which(temp=='CFTR no_exp in Cancer region')) + length(which(temp=='CFTR exp in Cancer region')),  
  cross.area = length(which(temp=='CFTR no_exp in Cancer region')), 
  category = c("CFTR no_exp", "cancer region"),
  fill = c('#E38AAE','#2B65EC'),
  alpha = 0.6,                     
  lty = "blank",                       
  cex = 2,                            
  cat.cex = 2,                         
  cat.pos = c(-20, 20)                 
)
grid.newpage()
grid.rect(gp = gpar(fill = "#FBE7A1", alpha = 0.8,col = NA))  
grid.draw(venn.plot)
grid.text(p_value_text, x = 0.5, y = 0.3, gp = gpar(fontsize = 16, col = "black"))
grid.text(N, x = 0.5, y = 0.9, gp = gpar(fontsize = 16, col = "black"))

###############################################################################
#CFTR is not expressed in the malignant region (A_GSM3036911)
load("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_A/A_GSM3036911_umap.Rdata")
threshold <- 0
cftr_expr <- FetchData(A_GSM3036911, vars = "CFTR")
A_GSM3036911$CFTR_exp_cancer <- ifelse(cftr_expr > threshold & A_GSM3036911$celltype!="Cancer region", "CFTR exp in no Cancer region", 
                                       ifelse(cftr_expr <= threshold & A_GSM3036911$celltype=="Cancer region", "CFTR no_exp in Cancer region", 
                                              ifelse(cftr_expr > threshold & A_GSM3036911$celltype=='Cancer region', "CFTR exp in Cancer region", 
                                                     "CFTR no_exp in no Cancer region")))

A_GSM3036911$CFTR_exp_cancer <- factor(A_GSM3036911$CFTR_exp_cancer,levels = c("CFTR exp in Cancer region","CFTR no_exp in no Cancer region","CFTR exp in no Cancer region","CFTR no_exp in Cancer region"))
table(A_GSM3036911$CFTR_exp_cancer)
N <- sum(table(A_GSM3036911$CFTR_exp_cancer))            
K <- length(which(A_GSM3036911$CFTR_exp_cancer=='CFTR no_exp in Cancer region')) + length(which(A_GSM3036911$CFTR_exp_cancer=='CFTR no_exp in no Cancer region'))         
M <- length(which(A_GSM3036911$CFTR_exp_cancer=='CFTR exp in Cancer region'))+ length(which(A_GSM3036911$CFTR_exp_cancer=='CFTR no_exp in Cancer region'))        
k <- length(which(A_GSM3036911$CFTR_exp_cancer=='CFTR no_exp in Cancer region'))                
# Hypergeometric test
p_value2 <- phyper(k - 1, M, N - M, K, lower.tail = F)
p_value_text <- paste0("p = ", formatC(p_value2, format = "e", digits = 2))
temp<-A_GSM3036911$CFTR_exp_cancer
library(VennDiagram)
venn.plot <- draw.pairwise.venn(
  area1 = length(which(temp=='CFTR no_exp in Cancer region')) + length(which(temp=='CFTR no_exp in no Cancer region')),    
  area2 = length(which(temp=='CFTR no_exp in Cancer region')) + length(which(temp=='CFTR exp in Cancer region')),  
  cross.area = length(which(temp=='CFTR no_exp in Cancer region')),    
  category = c("CFTR no_exp", "cancer region"),
  fill = c('#E38AAE','#2B65EC'), 
  alpha = 0.6,                         
  lty = "blank",                      
  cex = 2,                            
  cat.cex = 2,                        
  cat.pos = c(-20, 20)                
)
grid.newpage()
grid.rect(gp = gpar(fill = "#FBE7A1", alpha = 0.8,col = NA))  
grid.draw(venn.plot)
grid.text(p_value_text, x = 0.5, y = 0.3, gp = gpar(fontsize = 16, col = "black"))
grid.text(N, x = 0.5, y = 0.9, gp = gpar(fontsize = 16, col = "black"))


#CFTR is not expressed in the malignant region (B_GSM3405534)
load("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B/B_GSM3405534_umap.Rdata")
cftr_expr <- FetchData(B_GSM3405534, vars = "CFTR")
threshold <- 0
B_GSM3405534$CFTR_exp_cancer <- ifelse(cftr_expr > threshold & B_GSM3405534$celltype!="Cancer region", "CFTR exp in no Cancer region", 
                                       ifelse(cftr_expr <= threshold & B_GSM3405534$celltype=="Cancer region", "CFTR no_exp in Cancer region", 
                                              ifelse(cftr_expr > threshold & B_GSM3405534$celltype=='Cancer region', "CFTR exp in Cancer region", 
                                                     "CFTR no_exp in no Cancer region")))

B_GSM3405534$CFTR_exp_cancer <- factor(B_GSM3405534$CFTR_exp_cancer,levels = c("CFTR exp in Cancer region","CFTR no_exp in no Cancer region","CFTR exp in no Cancer region","CFTR no_exp in Cancer region"))
table(B_GSM3405534$CFTR_exp_cancer)
N <- sum(table(B_GSM3405534$CFTR_exp_cancer))            
K <- length(which(B_GSM3405534$CFTR_exp_cancer=='CFTR no_exp in Cancer region')) + length(which(B_GSM3405534$CFTR_exp_cancer=='CFTR no_exp in no Cancer region'))          
M <- length(which(B_GSM3405534$CFTR_exp_cancer=='CFTR exp in Cancer region'))+ length(which(B_GSM3405534$CFTR_exp_cancer=='CFTR no_exp in Cancer region'))         
k <- length(which(B_GSM3405534$CFTR_exp_cancer=='CFTR no_exp in Cancer region'))          
# Hypergeometric test
p_value2 <- phyper(k - 1, M, N - M, K, lower.tail = F)
p_value_text <- paste0("p = ", formatC(p_value2, format = "e", digits = 2))
temp<-B_GSM3405534$CFTR_exp_cancer
library(VennDiagram)
venn.plot <- draw.pairwise.venn(
  area1 = length(which(temp=='CFTR no_exp in Cancer region')) + length(which(temp=='CFTR no_exp in no Cancer region')),   
  area2 = length(which(temp=='CFTR no_exp in Cancer region')) + length(which(temp=='CFTR exp in Cancer region')),  
  cross.area = length(which(temp=='CFTR no_exp in Cancer region')),    
  category = c("CFTR no_exp", "cancer region"),
  fill = c('#E38AAE','#2B65EC'),  
  alpha = 0.6,                        
  lty = "blank",                       
  cex = 2,                            
  cat.cex = 2,                        
  cat.pos = c(-20, 20)                 
)
grid.newpage()
grid.rect(gp = gpar(fill = "#FBE7A1", alpha = 0.8,col = NA))  
grid.draw(venn.plot)
grid.text(p_value_text, x = 0.5, y = 0.3, gp = gpar(fontsize = 16, col = "black"))
grid.text(N, x = 0.5, y = 0.9, gp = gpar(fontsize = 16, col = "black"))
