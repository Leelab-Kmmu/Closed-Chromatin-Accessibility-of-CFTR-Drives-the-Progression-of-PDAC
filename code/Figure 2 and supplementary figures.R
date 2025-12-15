############Figure 1
library(ArchR)
##################################################################################
###########################Figure S1A Plot Starting###############################
##################################################################################
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
proj@cellColData@listData[["log10(nFrags)"]]<-log10(proj@cellColData@listData[["nFrags"]])
df<-getCellColData(proj,select =c("log10(nFrags)","TSSEnrichment"))

p<-ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99))
)
p + geom_hline(yintercept =4,lty= "dashed")+ geom_vline(xintercept = 3, ity = "dashed")
#violin
plotGroups(
  ArchRProj = proj,
  groupBy= "Sample",
  colorBy ="cellcolData",
  name = "BlacklistRatio",
  plotAs ="violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  palette = c("#88C5EC", "#25ABAE", "#6F4598", "#C386B8", "#D1669B", "#EBA05F", "#79BEA5", "#D1669B")
)
plotGroups(
  ArchRProj = proj,
  groupBy= "Sample",
  colorBy ="cellcolData",
  name = "NucleosomeRatio",
  plotAs ="violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  palette = c("#88C5EC", "#25ABAE", "#6F4598", "#C386B8", "#D1669B", "#EBA05F", "#79BEA5", "#D1669B")
)
plotFragmentSizes(ArchRProj=proj)
plotTSSEnrichment(ArchRProj=proj)

#Dimension Reduction and Clustering
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix",name = "IterativeLSI",varFeatures = 30000,clusterParams = list(resolution = c(10), sampleCells = 10000, maxClusters = 20, n.start= 10),force = TRUE)
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI",maxClusters = 50,force = TRUE)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force = TRUE)
#proj_umap <- saveArchRProject(ArchRProj = proj)
##################################################################################
###########################Figure S1B Plot Starting###############################
##################################################################################
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample",embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p2, type = "h")#Cluster Umap
##################################################################################
###########################Figure 2B (Top right) Plot Starting#########################
##################################################################################
ggAlignPlots(p1, type = "h")#Sample Umap

##################################################################################
###########################Figure S1B (Top right) Plot Starting########################
##################################################################################
temp_type<-proj@cellColData@rownames
temp_type<-unlist(lapply(strsplit(temp_type,"#"),function(x){return(x[[1]])}))
temp_type[temp_type=='ATAC1']<-'Tumor'
temp_type[temp_type=='ATAC2']<-'Tumor'
temp_type[temp_type=='ATAC3']<-'Tumor'
temp_type[temp_type=='ATAC_n']<-'Normal'
temp_type[temp_type=='ATAC_p1']<-'Normal'
temp_type[temp_type=='ATAC_p2']<-'Normal'
temp_type[temp_type=='ATAC_p3']<-'Normal'
temp_type[temp_type=='ATAC_p4']<-'Normal'
proj@cellColData@listData$sampletype<-temp_type
p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "sampletype", embedding = "UMAP",
                    pal=c('Normal'='#447CAA','Tumor'='#841B1F'))+
  theme(
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14) 
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 6))  
  )
ggAlignPlots(p4, type = "h")

# Cell annotation
proj@cellColData@listData$celltype<-proj@cellColData@listData$Clusters
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C1']<-'B'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C2']<-'B'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C3']<-'B'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C4']<-'T'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C5']<-'T'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C6']<-'NKT'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C7']<-'T'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C8']<-'ENDOCRINE'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C9']<-'DUC_n'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C10']<-'ACINAR'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C11']<-'ACINAR'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C12']<-'DUC'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C13']<-'DUC'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C14']<-'DUC_n'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C15']<-'DC'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C16']<-'MAC'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C17']<-'FIBROBLAST'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C18']<-'ENDOTHELIAL'
proj@cellColData@listData$celltype[proj@cellColData@listData$celltype=='C19']<-'FIBROBLAST'

##################################################################################
###########################Figure 2B Plot Starting################################
##################################################################################
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "celltype", embedding = "UMAP",
                    pal=c('B'="#25ABAE",'T'= "#E3837B",'NKT'="#D1669B",'ENDOCRINE'='#E06AA4','DUC_n'="#C386B8",'ACINAR'="#88C5EC",'DUC'="#6F4598",
                          'DC'="#999B24",'MAC'="#79B85B",'FIBROBLAST'="#79BEA5",'ENDOTHELIAL'="#EBA05F")) +
  theme(
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14)  
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 6)) 
  )
p3

# marker gene
markerGenes  <- c('CD3D',
                  'FCGR3A',
                  'CD19','MS4A1','CD79A','CD79B','IGLL5',
                  'KRT19','MKI67',
                  'CDH5',
                  'DCN','LUM',
                  'CHGB',
                  'PRSS1',
                  'EPCAM',
                  'AIF1','ITGAX','CD1C',
                  'CD68')
##################################################################################
###########################Figure 2C Plot Starting################################
##################################################################################
#Marker gene dot plot drawing
detach("package:plyr", unload=TRUE)
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
colorMat <- ArchR:::.getMatrixValues(ArchRProj = proj, name = markerGenes,matrixName = "GeneScoreMatrix", log2Norm = FALSE)
temp_mat<-data.frame(t(colorMat))

cell_annotations <- data.frame(
  Cell = rownames(temp_mat),
  CellType = proj@cellColData@listData[["celltype"]],
  Cluster = proj@cellColData@listData[["Clusters"]]
)
long_data <- temp_mat %>%
  rownames_to_column(var = "Cell") %>%
  pivot_longer(cols = -Cell, names_to = "Gene", values_to = "Expression") %>%
  left_join(cell_annotations, by = "Cell")

# Calculate the average normalised expression values and expression ratios
gene_data_avg <- long_data %>%
  group_by(Gene, CellType) %>%
  summarise(
    AverageStdExpression = mean(Expression, na.rm = TRUE),
    CellCount = dplyr::n(),  
    
    ExpressionProportion = sum(Expression > 0) / length(Expression) * 100, 
    .groups = 'drop'
  )
# z-score
gene_data_avg <- gene_data_avg %>%
  group_by(Gene) %>%
  mutate(StdExpression = (AverageStdExpression - mean(AverageStdExpression, na.rm = TRUE)) / sd(AverageStdExpression, na.rm = TRUE)) %>%
  ungroup()
#Set Gene as a factor and arrange it according to the order of markerGenes
gene_data_avg$Gene <- factor(gene_data_avg$Gene, levels = markerGenes)
gene_data_avg$CellType<-factor(gene_data_avg$CellType,levels =c('T','NKT','B','DUC','ENDOTHELIAL','FIBROBLAST','ENDOCRINE','ACINAR','DUC_n','DC','MAC'))
p1 <- ggplot(gene_data_avg, aes(x = Gene, y = CellType)) +
  geom_point(aes(size = ExpressionProportion, color = StdExpression)) + 
  scale_color_gradient2(low = "steelblue", mid = "lightgrey", high = "darkgoldenrod1", midpoint = 1) +  
  theme_minimal() +  
  labs(size = "ExpressionProportion", color = "Gene accessibility score") +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid = element_blank(), 
    panel.border = element_blank(), 
    axis.line = element_line() 
  )
p1
##################################################################################
###########################Figure S1C Plot Starting###############################
##################################################################################
gene_data_avg <- long_data %>%
  group_by(Gene, Cluster) %>%
  summarise(
    AverageStdExpression = mean(Expression, na.rm = TRUE),  
    CellCount = dplyr::n(), 
    ExpressionProportion = sum(Expression > 0) / length(Expression) * 100,  #
    .groups = 'drop'
  )
#z-score 
gene_data_avg <- gene_data_avg %>%
  group_by(Gene) %>%
  mutate(StdExpression = (AverageStdExpression - mean(AverageStdExpression, na.rm = TRUE)) / sd(AverageStdExpression, na.rm = TRUE)) %>%
  ungroup()

#ploting
#Set Gene as a factor and arrange it according to the order of markerGenes 
gene_data_avg$Gene <- factor(gene_data_avg$Gene, levels = markerGenes)
gene_data_avg$Cluster<-factor(gene_data_avg$Cluster,levels = c('C4','C5','C7','C6','C1','C2','C3','C12','C13','C18','C17','C19',                                                               'C8','C10','C11','C9','C14','C15','C16'))
p2 <- ggplot(gene_data_avg, aes(x = Gene, y = Cluster)) +
  geom_point(aes(size = ExpressionProportion, color = StdExpression)) +  
  scale_color_gradient2(low = "steelblue", mid = "lightgrey", high = "darkgoldenrod1", midpoint = 1) + 
  theme_minimal() +  
  labs(size = "ExpressionProportion", color = "Gene accessibility score") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    panel.grid = element_blank(),  
    panel.border = element_blank(),  
    axis.line = element_line()  
  )
p2
##################################################################################
###########################Figure 2D Plot Starting################################
##################################################################################
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
proj <- addImputeWeights(proj)
plotEmbedding(
  proj, 
  colorBy = "GeneScoreMatrix", 
  name = "CFTR", 
  embedding = "UMAP",
  pal = paletteContinuous("solarExtra")
)
##################################################################################
###########################Figure 2E Plot Starting################################
##################################################################################
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
proj <- addImputeWeights(proj)
library(presto)
library("hexbin")
colorMat <- ArchR:::.getMatrixValues(ArchRProj = proj, name = markerGenes,matrixName = "GeneScoreMatrix", log2Norm = FALSE)
colorMat <- imputeMatrix(mat = as.matrix(colorMat),imputeWeights = getImputeWeights(proj))
#save(colorMat,file='/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/SPOT_test/colorMat.Rdata')

load('/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/SPOT_test/colorMat.Rdata')
#Boxplot of CFTR gene score in different cell types
marker_TSS<-colorMat
CFTR_tss_150_data<-marker_TSS['CFTR',]
df<-as.data.frame(CFTR_tss_150_data)
colnames(df)<-c("CFTR")
df$type<-proj@cellColData@listData[["celltype"]]
df_melt<-reshape::melt(df,id.vars=c("type"))
df_melt$combine_name<-factor(paste(df_melt$variable,df_melt$type,sep="_"))
df_melt$combine_name <- factor(df_melt$combine_name,levels = c("CFTR_DUC_n", "CFTR_DUC", "CFTR_FIBROBLAST", "CFTR_T", "CFTR_NKT", "CFTR_MAC",
                                                               'CFTR_ENDOTHELIAL','CFTR_ENDOCRINE','CFTR_DC','CFTR_B','CFTR_ACINAR'))
color = c('B'="#25ABAE",'T'= "#E3837B",'NKT'="#D1669B",'ENDOCRINE'='#E06AA4','DUC_n'="#C386B8",'ACINAR'="#88C5EC",'DUC'="#6F4598",
          'DC'="#999B24",'MAC'="#79B85B",'FIBROBLAST'="#79BEA5",'ENDOTHELIAL'="#EBA05F")
#plot 
comparisons<- list(c("CFTR_DUC_n","CFTR_DUC"),
                   c("CFTR_DUC_n","CFTR_FIBROBLAST"),
                   c("CFTR_DUC_n","CFTR_T"),
                   c("CFTR_DUC_n","CFTR_NKT"),
                   c("CFTR_DUC_n","CFTR_MAC"),
                   c("CFTR_DUC_n","CFTR_ENDOTHELIAL"),
                   c("CFTR_DUC_n","CFTR_ENDOCRINE"),
                   c("CFTR_DUC_n","CFTR_DC"),
                   c("CFTR_DUC_n","CFTR_B"),
                   c("CFTR_DUC_n","CFTR_ACINAR"))

library(ggplot2)
library(ggpubr)

p <- df_melt %>% ggplot(aes(x=combine_name, y=value, fill=type)) +
  geom_boxplot(width=0.5, position=position_dodge(width=1), outlier.colour = NA, color="black", size=0.3) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12), 
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = rel(1)),  
    panel.border = element_rect(color = "black", size = 1)  
  ) +
  scale_fill_manual(values = color) +
  labs(y = "Gene accessibility score") +  
  scale_x_discrete(labels = c(
    "CFTR_DUC_n" = "Normal ductal", "CFTR_DUC" = "Tumor ductal", "CFTR_FIBROBLAST" = "Fibroblast", 
    "CFTR_T" = "T", "CFTR_NKT" = "NKT", "CFTR_MAC" = "Macrophage",
    "CFTR_ENDOTHELIAL" = "Endothelial", "CFTR_ENDOCRINE" = "Endocrine",
    "CFTR_DC" = "DC", "CFTR_B" = "B", "CFTR_ACINAR" = "Acinar")) 
p+stat_compare_means(comparisons = comparisons,label= "p.signif")
##################################################################################
###########################Figure 2F Plot Starting################################
##################################################################################
#Box plot--healthy,PDAC,STAD,other cancers, Noncancer,and gastritis
# Import data
setwd('/400T/wangmeiheng/cfDNA/tss_coverage')
load('PDAC_sample_tss_150.Rdata')#39
colnames(temp_tss_150)[1:3]<-c('ATAC_1836591','ATAC_T001835227','ATAC_T001837519')
new_temp_tss_150<-temp_tss_150
load('healthy_sample_tss_150.Rdata')#100
new_temp_tss_150<-cbind(new_temp_tss_150,temp_tss_150)
load('gastric_sample_tss_150.Rdata')#101
new_temp_tss_150<-cbind(new_temp_tss_150,temp_tss_150)
load('TGY_sample_tss_150.Rdata')#38
new_temp_tss_150<-cbind(new_temp_tss_150,temp_tss_150)
load('NGY_sample_tss_150.Rdata')#58
new_temp_tss_150<-cbind(new_temp_tss_150,temp_tss_150)
load('GGY_sample_tss_150.Rdata')#17
new_temp_tss_150<-cbind(new_temp_tss_150,temp_tss_150)

temp_tss_150<-new_temp_tss_150#353
temp_tss_150[is.na(temp_tss_150)]<-0.01
temp_tss_150[temp_tss_150==0]<-0.01
tss_150_file<-temp_tss_150
#save(tss_150_file,file='tss_150_sample353.Rdata')
load('tss_pro_table.Rdata')#70641     8
load('tss_150_sample353.Rdata')

CFTR_5<-tss_pro_table[which(tss_pro_table$gene_name=="CFTR"),]
#Ploting
pdf("CFTR_5_cfDNA_TSScoverage.pdf",width = 8, height = 8)
for (i in 2) {
  CFTR_tss_150_data<-tss_150_file[which(tss_pro_table$gene_name=="CFTR")[i],]
  df<-as.data.frame(t(CFTR_tss_150_data))
  colnames(df)<-c("CFTR")
  df$type<-c(rep("PDAC",39),rep("Normal",100),rep("gastric",101),rep("TGY",38),rep("NGY",58),rep("GGY",17))
  df_melt<-reshape::melt(df,id.vars=c("type"))
  df_melt$combine_name<-factor(paste(df_melt$variable,df_melt$type,sep="_"))
  df_melt$combine_name <- factor(df_melt$combine_name, 
                                 levels = c("CFTR_Normal","CFTR_PDAC", "CFTR_gastric", "CFTR_TGY", "CFTR_NGY", "CFTR_GGY"))
  color = c("#999B24", "#25ABAE", "#C386B8","#88C5EC","#6F4598", "#EBA05F" )
  #plot 
  comparisons<- list(c("CFTR_PDAC","CFTR_Normal"),
                     c("CFTR_PDAC","CFTR_gastric"),
                     c("CFTR_PDAC","CFTR_TGY"),
                     c("CFTR_PDAC","CFTR_NGY"),
                     c("CFTR_PDAC","CFTR_GGY"))
  
  library(ggplot2)
  library(ggpubr)
  
  
  p<- df_melt %>% ggplot(aes(x=combine_name,y=-log2(value),fill=type))+
    geom_boxplot(width=0.5,position=position_dodge(width=1),outlier.colour = NA)+
    theme_bw()+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),  
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 12),  
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(1)),  
          panel.border = element_rect(color = "black", size = 1)
    )+scale_fill_manual(values = color)
  
  p1<- p +stat_compare_means(comparisons = comparisons[1], label.y =3,method="wilcox.test",size = 3.5,label= "p.signif")
  p2<- p1 +stat_compare_means(comparisons = comparisons[2], label.y =3.5,method="wilcox.test",size = 3.5,label= "p.signif")
  p3<- p2 +stat_compare_means(comparisons = comparisons[3], label.y =4,method="wilcox.test",size = 3.5,label= "p.signif")
  p4<- p3 +stat_compare_means(comparisons = comparisons[4], label.y =4.5,method="wilcox.test",size = 3.5,label= "p.signif")
  p5<- p4 +stat_compare_means(comparisons = comparisons[5], label.y =5,method="wilcox.test",size = 3.5,label= "p.signif")
  print(p5)
}
dev.off()
##################################################################################
###########################Figure S1D Plot Starting###############################
##################################################################################
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
bar.df<-data.frame(proj@cellColData@listData[["celltype"]],proj@cellColData@listData[["sampletype"]],sapply(strsplit(proj@cellColData@rownames,'#'),function(x){x[[1]]}))
colnames(bar.df)<-c('celltype','sampletype','sample')
# save(bar.df,file='bar.df.RData')
load("/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/SPOT_test/bar.df.RData")
# The bar chart shows the proportion of cells
bar.df <- mutate(bar.df,name=factor(bar.df$sample, levels=c("ATAC_p1", "ATAC_p2", "ATAC_p3","ATAC_p4", "ATAC_n",  "ATAC1", "ATAC2", "ATAC3")))
# Assume df is your DataFrame and group is the column name you want to modify
bar.df$group <- factor(bar.df$sample, 
                       levels = c("ATAC_p1", "ATAC_p2", "ATAC_p3", "ATAC_p4", "ATAC_n", "ATAC1", "ATAC2", "ATAC3"),
                       labels = c("GSM5589391","GSM5589392","GSM5589393","GSM5589394",  "ATAC_n","ATAC1", "ATAC2", "ATAC3"))

text.df <- as.data.frame(table(bar.df$name)) 
color_cluster = c('B'="#25ABAE",'T'= "#E3837B",'NKT'="#D1669B",'ENDOCRINE'='#E06AA4','DUC_n'="#C386B8",'ACINAR'="#88C5EC",'DUC'="#6F4598",
                  'DC'="#999B24",'MAC'="#79B85B",'FIBROBLAST'="#79BEA5",'ENDOTHELIAL'="#EBA05F")

ggplot(bar.df,aes(x=group))+
  geom_bar(aes(fill=celltype),position = "fill",width = .7)+
  scale_x_discrete("")+
  scale_y_continuous("Total cell proportion",expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual("Cell types",values = color_cluster)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+ 
  scale_fill_manual("Cell types", 
                    values = color_cluster, 
                    labels = c('B' = "B Cell", 
                               'T' = "T Cell", 
                               'NKT' = "NKT Cell", 
                               'ENDOCRINE' = "Endocrine Cell", 
                               'DUC_n' = "Normal Ductal", 
                               'ACINAR' = "Acinar Cell", 
                               'DUC' = "Tumor Ductal", 
                               'DC' = "Dendritic Cell", 
                               'MAC' = "Macrophage", 
                               'FIBROBLAST' = "Fibroblast", 
                               'ENDOTHELIAL' = "Endothelial Cell")) +  
  coord_flip() 
ggplot(bar.df,aes(x=group))+
  geom_bar(aes(fill=celltype),position = "stack",width = .7)+
  scale_x_discrete("",position = "top")+ 
  scale_fill_manual(values = color_cluster)+
  scale_y_reverse("Total cell number",expand = c(0,0),position = "right",limits=c(5000,0))+ 
  theme(
    panel.grid = element_blank(),
    legend.position = "none", 
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+
  
  coord_flip()

