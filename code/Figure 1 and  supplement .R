############Figure 1
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
##########################2violin
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

###降维与聚类
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
###########################Figure 1B (右上) Plot Starting#########################
##################################################################################
ggAlignPlots(p1, type = "h")#Sample Umap

##################################################################################
###########################Figure S1B (右上) Plot Starting########################
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
    legend.text = element_text(size = 12),  # 调整图例标签文本的大小
    legend.title = element_text(size = 14)  # 调整图例标题的大小
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 6))  # 调整图例中圆圈的大小
  )
ggAlignPlots(p4, type = "h")
# plotPDF(p4, name = "PDAC-UMAP-sampletype.pdf",
#         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)


#cell type
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
###########################Figure 1B Plot Starting################################
##################################################################################
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "celltype", embedding = "UMAP",
                    pal=c('B'="#25ABAE",'T'= "#E3837B",'NKT'="#D1669B",'ENDOCRINE'='#E06AA4','DUC_n'="#C386B8",'ACINAR'="#88C5EC",'DUC'="#6F4598",
                          'DC'="#999B24",'MAC'="#79B85B",'FIBROBLAST'="#79BEA5",'ENDOTHELIAL'="#EBA05F")) +
  theme(
    legend.text = element_text(size = 12),  # 调整图例标签文本的大小
    legend.title = element_text(size = 14)  # 调整图例标题的大小
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 6))  # 调整图例中圆圈的大小
  )
p3
# plotPDF(p3, name = "Plot-UMAP-cell.pdf",
#         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
# proj_umap <- saveArchRProject(ArchRProj = proj)

#######识别marker gene
markerGenes  <- c('CD3D',#T
                  'FCGR3A',#NK
                  'CD19','MS4A1','CD79A','CD79B','IGLL5',#B
                  'KRT19','MKI67',#恶性
                  'CDH5',#内皮
                  'DCN','LUM',#成纤维
                  'CHGB',#内分泌
                  'PRSS1',#胰腺腺泡细胞
                  'EPCAM',#上皮
                  'AIF1','ITGAX','CD1C',#DC
                  'CD68')#巨噬
##################################################################################
###########################Figure 1C Plot Starting################################
##################################################################################
#########################################marker gene点图绘制
detach("package:plyr", unload=TRUE)
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
colorMat <- ArchR:::.getMatrixValues(ArchRProj = proj, name = markerGenes,matrixName = "GeneScoreMatrix", log2Norm = FALSE)
temp_mat<-data.frame(t(colorMat))
#1整理dotplot绘图数据
cell_annotations <- data.frame(
  Cell = rownames(temp_mat),
  CellType = proj@cellColData@listData[["celltype"]],
  Cluster = proj@cellColData@listData[["Clusters"]]
)
long_data <- temp_mat %>%
  rownames_to_column(var = "Cell") %>%
  pivot_longer(cols = -Cell, names_to = "Gene", values_to = "Expression") %>%
  left_join(cell_annotations, by = "Cell")

# 2计算平均标准化表达值和表达比例
gene_data_avg <- long_data %>%
  group_by(Gene, CellType) %>%
  summarise(
    AverageStdExpression = mean(Expression, na.rm = TRUE),  # 使用标准化后的表达值
    CellCount = dplyr::n(),  # 计算细胞数量
    
    ExpressionProportion = sum(Expression > 0) / length(Expression) * 100,  # 计算表达的细胞比例
    .groups = 'drop'
  )
# 3标准化表达值z-score）
gene_data_avg <- gene_data_avg %>%
  group_by(Gene) %>%
  mutate(StdExpression = (AverageStdExpression - mean(AverageStdExpression, na.rm = TRUE)) / sd(AverageStdExpression, na.rm = TRUE)) %>%
  ungroup()

# 将 Gene 设置为因子，并按照 markerGenes 的顺序排列
gene_data_avg$Gene <- factor(gene_data_avg$Gene, levels = markerGenes)
gene_data_avg$CellType<-factor(gene_data_avg$CellType,levels =c('T','NKT','B','DUC','ENDOTHELIAL','FIBROBLAST','ENDOCRINE','ACINAR','DUC_n','DC','MAC'))
p1 <- ggplot(gene_data_avg, aes(x = Gene, y = CellType)) +
  geom_point(aes(size = ExpressionProportion, color = StdExpression)) +  # 点大小表示表达比例，颜色表示表达值
  scale_color_gradient2(low = "steelblue", mid = "lightgrey", high = "darkgoldenrod1", midpoint = 1) +  # 设置颜色渐变
  theme_minimal() +  # 使用简约主题
  labs(size = "ExpressionProportion", color = "Gene accessibility score") +  # 添加图例标签
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # x轴标签旋转
    panel.grid = element_blank(),  # 移除网格线
    panel.border = element_blank(),  # 移除边框
    axis.line = element_line()  # 添加坐标轴线
  )
p1
##################################################################################
###########################Figure S1C Plot Starting###############################
##################################################################################
gene_data_avg <- long_data %>%
  group_by(Gene, Cluster) %>%
  summarise(
    AverageStdExpression = mean(Expression, na.rm = TRUE),  # 使用标准化后的表达值
    CellCount = dplyr::n(),  # 计算细胞数量
    
    ExpressionProportion = sum(Expression > 0) / length(Expression) * 100,  # 计算表达的细胞比例
    .groups = 'drop'
  )
#标准化表达值z-score 
gene_data_avg <- gene_data_avg %>%
  group_by(Gene) %>%
  mutate(StdExpression = (AverageStdExpression - mean(AverageStdExpression, na.rm = TRUE)) / sd(AverageStdExpression, na.rm = TRUE)) %>%
  ungroup()

#绘图
# 将 Gene 设置为因子，并按照 markerGenes 的顺序排列
gene_data_avg$Gene <- factor(gene_data_avg$Gene, levels = markerGenes)
gene_data_avg$Cluster<-factor(gene_data_avg$Cluster,levels = c('C4','C5','C7','C6','C1','C2','C3','C12','C13','C18','C17','C19',
                                                               'C8','C10','C11','C9','C14','C15','C16'))
p2 <- ggplot(gene_data_avg, aes(x = Gene, y = Cluster)) +
  geom_point(aes(size = ExpressionProportion, color = StdExpression)) +  # 点大小表示表达比例，颜色表示表达值
  scale_color_gradient2(low = "steelblue", mid = "lightgrey", high = "darkgoldenrod1", midpoint = 1) +  # 设置颜色渐变
  theme_minimal() +  # 使用简约主题
  labs(size = "ExpressionProportion", color = "Gene accessibility score") +  # 添加图例标签
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # x轴标签旋转
    panel.grid = element_blank(),  # 移除网格线
    panel.border = element_blank(),  # 移除边框
    axis.line = element_line()  # 添加坐标轴线
  )
p2
##################################################################################
###########################Figure 1D Plot Starting################################
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
###########################Figure 1E Plot Starting################################
##################################################################################
proj <- loadArchRProject(path = "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments/singlecellATAC")
proj <- addImputeWeights(proj)
# 加载包
library(presto)
library("hexbin")
colorMat <- ArchR:::.getMatrixValues(ArchRProj = proj, name = markerGenes,matrixName = "GeneScoreMatrix", log2Norm = FALSE)
colorMat <- imputeMatrix(mat = as.matrix(colorMat),imputeWeights = getImputeWeights(proj))
#save(colorMat,file='/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/SPOT_test/colorMat.Rdata')

load('/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/SPOT_test/colorMat.Rdata')

################################CFTR gene score在不同细胞类型boxplot
marker_TSS<-colorMat
CFTR_tss_150_data<-marker_TSS['CFTR',]

df<-as.data.frame(CFTR_tss_150_data)
colnames(df)<-c("CFTR")
df$type<-proj@cellColData@listData[["celltype"]]
df_melt<-reshape::melt(df,id.vars=c("type"))
df_melt$combine_name<-factor(paste(df_melt$variable,df_melt$type,sep="_"))
#levels设置横坐标顺序
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
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),  # 调整 x 轴标签的倾斜角度
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 12),  # 调整 y 轴标签
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次网格线
    axis.title.x = element_blank(),  # 移除 x 轴标题
    axis.title.y = element_text(size = rel(1)),  # 调整 y 轴标题字体大小
    panel.border = element_rect(color = "black", size = 1)  # 设置完整的黑色边框
  ) +
  scale_fill_manual(values = color) +
  labs(y = "Gene accessibility score") +  # 更改 y 轴标题
  scale_x_discrete(labels = c(
    "CFTR_DUC_n" = "Normal ductal", "CFTR_DUC" = "Tumor ductal", "CFTR_FIBROBLAST" = "Fibroblast", 
    "CFTR_T" = "T", "CFTR_NKT" = "NKT", "CFTR_MAC" = "Macrophage",
    "CFTR_ENDOTHELIAL" = "Endothelial", "CFTR_ENDOCRINE" = "Endocrine",
    "CFTR_DC" = "DC", "CFTR_B" = "B", "CFTR_ACINAR" = "Acinar"))  # 修改 x 轴标签
p+stat_compare_means(comparisons = comparisons,label= "p.signif")
##################################################################################
###########################Figure 1E Plot Starting################################
##################################################################################
######################################箱式图--healthy,PDAC,STAD,other cancers, Noncancer,and gastritis
###########导入数据
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
###################################plot
pdf("CFTR_5_cfDNA_TSScoverage.pdf",width = 8, height = 8)
for (i in 2) {
  CFTR_tss_150_data<-tss_150_file[which(tss_pro_table$gene_name=="CFTR")[i],]
  df<-as.data.frame(t(CFTR_tss_150_data))
  colnames(df)<-c("CFTR")
  df$type<-c(rep("PDAC",39),rep("Normal",100),rep("gastric",101),rep("TGY",38),rep("NGY",58),rep("GGY",17))
  df_melt<-reshape::melt(df,id.vars=c("type"))
  df_melt$combine_name<-factor(paste(df_melt$variable,df_melt$type,sep="_"))
  # 设置横坐标顺序
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
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),  # 调整 x 轴标签的倾斜角度
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 12),  # 调整 y 轴标签
          panel.grid.major = element_blank(),  # 移除主网格线
          panel.grid.minor = element_blank(),  # 移除次网格线
          axis.title.x = element_blank(),  # 移除 x 轴标题
          axis.title.y = element_text(size = rel(1)),  # 调整 y 轴标题字体大小
          panel.border = element_rect(color = "black", size = 1)  # 设置完整的黑色边框
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
#save(bar.df,file='bar.df.RData')
load("/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/SPOT_test/bar.df.RData")
#条形图展示细胞占比
bar.df <- mutate(bar.df,name=factor(bar.df$sample, levels=c("ATAC_p1", "ATAC_p2", "ATAC_p3","ATAC_p4", "ATAC_n",  "ATAC1", "ATAC2", "ATAC3")))
# 假设 df 是你的数据框，并且 group 是你想修改的列名
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
  )+ # 更改x轴和y轴的标题名称
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
                               'ENDOTHELIAL' = "Endothelial Cell")) +  # 修改填充标签
  coord_flip() #让条形图横过来


ggplot(bar.df,aes(x=group))+
  geom_bar(aes(fill=celltype),position = "stack",width = .7)+
  scale_x_discrete("",position = "top")+  #x轴在上
  scale_fill_manual(values = color_cluster)+
  scale_y_reverse("Total cell number",expand = c(0,0),position = "right",limits=c(5000,0))+  #y轴在右
  theme(
    panel.grid = element_blank(),
    legend.position = "none", # 不加图例
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+
  
  coord_flip()


###############################################################################################################################
###########################################Figure 1SS  A#########################################################################
###############################################################################################################################
library(maftools)
library(maftools)
library(mclust)
library(NMF)
library(pheatmap)
#library(barplot3d)
setwd("/400T/wangmeiheng/TCGA/PDAC_mutation")
#注意--数据读取！
##all_mut的数据读取的原数据没放入服务器
files <- list.files(pattern = '*.gz',recursive = TRUE)

files=files[seq(1,length(files),2)]
#取奇数列----原因目前只能通过偶数列均无法读取来解答
all_mut <- data.frame()
for (file in files) {
  mut <- read.delim(file,skip = 7,, header = T, fill = TRUE,sep = "\t")#skip = 7,
  all_mut <- rbind(all_mut,mut)
}

all_mut <- read.maf(all_mut)

a <- all_mut@data %>%
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>%
  as.data.frame() %>%
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

gene <- as.character(unique(a$Hugo_Symbol))
sample <- as.character(unique(a$Tumor_Sample_Barcode))

mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))

for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}

for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

#所有样本突变情况汇总/排序
gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))
gene_top <- gene_count$gene[1:20] 



library(RColorBrewer)
# 生成8个暖色调（使用YlOrRd系列，偏黄橙红）
vc_cols <- brewer.pal(n = 8, name = "YlOrRd")  # 暖色调：黄→橙→红渐变

# 或使用Oranges系列（纯橙色系，更柔和）
# vc_cols <- brewer.pal(n = 8, name = "Oranges")

# 命名与突变类型对应（保持原顺序）
# 手动定义高对比度颜色（命名需与突变类型完全匹配）
vc_cols <- c(
  "Missense_Mutation" = "#FF2000",    # 橙色（错义突变）
  "Nonsense_Mutation" = "#90FF00",    # 红色（无义突变）
  "Frame_Shift_Del" = "#800080",      # 紫色（移码缺失）
  "Frame_Shift_Ins" = "#0000FF",      # 深蓝色（移码插入）
  "In_Frame_Ins" = "#FFA000",         # 绿色（框内插入）
  "In_Frame_Del" = "#FFFF00",         # 黄色（框内缺失）
  "Splice_Site" = "#FF00FF",          # 品红（剪接位点突变）
  "Multi_Hit" = "#A52A2A"             # 棕色（多 hit 突变）
)
# 查看颜色
print(vc_cols)
# 输出示例（YlOrRd）：
# Frame_Shift_Del   Missense_Mutation  Nonsense_Mutation        Multi_Hit 
#      "#FFFFCC"           "#FFEDA0"           "#FED976"           "#FEB24C" 
# Frame_Shift_Ins    In_Frame_Ins     Splice_Site      In_Frame_Del 
#      "#FD8D3C"           "#FC4E2A"           "#E31A1C"           "#B10026" 
par(family = "Arial")  # 设置字体为 Arial（需系统安装该字体）
oncoplot(
  maf = all_mut, 
  colors = vc_cols,
  top = 30,
  fontSize = 0.6,  # 基础字体大小（会影响大部分文本）
  showTumorSampleBarcodes = FALSE
)










###############################################################################################################################
###########################################Figure 1SS  C#########################################################################
###############################################################################################################################

# 绘制oncostrip
oncostrip(
  maf = all_mut,
  genes = c('KRAS', "TP53", "CFTR"),
  colors = vc_cols
  # 如需其他参数可在此添加（如sort = TRUE等）
)



###############################################################################################################################
###########################################Figure 1SS  B#########################################################################
###############################################################################################################################

plotmafSummary(maf = all_mut, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,color = vc_cols)






###############################################################################################################################
###########################################Figure 1SS  D/E#########################################################################
###############################################################################################################################


library(BSgenome)
#head(available.genomes())

#BiocManager::install("MutationalPatterns", update = FALSE)
###hg19数据暂时无法下载
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome.Hsapiens.UCSC.hg19)
library(MutationalPatterns)
# 3. 读取VCF文件（关键：定义ref_genome为加载的基因组对象）
# 2. 读取VCF文件（修改为你的路径）
vcf_files <- list.files(
  path = "/400T/ckn/yxa/vcf",  
  pattern = ".vcf",  # 或.vcf.gz
  full.names = TRUE,
  recursive = TRUE  # 递归查找子文件夹
)
sample_name=basename(vcf_files)
vcfs <- read_vcfs_as_granges(vcf_files, sample_name, ref_genome)


type_occurrences <- mut_type_occurrences(vcfs, ref_genome)

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(mut_mat, condensed = TRUE)
#plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE)
#plot_compare_profiles(type_occurrences,profile_names = tissue)

####整合18个样本突变位点信息
mut_mat=as.data.frame(mut_mat)
mut_mat$PDAC=apply(mut_mat,1,sum)
mut_mat1=as.data.frame(mut_mat[,19])%>%as.matrix()
colnames(mut_mat1)='PDAC'
row.names(mut_mat1)=row.names(mut_mat)
plot_96_profile(mut_mat1, condensed = TRUE)
#比较两个样本间突变频谱分布的差异
#plot_compare_profiles(mut_mat[,1], mut_mat[,2], condensed = TRUE)

