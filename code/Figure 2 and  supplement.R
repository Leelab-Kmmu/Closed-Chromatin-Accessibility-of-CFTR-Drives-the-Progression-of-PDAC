#R包引用
{
  library(devtools)
  library(Seurat)
  library(tidyverse)
  #library(SingleR)
  library(patchwork)
  library(celldex)
  library(rstan)
  library(monocle)
  #library(devtools)
  library(dyno)
  #library(monocle3)
  library(Matrix)
  library(dynmethods)
  library(dynwrap)
  library(hdf5r)
  #library(clusterProfiler)
  library(ggpubr) # 继承ggplot语法
  library(patchwork) # 拼图包
  library(ggsci)
  #remotes::install_github(repo = "samuel-marsh/scCustomize")
  library(scCustomize)
  library(ggplot2)
  library(RImagePalette)
  library(png)
  library(RImagePalette)
  library(imager)
  library(scales)
  library(ggplot2)
  library(ggprism)
  library(reshape)
  library(ggalluvial)
  library(tidydr)
  library(ggrepel)
}
###############################################################################################################################
###########################################Figure 2A###########################################################################
###############################################################################################################################
#figure -1--------umap图谱绘制
load("/400T/ckn/yxa/PDAC_sc_data.Rdata")
#save(scRNA,file="/public3/home/wangmeiheng/data_result/scRNA/yixian_scRNA.Rdata")

#color=c('#b35275','#3f9ac3')
#color2=c("#f08984","#3f9ac3","#28b8b8","#7f88b5","#b993aa","#89b456","#74c8b1","#d9a960",'#e3cf88')##c4bf82
color2=c("#ea888a","#cb73a9","#8cb255","#28b8b8","#82c1b8","#eac8b5","#75c8b2","#e6a3a1",'#e3cf88')
scRNA
#scRNA<-RunTSNE(scRNA,pc.num=1:18)

#DimPlot(scRNA, reduction = 'umap', group.by = "sample",
#        cols = color2,
#        pt.size = 1.5,
#        label = T,label.box = F
#) 

#提取umap二维数据
df=scRNA@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=scRNA@meta.data$clean_celltype)

head(df)
colnames(df)<-c("umap_1","umap_2","clean_celltype")

p_umap=ggplot(df, aes(umap_1, umap_2, color=clean_celltype,fill = clean_celltype))+
  geom_point(size = 0.9) + #细胞亚群点的大小
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        #panel.grid = element_blank(),
        # panel.grid.major = element_blank(), #网格线
        # panel.grid.minor = element_blank(), #网格线
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),#背景色
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #设置legend标签的大小
        legend.key.size=unit(0.6,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 4)))+ #图例点的大小
  scale_color_manual(values=color2)+theme_dr(xlength = 0.2, ylength = 0.2,
                                             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())

p_umap

###############################################################################################################################
###########################################Figure 2  A(右上)#########################################################################
###############################################################################################################################

#color=c('#b35275','#3f9ac3')
#color2=c("#f08984","#3f9ac3","#28b8b8","#7f88b5","#b993aa","#89b456","#74c8b1","#d9a960",'#e3cf88')##c4bf82
color2=c("#ea888a","#cb73a9","#8cb255","#28b8b8","#82c1b8","#eac8b5","#75c8b2","#e6a3a1",'#e3cf88')
scRNA
#scRNA<-RunTSNE(scRNA,pc.num=1:18)

#DimPlot(scRNA, reduction = 'umap', group.by = "sample",
#        cols = color2,
#        pt.size = 1.5,
#        label = T,label.box = F
#) 

#提取umap二维数据
df=scRNA@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=scRNA@meta.data$sample)

head(df)
colnames(df)<-c("umap_1","umap_2","sample")

p_umap=ggplot(df, aes(umap_1, umap_2, color=sample,fill = sample))+
  geom_point(size = 0.9) + #细胞亚群点的大小
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        #panel.grid = element_blank(),
        # panel.grid.major = element_blank(), #网格线
        # panel.grid.minor = element_blank(), #网格线
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),#背景色
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #设置legend标签的大小
        legend.key.size=unit(0.6,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 4)))+ #图例点的大小
  scale_color_manual(values=color2)+theme_dr(xlength = 0.2, ylength = 0.2,
                                             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())

p_umap



###############################################################################################################################
###########################################Figure 2S  B#########################################################################
###############################################################################################################################

#color=c('#b35275','#3f9ac3')
#color2=c("#f08984","#3f9ac3","#28b8b8","#7f88b5","#b993aa","#89b456","#74c8b1","#d9a960",'#e3cf88')##c4bf82
# 方法1：使用RColorBrewer生成渐变/分类颜色（推荐，颜色美观）
library(RColorBrewer)
# 获取亚群数量
n_clusters <- length(unique(df$seurat_clusters))
# 生成颜色向量（根据亚群数量选择配色方案）
if (n_clusters <= 12) {
  # 亚群数≤12：用Set3（颜色区分度高）
  color2 <- brewer.pal(n_clusters, "Set3")
} else {
  # 亚群数>12：用hcl颜色空间生成连续且区分度高的颜色
  color2 <- hcl(h = seq(0, 360, length.out = n_clusters + 1), c = 100, l = 65)[-1]
}

scRNA
#scRNA<-RunTSNE(scRNA,pc.num=1:18)

#DimPlot(scRNA, reduction = 'umap', group.by = "sample",
#        cols = color2,
#        pt.size = 1.5,
#        label = T,label.box = F
#) 

#提取umap二维数据
df=scRNA@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=scRNA@meta.data$seurat_clusters)

head(df)
colnames(df)<-c("umap_1","umap_2","seurat_clusters")

p_umap=ggplot(df, aes(umap_1, umap_2, color=seurat_clusters,fill = seurat_clusters))+
  geom_point(size = 0.9) + #细胞亚群点的大小
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        #panel.grid = element_blank(),
        # panel.grid.major = element_blank(), #网格线
        # panel.grid.minor = element_blank(), #网格线
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),#背景色
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #设置legend标签的大小
        legend.key.size=unit(0.6,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 4)))+ #图例点的大小
  scale_color_manual(values=color2)+theme_dr(xlength = 0.2, ylength = 0.2,
                                             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())

p_umap



###############################################################################################################################
###########################################Figure 2S  D#########################################################################
###############################################################################################################################


#color=c('#b35275','#3f9ac3')
#color2=c("#f08984","#3f9ac3","#28b8b8","#7f88b5","#b993aa","#89b456","#74c8b1","#d9a960",'#e3cf88')##c4bf82
color2=c("#ea888a","#cb73a9","#8cb255","#28b8b8","#82c1b8","#eac8b5","#75c8b2","#e6a3a1",'#e3cf88')
scRNA
#scRNA<-RunTSNE(scRNA,pc.num=1:18)

#DimPlot(scRNA, reduction = 'umap', group.by = "sample",
#        cols = color2,
#        pt.size = 1.5,
#        label = T,label.box = F
#) 

#提取umap二维数据
df=scRNA@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=scRNA@meta.data$tissue)

head(df)
colnames(df)<-c("umap_1","umap_2","tissue")

p_umap=ggplot(df, aes(umap_1, umap_2, color=tissue,fill = tissue))+
  geom_point(size = 0.9) + #细胞亚群点的大小
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        #panel.grid = element_blank(),
        # panel.grid.major = element_blank(), #网格线
        # panel.grid.minor = element_blank(), #网格线
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),#背景色
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #设置legend标签的大小
        legend.key.size=unit(0.6,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 4)))+ #图例点的大小
  scale_color_manual(values=color2)+theme_dr(xlength = 0.2, ylength = 0.2,
                                             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())

p_umap




###############################################################################################################################
###########################################Figure 2  C#########################################################################
###############################################################################################################################


FeaturePlot(scRNA, "CFTR", cols = c('grey',"#E54924"))+theme_dr(xlength = 0.2, ylength = 0.2,
                                                                arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
)+theme(panel.grid = element_blank())




###############################################################################################################################
###########################################Figure 2  B#########################################################################
###############################################################################################################################

# 定义marker基因
markerGenes  <- c('CD3D',#T_cell
                  'CD19','MS4A1','CD79A','CD79B','IGLL5',#B_cell
                  'MKI67',#恶性---Malignant_cell
                  'CDH5',#内皮---Endothelial_cell
                  'DCN','LUM',#成纤维---fibroblast_cell
                  'CHGB',#内分泌---Endocrin_cell
                  'PRSS1',#胰腺腺泡细胞---pancreatic_cell
                  'EPCAM',#上皮---Epithelial_cell
                  'ITGAX','KRT19','AMBP','TFF2',#DUC_cell --monocyte-derived
                  #'CD1C',#树突状细胞
                  'FCGR3A',#NK_cell
                  'AIF1','CD68')#巨噬---macrophage_cell

# 1. 提取细胞类型并自定义排列顺序
# 查看当前所有细胞类型（获取名称列表）
all_celltypes <- unique(scRNA@meta.data$clean_celltype)
print("当前细胞类型列表：")
print(all_celltypes)

# 2. 手动定义细胞类型的排列顺序（根据你的需求修改以下向量）
# 示例：按免疫细胞→基质细胞→恶性细胞的顺序排列
#custom_order <- c(
#  "T_cell", "B_cell", "Malignant_cell", # 免疫细胞
#  "Endothelial_cell","fibroblast_cell",   # 基质细胞
#  "Endocrin_cell","Epithelial_cell", "DUC_cell" , "pancreatic_cell","macrophage_cell" # 上皮/腺体相关
#    # 恶性细胞
#)

custom_order <- c(
  "T", "B", "Malignant", # 免疫细胞
  "Endothelial","fibroblast",   # 基质细胞
  "Endocrin","Epithelial", "DUC" , "pancreatic","macrophage" # 上皮/腺体相关
  # 恶性细胞
)

# 3. 确保自定义顺序包含所有细胞类型（避免遗漏）
# 检查是否有遗漏的细胞类型并补充到末尾
missing_types <- setdiff(all_celltypes, custom_order)
if (length(missing_types) > 0) {
  custom_order <- c(custom_order, missing_types)
  warning(paste("补充了未定义的细胞类型：", paste(missing_types, collapse = ", ")))
}

# 4. 重新设置细胞类型的因子水平（关键步骤）
scRNA@meta.data$clean_celltype <- factor(
  scRNA@meta.data$clean_celltype,
  levels = custom_order  # 按自定义顺序排列
)

# 5. 绘制DotPlot（此时y轴将按custom_order排列）
p1 <- DotPlot(
  scRNA,
  features = markerGenes,
  cols = c('grey', "#E54924"),
  group.by = "clean_celltype"  # 使用重新排序后的celltype
) +
  # 自定义颜色渐变，设置多个红色系节点实现阶梯下降
  scale_color_gradientn(
    colors = c("grey", "#FFCCCC", "#FF9999", "#FF6666", "#FF3333", "#E54924"),
    values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1))
  ) +
  theme(
    axis.text.x  = element_text(color = "black", size = 12, angle = 45, vjust = 0.5, hjust = 0.5),
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(2, "mm"),
    axis.line = element_blank()
  ) +
  labs(x = "", y = "")

print(p1)



###############################################################################################################################
###########################################Figure 2S  C#########################################################################
###############################################################################################################################

#dimplot图绘制
markerGenes  <- c('CD3D',#T_cell
                  'FCGR3A',#NK_cell
                  'CD19','MS4A1','CD79A','CD79B','IGLL5',#B_cell
                  'MKI67',#恶性---Malignant_cell
                  'CDH5',#内皮---Endothelial_cell
                  'DCN','LUM',#成纤维---fibroblast_cell
                  'CHGB',#内分泌---Endocrin_cell
                  'PRSS1',#胰腺腺泡细胞---pancreatic_cell
                  'EPCAM',#上皮---Epithelial_cell
                  'ITGAX','KRT19','AMBP','TFF2',#DUC_cell --monocyte-derived
                  'CD1C',#树突状细胞
                  'AIF1','CD68')#巨噬---macrophage_cell

p2=DotPlot(scRNA,
           features=markerGenes,cols = c('grey',"#E54924"),group.by = "seurat_clusters")+
  theme(
    axis.text.x  = element_text(color="black",size=12,angle = 45,vjust = 0.5, hjust=0.5),
    panel.border = element_rect(color="black"), #面板边框
    panel.spacing = unit(2, "mm"), #面板间距
    axis.line = element_blank(),
  )+labs(x="", y="") + coord_flip(); p2


###############################################################################################################################
###########################################Figure 2  E#########################################################################
###############################################################################################################################


df=scRNA@meta.data

my_comparisons <- list(c("DUC", "Endothelial"), 
                       c("DUC", "macrophage"),
                       c("DUC", "fibroblast"),
                       c("DUC", "T"),
                       c("DUC", "Endocrin"),
                       c("DUC", "B"),
                       c("DUC", "Epithelial"),
                       c("DUC", "Malignant"))



p <- ggplot(df, aes(x=clean_celltype, y=CFTR,color=clean_celltype)) +
  # 绘制箱线图
  geom_boxplot(aes(fill=clean_celltype),
               alpha=0.1)+ # 设置透明度
  # 绘制散点
  #geom_jitter()+
  # 设置颜色
  scale_color_manual(values = color2)+
  scale_fill_manual(values = color2)+
  # 设置主题
  theme_bw()+
  # 去除网格线
  theme(panel.grid = element_blank())+ylim(0.001,4.3)+
  stat_compare_means(method = "kruskal.test",
                     label = "p.format",
                     label.x = 0.8, 
                     label.y = 4.3,
                     show.legend = F) +
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(2.8,2.95,3.15,3.3,3.45, 3.6,3.75,3.9,4.05,4.2),
                     method="wilcox.test",
                     label="p.signif"
  )

p


###############################################################################################################################
###########################################Figure 2  D#########################################################################
###############################################################################################################################


library(ggpubr)
#mark='CFTR'
table(scRNA$tissue)
mark=c('CFTR','ZNF616','NFIA','MDFI','ZNF667','MYO10','ZFPL1','ENSA','PIFO','GGA1')
martix=scRNA@assays$RNA@data[intersect(row.names(scRNA),mark),]
#martix=as.data.frame(matrix)
color=c('#b35275','#3f9ac3',"pink")#c('#b35275','#3f9ac3',"pink")
#c("blue","red",'green')
plot.info_all <- martix %>%t()%>%as.data.frame()

plot.info_all=log2(plot.info_all+1)

plot.info_all$tissue=scRNA$tissue
plot.info_all=na.omit(plot.info_all)
plot.info_all=plot.info_all[,c('CFTR','tissue')] 

#colnames(plot.info_all)='CFTR'
#plot.info_all$tissue[which(is.na(plot.info_all),arr.ind = T)]='unknown'
library(reshape2)
df_long <- melt(plot.info_all, id.vars = "tissue", measure.vars = c(intersect(row.names(scRNA),'CFTR')))
#df_long
df_long$combine_name<-factor(paste(df_long$variable,df_long$tissue,sep="_"),levels = c('CFTR_normal','CFTR_tumor','CFTR_Metastatic'))
df_long$tissue<-factor(df_long$tissue,levels = c('normal','tumor','Metastatic'))
my_comparisons <- list(c("CFTR_normal", "CFTR_tumor"), 
                       c("CFTR_normal", "CFTR_Metastatic"),
                       c("CFTR_tumor", "CFTR_Metastatic"))
p=ggboxplot(df_long, x="combine_name", y="value", color = "tissue", 
            ylab="Gene expression",
            xlab="",
            legend.title='type',
            palette =color,
            width=0.4, add = "",
            bxp.errorbar=T,#显示误差条
            bxp.errorbar.width=0.1, #误差条大小
            size=0.5, #箱型图边线的粗细
            outlier.shape=NA)+#rotate_x_text(20)+#jitter
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(2.8,3.0,2.8),
                     method="wilcox.test",
                     label="p.signif")+
  ylim(0.001,4) # Add horizontal line at base mean
p




###############################################################################################################################
###########################################Figure 2  F#########################################################################
###############################################################################################################################

custom_palette <- c('pink','#b35275',"#3f9ac3")

#custom_palette <- c('#b35275',"#3f9ac3")
df=scRNA@meta.data

p <- ggplot(
  subset(df, CFTR != 0), 
  #df,# 过滤 CFTR = 0 的行
  aes(x = clean_celltype, y = CFTR, fill = tissue)
) +
  geom_boxplot(color = "black") +  # 箱体添加黑色边框
  labs(
    x = "Cell Type",
    y = "CFTR Expression",
    fill = "Tissue"
  ) +
  # 添加组内比较的显著性（tissue间比较）
  stat_compare_means(aes(group = tissue),  # 按tissue分组比较
                     method = "kruskal.test",  
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.15, 1), symbols = c("***", "**", "*", "ns")),
                     label = "p.signif",  # 显示显著性符号（*、**等）
                     hide.ns = T) +  # 隐藏不显著（ns）的标注
  scale_fill_manual(values = custom_palette) +  # 应用自定义颜色
  theme_minimal() +  # 简约主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",  # 图例放在顶部
    # 移除网格线但保留坐标轴
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 确保坐标轴可见（theme_minimal默认会淡化轴线，这里加强显示）
    axis.line = element_line(color = "black"),  # 添加强制的轴线
    axis.ticks = element_line(color = "black")  # 添加强制的刻度线
  )

print(p)







###############################################################################################################################
###########################################Figure 2S  E#########################################################################
###############################################################################################################################


df=scRNA@meta.data
library(ggplot2)  # 绘图核心包
library(dplyr)    # 数据预处理

# 统计每个 celltype 下不同 tissue 的细胞数量
count_df <- df %>%
  group_by(clean_celltype, tissue) %>%
  summarise(count = n()) %>%
  ungroup()

# 提取 celltype 下划线前的前缀（如 DUC_cell → DUC）
count_df$celltype_short <- str_split_fixed(count_df$clean_celltype, "_", 2)[, 1]

ggplot(count_df, aes(x = celltype_short, y = count, fill = tissue)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Cell Count by Celltype and Tissue",
    x = "Cell Type",
    y = "Cell Count",
    fill = "Tissue"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

count_df <- df %>%
  group_by(clean_celltype, sample) %>%
  summarise(count = n()) %>%
  ungroup()

# 提取 celltype 下划线前的前缀（如 DUC_cell → DUC）
count_df$celltype_short <- str_split_fixed(count_df$clean_celltype, "_", 2)[, 1]



ggplot(count_df, aes(x = celltype_short, y = count, fill = sample)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Cell Count by Celltype and Sample",
    x = "Celltype",
    y = "Cell Count",
    fill = "Sample"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
























