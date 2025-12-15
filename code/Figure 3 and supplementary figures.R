#R package
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
  library(ggpubr)
  library(patchwork)
  library(ggsci)
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
###########################################Figure 3A###########################################################################
###############################################################################################################################
#load("/400T/ckn/yxa/scRNA_PDAC.Rdata")
#save(scRNA,file="/400T/wangmeiheng/scRNA/scRNA.Rdata")
load("/400T/wangmeiheng/scRNA/scRNA.Rdata")
#color=c('#b35275','#3f9ac3')
#color2=c("#f08984","#3f9ac3","#28b8b8","#7f88b5","#b993aa","#89b456","#74c8b1","#d9a960",'#e3cf88')##c4bf82
color2=c("#ea888a","#cb73a9","#8cb255","#28b8b8","#82c1b8","#eac8b5","#75c8b2","#e6a3a1",'#e3cf88')
scRNA
#Extract UMAP two-dimensional data
df=scRNA@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=scRNA@meta.data$clean_celltype)

head(df)
colnames(df)<-c("umap_1","umap_2","clean_celltype")

p_umap=ggplot(df, aes(umap_1, umap_2, color=clean_celltype,fill = clean_celltype))+
  geom_point(size = 0.9) + 
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(),
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'),
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), 
        legend.key.size=unit(0.6,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 4)))+ 
  scale_color_manual(values=color2)+theme_dr(xlength = 0.2, ylength = 0.2,
                                             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())

p_umap

###############################################################################################################################
###########################################Figure 3  A(Top right)##############################################################
###############################################################################################################################

color2=c("#ea888a","#cb73a9","#8cb255","#28b8b8","#82c1b8","#eac8b5","#75c8b2","#e6a3a1",'#e3cf88')
scRNA
#Extract UMAP two-dimensional data
df=scRNA@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=scRNA@meta.data$sample)

head(df)
colnames(df)<-c("umap_1","umap_2","sample")

p_umap=ggplot(df, aes(umap_1, umap_2, color=sample,fill = sample))+
  geom_point(size = 0.9) + 
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), 
        legend.key.size=unit(0.6,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 4)))+ 
  scale_color_manual(values=color2)+theme_dr(xlength = 0.2, ylength = 0.2,
                                             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())

p_umap

###############################################################################################################################
###########################################Figure S2 A#########################################################################
###############################################################################################################################
selected_samples <- c("CRA001160","GSE111672","GSE141017","GSE154778","GSE162708","GSE165399","GSE148673","GSE158356")
selected_samples1 <- c("GSE154763")
VlnPlot(scRNA, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample", 
        pt.size = 0,       
        ncol = 3) +         
  NoLegend()
VlnPlot(subset(scRNA, subset = sample %in% selected_samples), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample", 
        pt.size = 0,          
        ncol = 3) +         
  NoLegend()
VlnPlot(subset(scRNA, subset = sample %in% selected_samples1), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample",  
        pt.size = 0,         
        ncol = 3) +          
  NoLegend()
table(scRNA$sample)

###############################################################################################################################
###########################################Figure S2  B########################################################################
###############################################################################################################################

library(RColorBrewer)
n_clusters <- length(unique(df$seurat_clusters))
if (n_clusters <= 12) {
  color2 <- brewer.pal(n_clusters, "Set3")
} else {
  color2 <- hcl(h = seq(0, 360, length.out = n_clusters + 1), c = 100, l = 65)[-1]
}

scRNA
#Extract UMAP 2D data
df=scRNA@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=scRNA@meta.data$seurat_clusters)

head(df)
colnames(df)<-c("umap_1","umap_2","seurat_clusters")

p_umap=ggplot(df, aes(umap_1, umap_2, color=seurat_clusters,fill = seurat_clusters))+
  geom_point(size = 0.9) + 
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), 
        legend.key.size=unit(0.6,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 4)))+
  scale_color_manual(values=color2)+theme_dr(xlength = 0.2, ylength = 0.2,
                                             arrow = grid::arrow(length = unit(0.15, "inches"), type = "closed")
  )+theme(panel.grid = element_blank())

p_umap

###############################################################################################################################
###########################################Figure S2  D########################################################################
###############################################################################################################################

color2=c("#ea888a","#cb73a9","#8cb255","#28b8b8","#82c1b8","#eac8b5","#75c8b2","#e6a3a1",'#e3cf88')
scRNA
#Extract UMAP two-dimensional data
df=scRNA@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=scRNA@meta.data$tissue)

head(df)
colnames(df)<-c("umap_1","umap_2","tissue")

p_umap=ggplot(df, aes(umap_1, umap_2, color=tissue,fill = tissue))+
  geom_point(size = 0.9) + 
  scale_color_manual(values=color2)+
  theme(panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'),
        plot.background=element_rect(fill="white"),
        legend.title = element_blank(), 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14),
        legend.key.size=unit(0.6,'cm') )+
  guides(fill= guide_legend(override.aes = list(size = 4)))+ 
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
###########################################Figure 3  B#########################################################################
###############################################################################################################################

#marker gene
markerGenes  <- c('CD3D',#T
                  'CD19','MS4A1','CD79A','CD79B','IGLL5',#B
                  'MKI67',#Malignant
                  'CDH5',#Endothelial
                  'DCN','LUM',#Fibroblast
                  'CHGB',#Endocrine
                  'PRSS1',#Pancreatic acinar cells
                  'EPCAM',#Epithelial
                  'ITGAX','KRT19','AMBP','TFF2',#DUC
                  'CD1C',#DCs
                  'FCGR3A',#NK_
                  'AIF1','CD68')#Macrophage

all_celltypes <- unique(scRNA@meta.data$clean_celltype)
print(all_celltypes)

custom_order <- c(
  "T", "B", "Malignant", 
  "Endothelial","fibroblast",  
  "Endocrin","Epithelial", "DUC" , "pancreatic","macrophage" 
)

missing_types <- setdiff(all_celltypes, custom_order)
scRNA@meta.data$clean_celltype <- factor(
  scRNA@meta.data$clean_celltype,
  levels = custom_order  
)

#ploting
p1 <- DotPlot(
  scRNA,
  features = markerGenes,
  cols = c('grey', "#E54924"),
  group.by = "clean_celltype" 
) +
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
###########################################Figure S2  C########################################################################
###############################################################################################################################

#dimplot
markerGenes  <- c('CD3D',
                  'FCGR3A',
                  'CD19','MS4A1','CD79A','CD79B','IGLL5',
                  'MKI67',
                  'CDH5',
                  'DCN','LUM',
                  'CHGB',
                  'PRSS1',
                  'EPCAM',
                  'ITGAX','KRT19','AMBP','TFF2',
                  'CD1C',
                  'AIF1','CD68')

p2=DotPlot(scRNA,
           features=markerGenes,cols = c('grey',"#E54924"),group.by = "seurat_clusters")+
  theme(
    axis.text.x  = element_text(color="black",size=12,angle = 45,vjust = 0.5, hjust=0.5),
    panel.border = element_rect(color="black"), 
    panel.spacing = unit(2, "mm"), 
    axis.line = element_blank(),
  )+labs(x="", y="") + coord_flip(); p2


###############################################################################################################################
###########################################Figure 3  E#########################################################################
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
  geom_boxplot(aes(fill=clean_celltype),
               alpha=0.1)+
  scale_color_manual(values = color2)+
  scale_fill_manual(values = color2)+
  theme_bw()+
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
###########################################Figure 3  D#########################################################################
###############################################################################################################################


library(ggpubr)
table(scRNA$tissue)
mark=c('CFTR','ZNF616','NFIA','MDFI','ZNF667','MYO10','ZFPL1','ENSA','PIFO','GGA1')
martix=scRNA@assays$RNA@data[intersect(row.names(scRNA),mark),]
color=c('#b35275','#3f9ac3',"pink")#c('#b35275','#3f9ac3',"pink")
plot.info_all <- martix %>%t()%>%as.data.frame()

plot.info_all=log2(plot.info_all+1)

plot.info_all$tissue=scRNA$tissue
plot.info_all=na.omit(plot.info_all)
plot.info_all=plot.info_all[,c('CFTR','tissue')] 
library(reshape2)
df_long <- melt(plot.info_all, id.vars = "tissue", measure.vars = c(intersect(row.names(scRNA),'CFTR')))
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
            bxp.errorbar=T,
            bxp.errorbar.width=0.1, 
            size=0.5, 
            outlier.shape=NA)+
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(2.8,3.0,2.8),
                     method="wilcox.test",
                     label="p.signif")+
  ylim(0.001,4) 
p


###############################################################################################################################
###########################################Figure 3  F#########################################################################
###############################################################################################################################

custom_palette <- c('pink','#b35275',"#3f9ac3")

#custom_palette <- c('#b35275',"#3f9ac3")
df=scRNA@meta.data

p <- ggplot(
  subset(df, CFTR != 0), 
  aes(x = clean_celltype, y = CFTR, fill = tissue)
) +
  geom_boxplot(color = "black") +  
  labs(
    x = "Cell Type",
    y = "CFTR Expression",
    fill = "Tissue"
  ) +
  stat_compare_means(aes(group = tissue),  
                     method = "kruskal.test",  
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.15, 1), symbols = c("***", "**", "*", "ns")),
                     label = "p.signif",
                     hide.ns = T) +  
  scale_fill_manual(values = custom_palette) +  
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black")  
  )

print(p)
###############################################################################################################################
###########################################Figure S2  E########################################################################
###############################################################################################################################


df=scRNA@meta.data
library(ggplot2)  
library(dplyr)    
count_df <- df %>%
  group_by(clean_celltype, tissue) %>%
  summarise(count = n()) %>%
  ungroup()
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



