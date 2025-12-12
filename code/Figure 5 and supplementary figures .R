library(WGCNA)
library(reshape2)
library(stringr)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(rtracklayer)
rm(list = ls())
options(stringsAsFactors = F)

#Import enrichment results
setwd("/400T/wangmeiheng/spatial_transcripts")
library(openxlsx)
tumor_scRNA<-read.xlsx("GO_scRNA_malignat.xlsx", sheet = 1,rowNames = TRUE)
tumor_GSM3036911<-read.table("GO_tumor_A_GSM3036911.txt",header = T)
tumor_GSM3405534<-read.table("GO_tumor_B_GSM3405534.txt",header = T)
tumor_spatial_A<-read.table("GO_tumor_spatial_A.txt",header = T)
tumor_spatial_B<-read.table("GO_tumor_spatial_B.txt",header = T)

matrices <- list(
  tumor_scRNA[,c(1,2,3,8,9,11)],
  tumor_GSM3036911[,c(1,2,3,8,9,11)],
  tumor_GSM3405534[,c(1,2,3,8,9,11)],
  tumor_spatial_A[,c(1,2,3,8,9,11)],
  tumor_spatial_B[,c(1,2,3,8,9,11)]
)
names(matrices)<-c("tumor_scRNA","tumor_GSM3036911",'tumor_GSM3405534','tumor_spatial_A','tumor_spatial_B')
names(matrices[[1]])<-c("ID","Description" ,"GeneRatio","pvalue","p.adjust","geneID")
names(matrices[[2]])<-c("ID","Description" ,"GeneRatio","pvalue","p.adjust","geneID")
names(matrices[[3]])<-c("ID","Description" ,"GeneRatio","pvalue","p.adjust","geneID")
names(matrices[[4]])<-c("ID","Description" ,"GeneRatio","pvalue","p.adjust","geneID")
names(matrices[[5]])<-c("ID","Description" ,"GeneRatio","pvalue","p.adjust","geneID")

Go_data<-do.call(rbind, matrices)
###############################spatialA&B
shared_pathway<-intersect(matrices$tumor_spatial_A$Description,matrices$tumor_spatial_B$Description)
shared_ECM_pathway<-grep("structure|coll|matrix|fibro",shared_pathway,value = T)
Go_data_df1<-Go_data[intersect(which(Go_data$Description %in% shared_ECM_pathway),grep("spatial_A",rownames(Go_data))) ,]
###############################GSM3036911 &GSM3405534
shared_pathway<-intersect(matrices$tumor_GSM3036911$Description,matrices$tumor_GSM3405534$Description)
shared_ECM_pathway<-grep("structure|coll|matrix|fibro",shared_pathway,value = T)
Go_data_df2<-Go_data[intersect(which(Go_data$Description %in% shared_ECM_pathway),grep("GSM3036911",rownames(Go_data))) ,]
###############################scRNA
shared_pathway<-matrices$tumor_scRNA$Description
shared_ECM_pathway<-grep("structure|coll|matrix|fibro",shared_pathway,value = T)
Go_data_df3<-Go_data[intersect(which(Go_data$Description %in% shared_ECM_pathway),grep("scRNA",rownames(Go_data))) ,]
#Enrichment results visualisation
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
###############################################################################################################################
###########################################Figure 5SC ###########################################################################
###############################################################################################################################
#spatialA&B Enrichment dotplot
df<-Go_data_df1
df$count<- as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))
df$Gene_ratio<-as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))/as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[2]]}))

ggplot(df, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),        
    panel.border = element_rect(color = "black", fill = NA, size = 0.5) 
  )
###############################################################################################################################
###########################################Figure 5SD ###########################################################################
###############################################################################################################################
#GSM3036911 &GSM3405534 Enrichment dotplot
df<-Go_data_df2
df$count<- as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))
df$Gene_ratio<-as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))/as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[2]]}))

ggplot(df, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + #
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),         
    panel.border = element_rect(color = "black", fill = NA, size = 0.5) 
  )

###############################################################################################################################
###########################################Figure 5SA #########################################################################
###############################################################################################################################
#The proportion of cells for each spot in ST data
load("/400T/wangmeiheng/spatial_transcripts/data_merge_new.RData")
#seurat_obj <- readRDS("/400T/wangmeiheng/FJ/ST/data.merge.rds")
seurat_obj<-data.merge
df.mat<-read.table("/400T/wangmeiheng/spatial_transcripts/result/SpatialAB_decon_mtrx.txt",header = T)
n=length(data.merge@images[["image"]]@coordinates[["imagerow"]])
colors <- c('#00BEC2','#3E4F9F','#A5CCE1','#E21A19','#974DA1','#F6B32B','#FA9897','#979797',"cyan", "brown","green")
all_scatterpie_plt <- suppressMessages(
  ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      data = df.mat[1:n,], 
      ggplot2::aes(x = y,
                   y = x),
      cols = setdiff(colnames(df.mat), c("x", "y")),
      color = NA,
      alpha = 0.5,
      pie_scale = 0.3) + scale_fill_manual(values = colors)+
    coord_cartesian()+
    ggplot2::scale_y_reverse()+
    cowplot::theme_half_open(11, rel_small = 1) +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed(ratio = 1,
                         xlim = NULL,
                         ylim = NULL,
                         expand = TRUE,
                         clip = "on"))
all_scatterpie_plt

all_scatterpie_plt2 <- suppressMessages(
  ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      data = df.mat[n+1:nrow(df.mat),], 
      ggplot2::aes(x = y,
                   y = x),
      cols = setdiff(colnames(df.mat), c("x", "y")),
      color = NA,
      alpha = 0.5,
      pie_scale = 0.3) + scale_fill_manual(values = colors)+
    coord_cartesian()+
    ggplot2::scale_y_reverse()+
    cowplot::theme_half_open(11, rel_small = 1) +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed(ratio = 1,
                         xlim = NULL,
                         ylim = NULL,
                         expand = TRUE,
                         clip = "on"))
all_scatterpie_plt2



###############################################################################################################################
###########################################Figure 5SB ###########################################################################
###############################################################################################################################
#The proportion of cells for each spot in ST data
seurat_obj<-readRDS("/400T/wangmeiheng/FJ/ST/data.mergeABD.rds")
df.mat<-read.table("/400T/wangmeiheng/spatial_transcripts/result/SpatialGEO_decon_mtrx.txt",header = T)
n1<-length(seurat_obj@images[["slice1"]]@coordinates[["imagerow"]])
n2<-length(seurat_obj@images[["slice1.2"]]@coordinates[["imagerow"]])
colors <- c('#00BEC2','#3E4F9F','#A5CCE1','#E21A19','#974DA1','#F6B32B','#FA9897','#979797',"cyan", "brown","green")

all_scatterpie_pltA <- suppressMessages(
  ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      data = temp_data[1:n1,], ##
      ggplot2::aes(x = y,
                   y = x),
      cols = colnames(mat),
      color = NA,
      alpha = 0.5,
      pie_scale = 0.8) + scale_fill_manual(values = colors)+
    coord_cartesian()+
    ggplot2::scale_y_reverse()+
    cowplot::theme_half_open(11, rel_small = 1) +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed(ratio = 1,
                         xlim = NULL,
                         ylim = NULL,
                         expand = TRUE,
                         clip = "on"))
all_scatterpie_pltA

all_scatterpie_pltB <- suppressMessages(
  ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      data = temp_data[429:652,], ##
      ggplot2::aes(x = y,
                   y = x),
      cols = colnames(mat),
      color = NA,
      alpha = 0.5,
      pie_scale = 1.3) + scale_fill_manual(values = colors)+
    coord_cartesian()+
    ggplot2::scale_y_reverse()+
    cowplot::theme_half_open(11, rel_small = 1) +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed(ratio = 1,
                         xlim = NULL,
                         ylim = NULL,
                         expand = TRUE,
                         clip = "on"))
all_scatterpie_pltB
all_scatterpie_pltD <- suppressMessages(
  ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      data = temp_data[653:nrow(temp_data),], ##
      ggplot2::aes(x = y,
                   y = x),
      cols = colnames(mat),
      color = NA,
      alpha = 0.5,
      pie_scale = 1) + scale_fill_manual(values = colors)+
    coord_cartesian()+
    ggplot2::scale_y_reverse()+
    cowplot::theme_half_open(11, rel_small = 1) +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed(ratio = 1,
                         xlim = NULL,
                         ylim = NULL,
                         expand = TRUE,
                         clip = "on"))
all_scatterpie_pltD

###############################################################################################################################
###########################################Figure 5A ##########################################################################
###############################################################################################################################
#scRNA Enrichment dotplot
df<-Go_data_df3
df$count<- as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))
df$Gene_ratio<-as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))/as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[2]]}))

ggplot(df, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),        
    panel.border = element_rect(color = "black", fill = NA, size = 0.5) 
  )
#save data
save(matrices,file = "Go_FJ_matrices.RData")
#bulk Enrichment dotplot
###############################################################################################################################
###########################################Figure 5SSA ###########################################################################
###############################################################################################################################
#TCGA Enrichment dotplot
library(reshape2)
library(stringr)
load("/400T/wangmeiheng/FJ/bulk/CFTR_group_difgene.RData")
go_down <- enrichGO(gene=gene_down_id$ENTREZID, OrgDb=org.Hs.eg.db,
                    ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
                    readable=TRUE)
View(go_down@result)
library(ggplot2)
df<-go_down@result[grep("matrix|structure|fibro",go_down@result$Description), c("Description" ,"GeneRatio","geneID","pvalue")]
df$count<- as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))
df$Gene_ratio<-as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))/as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[2]]}))
df1<-df[which(df$pvalue<0.05),]
options(digits = 3)   
df1$pvalue

df1$pvalue<-round(df1$pvalue,3)
#Go Enrichment
ggplot(df1, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + 
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),       
    panel.border = element_rect(color = "black", fill = NA, size = 0.5) 
  )
###############################################################################################################################
###########################################Figure 5SSB ###########################################################################
###############################################################################################################################
#mouse Enrichment dotplot
load("/400T/wangmeiheng/FJ/bulk/house/mouse_difgene.RData")
go_down <- enrichGO(gene=gene_down_id$ENTREZID, OrgDb=org.Mm.eg.db,
                    ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=1,
                    readable=TRUE)

View(go_down@result)
library(ggplot2)
df<-go_down@result[grep("matrix|structure|fibroblast",go_down@result$Description), c("Description" ,"GeneRatio","geneID","pvalue")]
df$count<- as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))
df$Gene_ratio<-as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))/as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[2]]}))
df1<-df[which(df$pvalue<0.05),]
options(digits = 3)  
df1$pvalue
df1$pvalue<-round(df1$pvalue,3)
#Go Enrichment
ggplot(df1, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + 
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),         
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )
###############################################################################################################################
###########################################Figure 5SSC ########################################################################
###############################################################################################################################
#Differential gene analysis results of multi-omics data--Volcano plot
#TCGA
library(ggVolcano)
load("/400T/wangmeiheng/FJ/bulk/deg_sig.RData")
allDEG=deg_results
allDEG$geneName <- rownames(allDEG)
allDEG <- na.omit(allDEG)
allDEG$regulate <- ifelse(allDEG$adj.P.Val > 0.05, "unchanged",
                          ifelse(allDEG$logFC > 1, "up-regulated",
                                 ifelse(allDEG$logFC < -1, "down-regulated", "unchanged")))
ggvolcano(data = allDEG,x = "logFC",y = "adj.P.Val",output = FALSE,label_number = 0,#label = "Genes",
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10P.Value",
          legend_position = "UR") 
#mouse 
load("/400T/wangmeiheng/FJ/bulk/house/house_allDEG.RData")
library(ggVolcano)
allDEG$geneName <- rownames(allDEG)
allDEG <- na.omit(allDEG)
allDEG$regulate <- ifelse(allDEG$adj.P.Val > 0.05, "unchanged",
                          ifelse(allDEG$logFC > 1, "up-regulated",
                                 ifelse(allDEG$logFC < -1, "down-regulated", "unchanged")))
ggvolcano(data = allDEG,x = "logFC",y = "adj.P.Val",output = FALSE,label_number = 0,#label = "Genes",
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10P.Value",
          legend_position = "UR") 


#spatial_A
allDEG <- read.csv("/400T/wangmeiheng/spatial_transcripts/result/A/spatial_A_markers.csv",row.names = 1)
allDEG$geneName <- rownames(allDEG)
allDEG <- na.omit(allDEG)
allDEG$regulate <- ifelse(allDEG$p_val_adj > 0.05, "unchanged",
                          ifelse(allDEG$avg_log2FC > 1, "up-regulated",
                                 ifelse(allDEG$avg_log2FC < -1, "down-regulated", "unchanged")))
ggvolcano(data = allDEG,x = "avg_log2FC",y = "p_val_adj",output = FALSE,label_number = 0,#label = "Genes",
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10P.Value",
          legend_position = "UR")
#spatial_B
allDEG <- read.csv("/400T/wangmeiheng/spatial_transcripts/result/B/spatial_B_markers.csv",row.names = 1)
allDEG$geneName <- rownames(allDEG)
allDEG <- na.omit(allDEG)
allDEG$regulate <- ifelse(allDEG$p_val_adj > 0.05, "unchanged",
                          ifelse(allDEG$avg_log2FC > 1, "up-regulated",
                                 ifelse(allDEG$avg_log2FC < -1, "down-regulated", "unchanged")))
ggvolcano(data = allDEG,x = "avg_log2FC",y = "p_val_adj",output = FALSE,label_number = 0,#label = "Genes",
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10P.Value",
          legend_position = "UR") 
#A_GSM3036911
allDEG <- read.csv("/400T/wangmeiheng/spatial_transcripts/result/PDAC_ST_A/GSM3036911_markers.csv",row.names = 1)
allDEG$geneName <- rownames(allDEG)
allDEG <- na.omit(allDEG)
allDEG$regulate <- ifelse(allDEG$p_val > 0.05, "unchanged",
                          ifelse(allDEG$avg_log2FC > 1, "up-regulated",
                                 ifelse(allDEG$avg_log2FC < -1, "down-regulated", "unchanged")))
ggvolcano(data = allDEG,x = "avg_log2FC",y = "p_val",output = FALSE,label_number = 0,#label = "Genes",
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10P.Value",
          legend_position = "UR") 

allDEG[which(rownames(allDEG)=="CFTR"),]
#B_GSM3405534
allDEG <- read.csv("/400T/wangmeiheng/spatial_transcripts/result/PDAC_ST_B/B_GSM3405534_markers.csv",row.names = 1)
allDEG$geneName <- rownames(allDEG)
allDEG <- na.omit(allDEG)
allDEG$regulate <- ifelse(allDEG$p_val > 0.05, "unchanged",
                          ifelse(allDEG$avg_log2FC > 1, "up-regulated",
                                 ifelse(allDEG$avg_log2FC < -1, "down-regulated", "unchanged")))
ggvolcano(data = allDEG,x = "avg_log2FC",y = "p_val",output = FALSE,label_number = 0,#label = "Genes",
          fills = c("#00AFBB", "#999999", "#FC4E07"),
          colors = c("#00AFBB", "#999999", "#FC4E07"),
          x_lab = "log2FC",
          y_lab = "-Log10P.Value",
          legend_position = "UR") 
allDEG[which(rownames(allDEG)=="CFTR"),]
###############################################################################################################################
###########################################Figure 5C ##########################################################################
###############################################################################################################################
library(ggcorrplot)
library(ggplot2)
library(scales)
library(dplyr)
library(ggplot2)
library(dplyr)
#Correlation between CFTR expression and the proportion of fibroblasts
#Import data
df.mat <- read.table("/400T/wangmeiheng/spatial_transcripts/result/SpatialAB_decon_mtrx.txt", header = TRUE)
load("/400T/wangmeiheng/spatial_transcripts/data_merge_new.RData")
seurat_obj<-data_merge_new
mat<-df.mat[,3:ncol(df.mat)]
Fib_frac <- seurat_obj$fibroblast_cell
cftr_val <- as.vector(seurat_obj@assays[["SCT"]]@counts["CFTR",])
# Filter spots where the proportion of fibroblasts is less than 0.1
plot_df <- data.frame(CFTR = cftr_val, Fib_fraction = Fib_frac) %>%
  filter(Fib_fraction > 0.1)
# Calculate correlation
cor_test <- cor.test(plot_df$CFTR, plot_df$Fib_fraction, method = "spearman")
cor_label <- paste0("ρ = ", round(cor_test$estimate, 3), "\nP = ", signif(cor_test$p.value, 3))

plot_summary <- plot_df %>%
  group_by(CFTR) %>%
  summarise(
    Fib_mean = max(Fib_fraction),
    Fib_sd = sd(Fib_fraction),
    .groups = "drop"
  )

# plot
ggplot() +
  geom_col(data = plot_summary, 
           aes(x = CFTR, y = Fib_mean),
           fill = "white", color = "grey", width = 0.5) +
  geom_errorbar(data = plot_summary,
                aes(x = CFTR, ymin = Fib_mean - Fib_sd, ymax = Fib_mean + Fib_sd),
                size = 0.3, width = 0.15, color = "black") +
  geom_jitter(data = plot_df,
              aes(x = CFTR, y = Fib_fraction),
              width = 0.2, color = "grey", size = 0.5, alpha = 0.7) +
  geom_smooth(data = plot_df,
              aes(x = as.numeric(CFTR), y = Fib_fraction),
              method = "lm", se = TRUE, color = "red", 
              size = 0.8, alpha = 0.7, fill = NA) +
  geom_text(aes(x = max(plot_df$CFTR) * 0.8, y = max(plot_df$Fib_fraction), 
                label = cor_label), 
            hjust = 1, vjust = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "CFTR expression", y = "Fibroblast proportion")
###############################################################################################################################
###########################################Figure 5D ###########################################################################
###############################################################################################################################
#Correlation between CFTR expression and the proportion of fibroblasts
#Import data
df.mat<-read.table("/400T/wangmeiheng/spatial_transcripts/result/SpatialGEO_decon_mtrx.txt",header = T)
seurat_obj<-readRDS("/400T/wangmeiheng/FJ/ST/data.mergeABD.rds")
df.mat<-df.mat[,3:ncol(df.mat)]
Fib_frac <- seurat_obj$fibroblast_cell
cftr_val <- as.vector(seurat_obj@assays[["SCT"]]@counts["CFTR",])
# Filter spots where the proportion of fibroblasts is less than 0.1
plot_df <- data.frame(CFTR = cftr_val, Fib_fraction = Fib_frac) %>%
  filter(Fib_fraction > 0.1)
# Calculate correlation
cor_test <- cor.test(plot_df$CFTR, plot_df$Fib_fraction, method = "spearman")
cor_label <- paste0("ρ = ", round(cor_test$estimate, 3), "\nP = ", signif(cor_test$p.value, 3))
plot_summary <- plot_df %>%
  group_by(CFTR) %>%
  summarise(
    Fib_mean = max(Fib_fraction),
    Fib_sd = sd(Fib_fraction),
    .groups = "drop"
  )

# plot
ggplot() +
  geom_col(data = plot_summary, 
           aes(x = CFTR, y = Fib_mean),
           fill = "white", color = "grey", width = 0.5) +
  geom_errorbar(data = plot_summary,
                aes(x = CFTR, ymin = Fib_mean - Fib_sd, ymax = Fib_mean + Fib_sd),
                size = 0.3, width = 0.15, color = "black") +
  geom_jitter(data = plot_df,
              aes(x = CFTR, y = Fib_fraction),
              width = 0.2, color = "grey", size = 0.5, alpha = 0.7) +
  geom_smooth(data = plot_df,
              aes(x = as.numeric(CFTR), y = Fib_fraction),
              method = "lm", se = TRUE, color = "red", 
              size = 0.8, alpha = 0.7, fill = NA) +
  geom_text(aes(x = max(plot_df$CFTR) * 0.8, y = max(plot_df$Fib_fraction), 
                label = cor_label), 
            hjust = 1, vjust = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "CFTR expression", y = "Fibroblast proportion")

