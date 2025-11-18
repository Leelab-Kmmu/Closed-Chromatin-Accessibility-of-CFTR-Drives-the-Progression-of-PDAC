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

#############################################################################ST-seq
setwd("/400T/wangmeiheng/spatial_transcripts")
library(openxlsx)
tumor_scRNA<-read.xlsx("GO_scRNA_malignat.xlsx", sheet = 1,rowNames = TRUE)
tumor_GSM3036911<-read.table("GO_tumor_A_GSM3036911.txt",header = T)
tumor_GSM3405534<-read.table("GO_tumor_B_GSM3405534.txt",header = T)
tumor_spatial_A<-read.table("GO_tumor_spatial_A.txt",header = T)
tumor_spatial_B<-read.table("GO_tumor_spatial_B.txt",header = T)

# 将矩阵放入列表
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
################################富集
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
###############################################################################################################################
###########################################Figure 5SC ###########################################################################
###############################################################################################################################
###########################################################spatialA&B
df<-Go_data_df1
df$count<- as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))
df$Gene_ratio<-as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))/as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[2]]}))

ggplot(df, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + # 蓝到红的渐变
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),          # 去掉网格线
    panel.border = element_rect(color = "black", fill = NA, size = 0.5) # 添加黑色边框
  )
###############################################################################################################################
###########################################Figure 5SD ###########################################################################
###############################################################################################################################
#######################################################GSM3036911 &GSM3405534
df<-Go_data_df2
df$count<- as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))
df$Gene_ratio<-as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))/as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[2]]}))

ggplot(df, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + # 蓝到红的渐变
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),          # 去掉网格线
    panel.border = element_rect(color = "black", fill = NA, size = 0.5) # 添加黑色边框
  )

###############################################################################################################################
###########################################Figure 5SA ###########################################################################
###############################################################################################################################
seurat_obj<-readRDS("/400T/wangmeiheng/FJ/ST/data.merge.rds")
df.mat<-read.table("/400T/wangmeiheng/spatial_transcripts/result/SpatialAB_decon_mtrx.txt",header = T)
n=nrow(seurat_obj@images$image@boundaries$centroids@coords)
colors <- c('#00BEC2','#3E4F9F','#A5CCE1','#E21A19','#974DA1','#F6B32B','#FA9897','#979797',"cyan", "brown","green")
all_scatterpie_plt <- suppressMessages(
  ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      data = df.mat[1:n,], ##
      ggplot2::aes(x = y,
                   y = x),
      cols = colnames(mat),
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
      data = df.mat[n+1:nrow(df.mat),], ##
      ggplot2::aes(x = y,
                   y = x),
      cols = colnames(mat),
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
###########################################Figure 5A ###########################################################################
###############################################################################################################################
##############################################scRNA
df<-Go_data_df3
df$count<- as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))
df$Gene_ratio<-as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))/as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[2]]}))

ggplot(df, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + # 蓝到红的渐变
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),          # 去掉网格线
    panel.border = element_rect(color = "black", fill = NA, size = 0.5) # 添加黑色边框
  )
########################################save data
save(matrices,file = "Go_FJ_matrices.RData")
###########################################################bulk
###############################################################################################################################
###########################################Figure 5SSA ###########################################################################
###############################################################################################################################
###############################TCGA
library(reshape2)
library(stringr)
# 1.读取数据
load("/400T/wangmeiheng/FJ/bulk/datExpr.RData")
# 2. 提取 CFTR 表达并分组
cftr_exp <- as.numeric(datExpr_all["CFTR", ])
cftr_group <- ifelse(cftr_exp > median(cftr_exp), "CFTR_High", "CFTR_Low")
# 3. 创建 design matrix
design <- model.matrix(~0 + factor(cftr_group))
colnames(design) <- levels(factor(cftr_group))
# 4. 差异分析
fit <- lmFit(datExpr_all, design)
contrast.matrix <- makeContrasts(CFTR_High - CFTR_Low, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2,trend = TRUE)
deg_results <- topTable(fit2, adjust="BH", number=Inf)
# 5. 提取差异基因
deg_sig <- deg_results[which(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1), ]
gene_up <- rownames(deg_sig[deg_sig$logFC > 1, ])
gene_down <- rownames(deg_sig[deg_sig$logFC < -1, ])
#save(deg_results,file = "/400T/wangmeiheng/FJ/bulk/deg_sig.RData")
# 6. 富集分析：转换为ENTREZ ID
gene_up_id <- bitr(gene_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_down_id <- bitr(gene_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
#save(gene_up_id, gene_down_id, file = "/400T/wangmeiheng/FJ/bulk/CFTR_group_difgene.RData")
load("/400T/wangmeiheng/FJ/bulk/CFTR_group_difgene.RData")
# 7. GO富集（下调）
go_down <- enrichGO(gene=gene_down_id$ENTREZID, OrgDb=org.Hs.eg.db,
                    ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
                    readable=TRUE)
View(go_down@result)
library(ggplot2)
df<-go_down@result[grep("matrix|structure|fibro",go_down@result$Description), c("Description" ,"GeneRatio","geneID","pvalue")]
df$count<- as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))
df$Gene_ratio<-as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[1]]}))/as.numeric(sapply(strsplit(df$GeneRatio,split='/'),function(x){x[[2]]}))
df1<-df[which(df$pvalue<0.05),]
options(digits = 3)   # 全局控制显示精度
df1$pvalue

df1$pvalue<-round(df1$pvalue,3)
ggplot(df1, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + # 蓝到红的渐变
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +  # 横坐标三位小数
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),          # 去掉网格线
    panel.border = element_rect(color = "black", fill = NA, size = 0.5) # 添加黑色边框
  )
###############################################################################################################################
###########################################Figure 5SSB ###########################################################################
###############################################################################################################################
###############################mouse
#### 1导入数据
mouse_data<-read.table("/400T/wangmeiheng/bulk_knock/result/gene_count_0.txt",header = T,row.names = 1)
#仅保留蛋白编码基因
# 读取 GTF 文件
gtf <- import("/400T/reference_genome/gencode.vM25.annotation.gff3")
# 查看前几行
head(gtf)

# 过滤非蛋白编码基因（gene_biotype 或 gene_type 字段）
protein_coding_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
library(stringr)
rownames(mouse_data) <- str_replace(rownames(mouse_data), "\\..*$", "")
protein_coding_genes$ID<- str_replace(protein_coding_genes$ID, "\\..*$", "")
# 提取基因ID或symbol
protein_coding_gene_ids <- unique(protein_coding_genes$ID)
mouse_data1<-mouse_data[protein_coding_gene_ids,]
mouse_data<-na.omit(mouse_data1)
#基因ID转换
id_to_gene <- setNames(protein_coding_genes$gene_name, protein_coding_genes$ID)
# 先去掉版本号
rownames(mouse_data) <- sub("\\..*$", "", rownames(mouse_data))

# 匹配 ID 与 gene_name
id_to_gene <- setNames(protein_coding_genes$gene_name, protein_coding_genes$ID)

# 提取对应 gene_name
gene_names <- id_to_gene[rownames(mouse_data)]

# 检查前几个结果
head(gene_names)
mouse_data$gene_name<-gene_names
mouse_data<-mouse_data[!duplicated(mouse_data$gene_name),]
rownames(mouse_data)<-mouse_data$gene_name
mouse_data<-mouse_data[,-13]
#save(mouse_data,file="/400T/wangmeiheng/FJ/bulk/house/mouse_data.RData")
load("/400T/wangmeiheng/FJ/bulk/house/mouse_data.RData")
# 2分组
cftr_group <-c(rep("RNAsh",6),rep("control",6))

# 3. 创建 design matrix
design <- model.matrix(~0 + factor(cftr_group))
colnames(design) <- levels(factor(cftr_group))
contrast.matrix <- makeContrasts(control-RNAsh, levels=design)
# 4. 差异分析
library(edgeR)

#Creates a DGEList
dge <- DGEList(counts = mouse_data,remove.zeros = T)

v <- voom(dge, design, plot=F) #会自动计算log(cpm)值
#拟合线性模型
fit <- lmFit(v, design)
#针对给定的对比计算估计系数和标准误差
fit2 <- contrasts.fit(fit,contrast.matrix)
#计算出t统计量，F统计量和差异表达倍数的对数
fit2 <- eBayes(fit2)

allDEG <- topTable(fit2, adjust="BH", number=Inf)
allDEG <- na.omit(allDEG)
padj = 0.05
foldChange= 1
diff_signif = allDEG[(allDEG$adj.P.Val < padj & abs(allDEG$logFC)>foldChange),]                    
diff_signif = diff_signif[order(diff_signif$logFC),]
# 转换
library(org.Mm.eg.db)   # Mus musculus (小鼠)
symbol_entrez <- bitr(rownames(diff_signif), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)

# 把 ENTREZID 加入 diff_signif
diff_signif$SYMBOL <- rownames(diff_signif)
diff_signif <- merge(diff_signif, symbol_entrez ,by="SYMBOL", all.x=TRUE)
# 提取差异基因
gene_up <- diff_signif[diff_signif$logFC > 1, "SYMBOL"]
gene_down <-diff_signif[diff_signif$logFC < -1, "SYMBOL"]
# 富集分析：转换为ENTREZ ID
gene_up_id <- bitr(gene_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
gene_down_id <- bitr(gene_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
#save(allDEG,file = "/400T/wangmeiheng/FJ/bulk/house/house_allDEG.RData")
#save(diff_signif ,gene_up_id, gene_down_id, file = "/400T/wangmeiheng/FJ/bulk/house/mouse_difgene.RData")
#################################################################GoGoGo
# GO富集
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
options(digits = 3)   # 全局控制显示精度
df1$pvalue

df1$pvalue<-round(df1$pvalue,3)
ggplot(df1, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + # 蓝到红的渐变
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +  # 横坐标三位小数
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    #panel.grid = element_blank(),          # 去掉网格线
    panel.border = element_rect(color = "black", fill = NA, size = 0.5) # 添加黑色边框
  )
###############################################################################################################################
###########################################Figure 5SSC ########################################################################
###############################################################################################################################
################################################ggVolcano
#######################################################TCGA
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
          legend_position = "UR") #标签位置为up right
#######################################################CFTR Knockdown
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
          legend_position = "UR") #标签位置为up right


#####################################################spatial_A
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
          legend_position = "UR") #标签位置为up right
#####################################################spatial_B
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
          legend_position = "UR") #标签位置为up right
#####################################################A_GSM3036911
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
          legend_position = "UR") #标签位置为up right

allDEG[which(rownames(allDEG)=="CFTR"),]
#####################################################B_GSM3405534
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
          legend_position = "UR") #标签位置为up right
allDEG[which(rownames(allDEG)=="CFTR"),]
###############################################################################################################################
###########################################Figure 5C ###########################################################################
###############################################################################################################################
library(ggcorrplot)
library(ggplot2)
library(scales)
library(dplyr)
library(ggplot2)
library(dplyr)
# 读取数据
df.mat <- read.table("/400T/wangmeiheng/spatial_transcripts/result/SpatialAB_decon_mtrx.txt", header = TRUE)
seurat_obj <- readRDS("/400T/wangmeiheng/FJ/ST/data.merge.rds")
mat<-df.mat[,3:ncol(mat)]
# 提取数据
Fib_frac <- seurat_obj$fibroblast_cell
cftr_val <- as.vector(seurat_obj@assays[["SCT"]]@counts["CFTR",])
# 创建数据框并过滤
plot_df <- data.frame(CFTR = cftr_val, Fib_fraction = Fib_frac) %>%
  filter(Fib_fraction > 0.1)
# 计算相关性
cor_test <- cor.test(plot_df$CFTR, plot_df$Fib_fraction, method = "spearman")
cor_label <- paste0("ρ = ", round(cor_test$estimate, 3), "\nP = ", signif(cor_test$p.value, 3))

# 创建汇总数据用于误差棒
plot_summary <- plot_df %>%
  group_by(CFTR) %>%
  summarise(
    Fib_mean = max(Fib_fraction),
    Fib_sd = sd(Fib_fraction),
    .groups = "drop"
  )

# plot
ggplot() +
  # 柱状图
  geom_col(data = plot_summary, 
           aes(x = CFTR, y = Fib_mean),
           fill = "white", color = "grey", width = 0.5) +
  # 误差棒
  geom_errorbar(data = plot_summary,
                aes(x = CFTR, ymin = Fib_mean - Fib_sd, ymax = Fib_mean + Fib_sd),
                size = 0.3, width = 0.15, color = "black") +
  # 散点
  geom_jitter(data = plot_df,
              aes(x = CFTR, y = Fib_fraction),
              width = 0.2, color = "grey", size = 0.5, alpha = 0.7) +
  # 趋势线
  geom_smooth(data = plot_df,
              aes(x = as.numeric(CFTR), y = Fib_fraction),
              method = "lm", se = TRUE, color = "red", 
              size = 0.8, alpha = 0.7, fill = NA) +
  # 相关性标注
  geom_text(aes(x = max(plot_df$CFTR) * 0.8, y = max(plot_df$Fib_fraction), 
                label = cor_label), 
            hjust = 1, vjust = 1) +
  # 主题和标签
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "CFTR expression", y = "Fibroblast proportion")
###############################################################################################################################
###########################################Figure 5D ###########################################################################
###############################################################################################################################
# 读取数据
df.mat<-read.table("/400T/wangmeiheng/spatial_transcripts/result/SpatialGEO_decon_mtrx.txt",header = T)
seurat_obj<-readRDS("/400T/wangmeiheng/FJ/ST/data.mergeABD.rds")
df.mat<-df.mat[,3:ncol(df.mat)]

# 提取数据
Fib_frac <- seurat_obj$fibroblast_cell
cftr_val <- as.vector(seurat_obj@assays[["SCT"]]@counts["CFTR",])
# 创建数据框并过滤
plot_df <- data.frame(CFTR = cftr_val, Fib_fraction = Fib_frac) %>%
  filter(Fib_fraction > 0.1)
# 计算相关性
cor_test <- cor.test(plot_df$CFTR, plot_df$Fib_fraction, method = "spearman")
cor_label <- paste0("ρ = ", round(cor_test$estimate, 3), "\nP = ", signif(cor_test$p.value, 3))

# 创建汇总数据用于误差棒
plot_summary <- plot_df %>%
  group_by(CFTR) %>%
  summarise(
    Fib_mean = max(Fib_fraction),
    Fib_sd = sd(Fib_fraction),
    .groups = "drop"
  )

# plot
ggplot() +
  # 柱状图
  geom_col(data = plot_summary, 
           aes(x = CFTR, y = Fib_mean),
           fill = "white", color = "grey", width = 0.5) +
  # 误差棒
  geom_errorbar(data = plot_summary,
                aes(x = CFTR, ymin = Fib_mean - Fib_sd, ymax = Fib_mean + Fib_sd),
                size = 0.3, width = 0.15, color = "black") +
  # 散点
  geom_jitter(data = plot_df,
              aes(x = CFTR, y = Fib_fraction),
              width = 0.2, color = "grey", size = 0.5, alpha = 0.7) +
  # 趋势线
  geom_smooth(data = plot_df,
              aes(x = as.numeric(CFTR), y = Fib_fraction),
              method = "lm", se = TRUE, color = "red", 
              size = 0.8, alpha = 0.7, fill = NA) +
  # 相关性标注
  geom_text(aes(x = max(plot_df$CFTR) * 0.8, y = max(plot_df$Fib_fraction), 
                label = cor_label), 
            hjust = 1, vjust = 1) +
  # 主题和标签
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "CFTR expression", y = "Fibroblast proportion")

