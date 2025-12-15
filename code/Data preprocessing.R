###############################################################################################################################
###########################################Figure 2 data processing############################################################
###############################################################################################################################
#########################ArchR
#########createArrowFiles, 
options(stringsAsFactors=F)
setwd('/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/fragments')
inputFiles<-as.character(c("/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/ATAC_1836591/fragments.tsv.gz",
                           "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/ATAC_T001835227/fragments.tsv.gz",
                           "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/ATAC_T001837519/fragments.tsv.gz",
                           "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/ATAC_T001837519_n/fragments.tsv.gz",
                           "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/GSM5589391_pancreas_SM-ADRUQ_rep1/GSM5589391_pancreas_SM-ADRUQ_rep1_fragments_final.tsv.gz",
                           "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/GSM5589392_pancreas_SM-IOBHS_rep1/GSM5589392_pancreas_SM-IOBHS_rep1_fragments_final.tsv.gz",
                           "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/GSM5589393_pancreas_SM-JF1NS_rep1/GSM5589393_pancreas_SM-JF1NS_rep1_fragments_final.tsv.gz",
                           "/400T/wangmeiheng/SPOT_ATAC/SPOT_ATAC/GSM5589394_pancreas_SM-JF1O6_rep1/GSM5589394_pancreas_SM-JF1O6_rep1_fragments_final.tsv.gz"))
#names(inputFiles)<-c("Tumor1","Tumor2",'Tumor3','Para5','Para1','Para2','Para3','Para4')


library(Cairo)
library(ArchR)#Version'1.0.2'
rhdf5::h5disableFileLocking()
set.seed(1)
addArchRThreads(threads = 23) 
addArchRGenome("hg19")
# ArrowFiles <- createArrowFiles(
#   inputFiles = inputFiles,
#   sampleNames = names(inputFiles),
#   minTSS = 0, 
#   minFrags = 0, 
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE,
#   force=TRUE
# )
# ArrowFiles<-as.character(c("Tumor1.arrow","Tumor2.arrow",'Tumor3.arrow','Para5.arrow','Para1.arrow','Para2.arrow','Para3.arrow','Para4.arrow'))
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
)
ArrowFiles<-as.character(c('ATAC1.arrow','ATAC2.arrow','ATAC3.arrow','ATAC_n.arrow','ATAC_p1.arrow','ATAC_p2.arrow','ATAC_p3.arrow','ATAC_p4.arrow'))

proj <- ArchRProject(ArrowFiles = ArrowFiles,copyArrows = F )
proj <- addIterativeLSI(proj, useMatrix = "TileMatrix", name = "IterativeLSI") 
proj <- addUMAP(proj, reducedDims = "IterativeLSI")
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "singlecellATAC",
  copyArrows = F #This is recommened so that you maintain an unaltered copy for later usage.
)
proj <- filterDoublets(ArchRProj = proj)##Remove double cells
#proj_umap <- saveArchRProject(ArchRProj = proj)
###############################################################################################################################
###########################################Figure 3 data processing############################################################
###############################################################################################################################
#Import data
load("/400T/ckn/胰腺癌/CRA001160_NT/sce.Rdata")
load("/400T/ckn/胰腺癌/GSE111672/sce.Rdata")
load("/400T/ckn/胰腺癌/GSE141017/sce.Rdata")
load("/400T/ckn/胰腺癌/GSE148673_NT_无/sce.Rdata")
load("/400T/ckn/胰腺癌/GSE154763_NT_无/sce.Rdata")
load("/400T/ckn/胰腺癌/GSE154778_TM/sce.Rdata")
load("/400T/ckn/胰腺癌/GSE158356/sce.Rdata")
load("/400T/ckn/胰腺癌/GSE162708_NT/sce.Rdata")
load("/400T/ckn/胰腺癌/GSE165399_NT/sce.Rdata")
#CRA001160
{
  library(dplyr)
  meta.data=CRA001160@meta.data
  meta.data$tissue <- ifelse(meta.data$tissue == "N", "normal", 
                             ifelse(meta.data$tissue == "T", "tumor", meta.data$tissue))
  
  table(meta.data$tissue)
  CRA001160$tissue=meta.data$tissue
  ############################################2
  #GSE111672--all-tumor
  table(GSE111672$sample)
  GSE111672$tissue='tumor'
  ############################################3
  #GSE141017
  table(GSE141017)
  GSE141017$tissue='tumor'
  ############################################4
  #GSE158356
  table(GSE158356)
  GSE158356$tissue='Metastatic'
  ############################################5
  #GSE148673
  table(GSE148673$cluster.pred)
  meta.data=GSE148673@meta.data
  # 假设你的数据框名为meta.data
  meta.data$tissue <- ifelse(meta.data$cluster.pred == "N", "normal", 
                             ifelse(meta.data$cluster.pred == "T", "tumor", meta.data$cluster.pred))
  table(meta.data$tissue);table(GSE148673$tissue)
  GSE148673$tissue=meta.data$tissue
  ############################################6
  #GSEGSE154763
  table(GSE154763$tissue)
  meta.data=GSE154763@meta.data
  # 假设你的数据框名为meta.data
  meta.data$tissue <- ifelse(meta.data$tissue == "N", "normal", 
                             ifelse(meta.data$tissue == "T", "tumor", meta.data$tissue))
  GSE154763$tissue=meta.data$tissue
  ############################################7
  #GSE154778
  table(GSE154778$tissue)
  meta.data=GSE154778@meta.data
  meta.data$tissue <- ifelse(meta.data$tissue == "Metastatic", "Metastatic", 
                             ifelse(meta.data$tissue == "primary_tumor", "tumor", meta.data$tissue))
  GSE154778$tissue=meta.data$tissue
  ############################################8
  #GSE162708
  table(GSE162708$tissue)
  ############################################9
  #GSE165399
  table(GSE165399$tissue)
  samples_name=c('CRA001160','GSE111672','GSE141017','GSE154778','GSE162708','GSE165399','GSE148673','GSE154763','GSE158356')
  scRNAlist <- list(CRA001160,GSE111672,GSE141017,GSE154778,GSE162708,GSE165399,GSE148673,GSE154763,GSE158356) 
  
  #Add sample name identifiers for each dataset
  for(i in 1:length(samples_name)){
    #Mitochondrial ratio
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
    #Ribosome ratio
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]") 
    scRNAlist[[i]][["sample"]] <- samples_name[i]
  }
  
  #Merge all scRNA-seq datasets
  scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)]) 
  table(scRNA@meta.data[["sample"]])
  #names(scRNAlist) <- samples_name
  scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)]) 
  scRNA 
  # An object of class Seurat  # 18818 features across 19738 samples within 1 assay 
  # Active assay: RNA (18818 features, 0 variable features) 
  table(scRNA$sample,scRNA$tissue) 
  # HNC01PBMC  HNC01TIL HNC10PBMC  HNC10TIL HNC20PBMC  
  #      1721      1298      1750      1383      1525  
  #  HNC20TIL     PBMC1     PBMC2   Tonsil1   Tonsil2  
  #      1148      2444      2436      3324      2709  
  save(scRNA,file = '/400T/wangmeiheng/scRNA/scRNA_raw.RData') 
  rm(CRA001160,GSE111672,GSE141017,GSE154778,GSE162708,GSE165399,GSE148673,GSE154763,GSE158356,scRNAlist)
  
  
}
{#Data quality control and cell annotation
  {
    #Calculate the proportion of Mitochondrial genes
    scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
    #Add the number of genes for each UMI in each cell to the metadata
    scRNA@meta.data$log10GenesPerUMI <- log10(scRNA@meta.data$nFeature_RNA) / log10(scRNA@meta.data$nCount_RNA)
    #Data quality control
    minGene=500
    minRNA=500
    log10 =0.83
    scRNA <- subset(scRNA, subset = nFeature_RNA > minGene  & nCount_RNA>minRNA &  log10GenesPerUMI>log10 &percent.mt < 10)
    dim(scRNA);
    #Searching for highly variable genes
    scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
    scale.genes <-  rownames(scRNA)
    scRNA <- ScaleData(scRNA, features = scale.genes)
    #PCA dimensionality reduction and principal component extraction
    scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
    head(scRNA)
    
    ElbowPlot(scRNA, ndims=40, reduction="pca")
    pc.num=1:40
    #Cell clustering
    scRNA <- FindNeighbors(scRNA, dims = pc.num) 
    scRNA <- FindClusters(scRNA, resolution = 0.6)
    table(scRNA@meta.data$seurat_clusters)
    #UMAP
    scRNA <- RunUMAP(scRNA, dims = pc.num)
    
    DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters",
            pt.size = 1.5,
            label = T,label.box = T
    )
    DimPlot(scRNA, reduction = "umap", group.by = "orig.ident",
            pt.size = 1.5,
            label = T,label.box = F
    )
    DimPlot(scRNA, reduction = "umap", group.by = "sample",
            pt.size = 1.5,
            label = T,label.box = F
    )
  }
  
  
  #removing Batch Effects
  {
    
    library(harmony)
    scRNA <- RunHarmony(scRNA,reduction = "pca",group.by.vars = "sample",reduction.save = "harmony")
    scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30,reduction.name = "umap")
    
  }
  
  #Cell annotation
  {
    markerGenes  <- c('CD3D',#T
                      'FCGR3A',#NK
                      'CD19','MS4A1','CD79A','CD79B','IGLL5',#B
                      'MKI67',#Malignant
                      'CDH5',#Endothelial
                      'DCN','LUM',#Fibroblast
                      'CHGB',#Endocrine
                      'PRSS1',#Pancreatic acinar cells
                      'EPCAM',#Epithelial
                      'ITGAX','KRT19','AMBP','TFF2',#Ductal
                      'CD1C',#DCs
                      'AIF1','CD68')#macrophage_cell
    library(stringr) 
    library(ggplot2)
    genes_to_check = unique(intersect(row.names(scRNA),markerGenes))
    genes_to_check
    P13 <- DotPlot(scRNA, features = genes_to_check,
                   assay='RNA'  )  + coord_flip()
    P14 <- VlnPlot(object = scRNA, features =genes_to_check,log =T )
    new.cluster.ids <- c("DUC_cell","DUC_cell",'Endothelial_cell','macrophage_cell','fibroblast_cell',#4
                         "T_cell","T_cell","fibroblast_cell","Epithelial_cell",'Endothelial_cell',#9
                         "Endocrin_cell","fibroblast_cell","DUC_cell","DUC_cell",'fibroblast_cell',#14
                         'Endocrin_cell','pancreatic_cell','DUC_cell','macrophage_cell','B_cell',#19
                         "DUC_cell","B_cell","macrophage_cell","Endothelial_cell",'macrophage_cell',#24
                         'fibroblast_cell','macrophage_cell','Epithelial_cell','Epithelial_cell','NK_cell',#29
                         "Malignant_cell","fibroblast_cell","macrophage_cell","Endothelial_cell",'DUC_cell',#34
                         'DUC_cell','T_cell','B_cell','DUC_cell','macrophage_cell',#39
                         "pancreatic_cell","macrophage_cell","Epithelial_cell","pancreatic_cell",'Endothelial_cell',#44
                         'DUC_cell','Epithelial_cell','DUC_cell','Epithelial_cell','macrophage_cell',#49
                         'B_cell','B_cell','T_cell'#52
    )
    scRNA@meta.data$celltype<- scRNA@meta.data$seurat_clusters
    levels(scRNA@meta.data$celltype) <- new.cluster.ids
    
    DotPlot(scRNA, features = unique(markerGenes),group.by = "celltype")+RotatedAxis()+
      scale_x_discrete("")+scale_y_discrete("")
    
    DimPlot(scRNA, reduction = "umap", group.by = "celltype",
            pt.size = 1.5,
            label = T,label.box = T
    )
    
  }
  
}

###############################################################################################################################
###########################################Figure 4 data processing############################################################
###############################################################################################################################
###########################################Microarray data download############################################################
library(dplyr)
library(sva)
library(GEOquery)
setwd("/400T/wangmeiheng/GEO")
rm(list=ls()) 
options(stringsAsFactors = F)
files<-list.files("/400T/wangmeiheng/GEO")
geo_ids <- c("GSE71729","GSE57495" ,"GSE62452","GSE28735","GSE183795","GSE91035","GSE71989")
###############################################################################1
GSE_Number <- geo_ids[1]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F) 
exp<-exprs(gset[[1]]) %>% as.data.frame() 
pdata<-pData(gset[[1]]) 
control_n<-grep("Normal-Pancreas",pdata$title)
tumor_n<-grep("Primary-Pancreas",pdata$title)
length(control_n)
length(tumor_n)
exp<-exp[,c(control_n,tumor_n)]
group_data<-c(rep("control",length(control_n)),rep("tumor",length(tumor_n)))
colnames(exp)<-paste(group_data,colnames(exp),sep="_")
save(exp,file=paste0(GSE_Number,"_exp.RData"))
###############################################################################1over

###############################################################################2
GSE_Number <- geo_ids[2]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F)
exp<-exprs(gset[[1]]) %>% as.data.frame()
gene.names <- read.table("/400T/wangmeiheng/GEO/GSE57495/GPL15048-tbl-1.txt",header = F,row.names = 1,sep = "\t")
exp$gene<-gene.names[,3]
exp<-exp[!duplicated(exp$gene),]
exp<-na.omit(exp)
rownames(exp)<-exp$gene
exp<-exp[,!names(exp) %in% "gene"]
pdata<-pData(gset[[1]]) 
length(tumor_n)
tumor_n<-grep("PDAC",pdata$title)
exp<-exp[,c(tumor_n)]
group_data<-c(rep("tumor",length(tumor_n)))
colnames(exp)<-paste(group_data,colnames(exp),sep="_")
save(exp,file=paste0(GSE_Number,"_exp.RData"))
################################################################################2over


###############################################################################3
GSE_Number <- geo_ids[3]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F)
exp<-exprs(gset[[1]]) %>% as.data.frame() 
gene.name<-read.table("/400T/wangmeiheng/GEO/GSE62452/GPL6244-tbl-1.txt",header = F,fill = T,sep='\t')
gene.name_s<-as.data.frame(sapply(strsplit(gene.name$V10,' '),function(x){x[3]}))
rownames(gene.name_s)<-gene.name$V1
colnames(gene.name_s)<-c("symbol")
gene.name_s<-gene.name_s[rownames(exp),]
exp$gene<-gene.name_s
exp<-exp[!duplicated(exp$gene),]
exp<-na.omit(exp)
rownames(exp)<-exp$gene
exp<-exp[,!names(exp) %in% "gene"]
pdata<-pData(gset[[1]]) 
control_n<-grep("Non-tumor",pdata$source_name_ch1)
tumor_n<-grep(" tumor",pdata$source_name_ch1)
length(control_n)
length(tumor_n)
exp<-exp[,c(control_n,tumor_n)]
group_data<-c(rep("control",length(control_n)),rep("tumor",length(tumor_n)))
colnames(exp)<-paste(group_data,colnames(exp),sep="_")
save(exp,file=paste0(GSE_Number,"_exp.RData"))
################################################################################3over

###############################################################################4
GSE_Number <- geo_ids[4]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F)
exp<-exprs(gset[[1]]) %>% as.data.frame()
gene.name<-read.table("/400T/wangmeiheng/GEO/GSE62452/GPL6244-tbl-1.txt",header = F,fill = T,sep='\t')
gene.name_s<-as.data.frame(sapply(strsplit(gene.name$V10,' '),function(x){x[3]}))
rownames(gene.name_s)<-gene.name$V1
colnames(gene.name_s)<-c("symbol")
gene.name_s<-gene.name_s[rownames(exp),]
exp$gene<-gene.name_s
exp<-exp[!duplicated(exp$gene),]
exp<-na.omit(exp)
rownames(exp)<-exp$gene
exp<-exp[,!names(exp) %in% "gene"]
pdata<-pData(gset[[1]]) 
control_n<-grep("nontumor",pdata$title)
tumor_n<-grep(" tumor",pdata$title)
length(control_n)
length(tumor_n)
exp<-exp[,c(control_n,tumor_n)]
group_data<-c(rep("control",length(control_n)),rep("tumor",length(tumor_n)))
colnames(exp)<-paste(group_data,colnames(exp),sep="_")
save(exp,file=paste0(GSE_Number,"_exp.RData"))
################################################################################4over

###############################################################################5
GSE_Number <- geo_ids[5]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F) 
exp<-exprs(gset[[1]]) %>% as.data.frame() 
gene.name<-read.table("/400T/wangmeiheng/GEO/GSE62452/GPL6244-tbl-1.txt",header = F,fill = T,sep='\t')
gene.name_s<-as.data.frame(sapply(strsplit(gene.name$V10,' '),function(x){x[3]}))
rownames(gene.name_s)<-gene.name$V1
colnames(gene.name_s)<-c("symbol")
gene.name_s<-gene.name_s[rownames(exp),]
exp$gene<-gene.name_s
exp<-exp[!duplicated(exp$gene),]
exp<-na.omit(exp)
rownames(exp)<-exp$gene
exp<-exp[,!names(exp) %in% "gene"]
pdata<-pData(gset[[1]]) 
control_n<-grep("Non-tumor",pdata$source_name_ch1)
tumor_n<-grep(" tumor",pdata$source_name_ch1)
length(control_n)
length(tumor_n)
exp<-exp[,c(control_n,tumor_n)]
group_data<-c(rep("control",length(control_n)),rep("tumor",length(tumor_n)))
colnames(exp)<-paste(group_data,colnames(exp),sep="_")
save(exp,file=paste0(GSE_Number,"_exp.RData"))
################################################################################5over


###############################################################################6
GSE_Number <- geo_ids[6]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F) 
exp<-exprs(gset[[1]]) %>% as.data.frame() 
gene.name<-read.table("/400T/wangmeiheng/GEO/GSE91035/GPL22763-tbl-1.txt",header = F,fill = T,sep='\t')
gene.name_s<-as.data.frame(sapply(strsplit(gene.name$V7,' '),function(x){x[1]}))
rownames(gene.name_s)<-gene.name$V1
colnames(gene.name_s)<-c("symbol")
gene.name_s<-gene.name_s[rownames(exp),]
exp$gene<-gene.name_s
exp<-exp[!duplicated(exp$gene),]
exp<-na.omit(exp)
rownames(exp)<-exp$gene
exp<-exp[,!names(exp) %in% "gene"]
pdata<-pData(gset[[1]]) 
control_n<-grep("normal pancreatic",pdata$title)
tumor_n<-grep("PDAC",pdata$title)
length(control_n)
length(tumor_n)
exp<-exp[,c(control_n,tumor_n)]
group_data<-c(rep("control",length(control_n)),rep("tumor",length(tumor_n)))
colnames(exp)<-paste(group_data,colnames(exp),sep="_")
save(exp,file=paste0(GSE_Number,"_exp.RData"))
################################################################################6over


############################################################################### 7
GSE_Number <- geo_ids[7]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F) 
exp<-exprs(gset[[1]]) %>% as.data.frame()
gene.name<-read.table("/400T/wangmeiheng/GEO/GSE71989/GPL570-tbl-1.txt",header = F,fill = T,sep='\t')
gene.name<-gene.name[grep("_at",gene.name$V1),]
gene.name_s<-as.data.frame(sapply(strsplit(gene.name$V11,' '),function(x){x[1]}))
rownames(gene.name_s)<-gene.name$V1
colnames(gene.name_s)<-c("symbol")
gene.name_s<-gene.name_s[rownames(exp),]
exp$gene<-gene.name_s
exp<-exp[!duplicated(exp$gene),]
exp<-na.omit(exp)
rownames(exp)<-exp$gene
exp<-exp[,!names(exp) %in% "gene"]
pdata<-pData(gset[[1]]) 
control_n<-grep("normal pancreatic",pdata$title)
tumor_n<-grep("PDAC",pdata$title)
length(control_n)
length(tumor_n)
exp<-exp[,c(control_n,tumor_n)]
group_data<-c(rep("control",length(control_n)),rep("tumor",length(tumor_n)))
colnames(exp)<-paste(group_data,colnames(exp),sep="_")
save(exp,file=paste0(GSE_Number,"_exp.RData"))
################################################################################7over
#####################################Bulk RNA-seq data download########################
library(GEOquery)
library(stringr)
library(reshape2)
library(limma)
library(readxl)
rm(list=ls()) 
options(stringsAsFactors = F) 
files<-list.files("/400T/wangmeiheng/GEO")
geo_ids <- c("GSE151580",  "GSE171485", "GSE114453")
###############################################################################1
GSE_Number <- geo_ids[1]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F)  
exp<-exprs(gset[[1]]) %>% as.data.frame()
exp<-read.table("GSE151580_rawCountTable.csv.gz",header=T,sep=',',row.names = 1)
keytypes(org.Hs.eg.db)
exp$gene<-mapIds(org.Hs.eg.db, 
                 keys = rownames(exp), 
                 column = "SYMBOL", 
                 keytype = "ENSEMBL", 
                 multiVals = "first")
exp<-exp[!duplicated(exp$gene),]
exp<-na.omit(exp)
rownames(exp)<-exp$gene
exp<-exp[,!names(exp) %in% "gene"]
control_n<-grep("N",colnames(exp))
tumor_n<-grep("PT",colnames(exp))
exp<-exp[,c(control_n,tumor_n)]
group_data<-c(rep("control",length(control_n)),rep("tumor",length(tumor_n)))
colnames(exp)<-paste(group_data,colnames(exp),sep="_")
save(exp,file=paste0(GSE_Number,"_exp.RData"))
################################################################################1over

###############################################################################2
GSE_Number <- geo_ids[2]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F)  
exp<-read.csv("GSE171485_PDAC-tissue-ExpressionMatrix.csv.gz")
exp<-exp[!duplicated(exp$gene_short_name),]
rownames(exp) <- exp[,2] 
exp<-exp[,c(5:ncol(exp))]
pdata<-pData(gset[[1]]) 
save(exp,file=paste0(GSE_Number,"_exp.RData"))
################################################################################2over

###############################################################################3
GSE_Number <- geo_ids[3]
print(GSE_Number)
setwd(file.path("/400T/wangmeiheng/GEO",GSE_Number))
#Data Download and Grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F) 
untar("GSE114453_RAW.tar")
fs <- list.files(pattern = "_annotated\\.txt\\.gz$")
temp_exp_list<-list()
for (i in 1:c(length(fs))) {
  temp_exp<-read.table(fs[i],header=T)
  temp_exp<-temp_exp[temp_exp$type=='protein_coding',c("type","name",'tpm')]
  temp_exp$gene<-unlist(lapply(str_split(temp_exp$name,"-"),function(x){x[[1]]}))
  temp_exp<-temp_exp[!duplicated(temp_exp$gene),]
  rownames(temp_exp)<-temp_exp$gene
  temp_exp_list[[i]]<-temp_exp
}
names(temp_exp_list)<-substr(fs[1:length(fs)],1,10)
pdata<-pData(gset[[1]]) 
pdata2<-pData(gset[[2]]) 
sample_info<-data.frame(c(rownames(pdata),rownames(pdata2)),c(pdata$source_name_ch1,pdata2$source_name_ch1))
colnames(sample_info)<-c("GSM_ID","type")
nor<-sample_info[which(sample_info$type=="Pancreas"),]
tum<-sample_info[which(sample_info$type=="Tumor"),]
which(names(temp_exp_list) %in% nor$GSM_ID)
which(names(temp_exp_list) %in% tum$GSM_ID)
temp_exp_list<-temp_exp_list[c(which(names(temp_exp_list) %in% nor$GSM_ID),
                               which(names(temp_exp_list) %in% tum$GSM_ID))]

id_lists <- lapply(temp_exp_list, function(df) df$gene)
common_ids <- Reduce(intersect, id_lists)
temp_exp_list<- lapply(temp_exp_list, function(df) df[common_ids,])
exp_list<-lapply(temp_exp_list, function(df) df$tpm)
exp<-Reduce(function(x, y) cbind(x, y), exp_list)
colnames(exp)<-names(exp_list)
rownames(exp)<-common_ids 
group_data<-c(rep("control",nrow(nor)),rep("tumor",nrow(tum)))
colnames(exp)<-paste(group_data,colnames(exp),sep="_")
save(exp,file=paste0(GSE_Number,"_exp.RData"))
################################################################################3over
#Standardise count and FPKM data types to TPM
rm(list=ls()) 
source("/home/wangmeiheng/code/Counts2TPM.R")
options(stringsAsFactors = F) 
files<-list.files("/400T/wangmeiheng/GEO")
geo_ids <- c( "GSE151580",  "GSE171485", "GSE114453")
for(i in 1:length(geo_ids) ){
  GSE_Number <- geo_ids[i]
  load(file.path("/400T/wangmeiheng/GEO",GSE_Number,paste0(GSE_Number,"_exp.RData")))
  assign(paste0(GSE_Number, "_exp"), exp)
}
geneid_efflen<-geneid_efflen[!duplicated(geneid_efflen$Symbol),]
geneid_efflen<-na.omit(geneid_efflen)
rownames(geneid_efflen)<-geneid_efflen$Symbol
###########GSE151580_exp#counts#################################1
GSE151580_exp<-GSE151580_exp[rownames(GSE151580_exp) %in% geneid_efflen$Symbol,]
geneid_efflen_temp<-geneid_efflen[rownames(GSE151580_exp),]
effLen<-geneid_efflen_temp$efflen
GSE151580_tpm <- apply(GSE151580_exp,2,function(x){Counts2TPM(x,effLen)})

##################GSE171485_exp#FPKM#############################2
GSE171485_exp<-GSE171485_exp[rownames(GSE171485_exp) %in% geneid_efflen$Symbol,]
geneid_efflen_temp<-geneid_efflen[rownames(GSE171485_exp),]
effLen<-geneid_efflen_temp$efflen
GSE171485_tpm <- apply(GSE171485_exp,2,FPKM2TPM)


#################GSE114453_exp#TPM##############################3
GSE114453_tpm<-GSE114453_exp
setwd("400T/wangmeiheng/GEO")
save(GSE248485_tpm,GSE252710_tpm,GSE255434_tpm,GSE151580_tpm,GSE171485_tpm,GSE114453_tpm,file="tpm.RData")
GSE151580_tpm_data<-log2(GSE151580_tpm+1)
marker_GSE151580_tpm_data<-GSE151580_tpm_data[markerGenes[markerGenes %in% rownames(GSE151580_tpm_data)],]

#Extract CFTR expression
markerGenes  <- c('CFTR')
GSE171485_tpm_data<-log2(GSE171485_tpm+1)
marker_GSE171485_tpm_data<-GSE171485_tpm_data[markerGenes[markerGenes %in% rownames(GSE171485_tpm_data)],]

GSE114453_tpm_data<-log2(GSE114453_tpm+1)
marker_GSE114453_tpm_data<-GSE114453_tpm_data[markerGenes[markerGenes %in% rownames(GSE114453_tpm_data)],]

save(marker_pbmc_data,GSE151580_tpm_data,GSE171485_tpm_data,GSE114453_tpm_data,file="marker_tpm.RData")

###############################################################################################################################
#################################################Remove batch effects##########################################################
###############################################################################################################################
library(dplyr)
library(sva)
rm(list=ls())
#######################################################################################Microarray data
setwd("/400T/wangmeiheng/GEO")
options(stringsAsFactors = F)  
files<-list.files("/400T/wangmeiheng/GEO",pattern = '^GSE')
geo_ids <- c("GSE71729","GSE57495" ,"GSE62452","GSE28735","GSE183795","GSE91035","GSE71989")

for(i in 1:length(geo_ids) ){
  GSE_Number <- geo_ids[i]
  load(file.path("/400T/wangmeiheng/GEO",GSE_Number,paste0(GSE_Number,"_exp.RData")))
  assign(paste0(GSE_Number, "_exp"), exp)
}

GSE71729_exp <- data.frame(Gene=rownames(GSE71729_exp), GSE71729_exp)
GSE57495_exp <- data.frame(Gene=rownames(GSE57495_exp), GSE57495_exp)
GSE62452_exp <- data.frame(Gene=rownames(GSE62452_exp), GSE62452_exp)
GSE28735_exp <- data.frame(Gene=rownames(GSE28735_exp), GSE28735_exp)
GSE183795_exp <- data.frame(Gene=rownames(GSE183795_exp),GSE183795_exp)
GSE91035_exp <- data.frame(Gene=rownames(GSE91035_exp), GSE91035_exp)
GSE71989_exp <- data.frame(Gene=rownames(GSE71989_exp), GSE71989_exp)
list_of_dataframes <- list( GSE71729_exp, GSE57495_exp,GSE62452_exp,GSE28735_exp,GSE183795_exp,GSE91035_exp,GSE71989_exp)#GSE49641_exp

#Merge all data frames using the Reduce function and full_join function
combined_exp <- Reduce(function(x, y) full_join(x, y, by = "Gene"), list_of_dataframes)
rownames(combined_exp)<-combined_exp$Gene
combined_exp<-na.omit(combined_exp)
combined_exp <- combined_exp[,-1]
dim(combined_exp)
print(head(combined_exp))
#Create batch vector
batch <- c( 
  rep("Batch1", ncol(GSE71729_exp) - 1),
  rep("Batch2", ncol(GSE57495_exp) - 1),
  rep("Batch3", ncol(GSE62452_exp) - 1),  
  rep("Batch4", ncol(GSE28735_exp) - 1),
  rep("Batch5", ncol(GSE183795_exp) - 1),
  rep("Batch6", ncol(GSE91035_exp) - 1),
  rep("Batch7", ncol(GSE71989_exp) - 1))

status <- sub('_.*$','',colnames(combined_exp))#array 787
table(status)
print(batch)
print(status)
#Correct batch effects
combined_exp_corrected <- ComBat(dat = as.matrix(combined_exp), batch = batch ,mod = model.matrix(~status))
combined_exp_corrected
rownames(combined_exp_corrected)<-rownames(combined_exp)
boxplot(colSums(combined_exp_corrected))
combined_exp_corrected<-apply(combined_exp_corrected, 2,function(x){x/sum(x)*10^5})
save(combined_exp_corrected,status,file="array_combined_exp_corrected.RData")
#######################################################################################Bulk RNA-seq
#Import data
setwd("/400T/wangmeiheng/TCGA/UCSC")
load("TCGA.RData")
load("normal_pancreas_data.RData")
exp_marker_PDAC<-tcga_PAAD_01[markerGenes,] #10 X 178
#marker_gtex_data_blood <-gtex_data_blood[markerGenes,] # 10 X 444
marker_normal_pancreas_data<-normal_pancreas_data[markerGenes,]#10 X 171
load("/400T/wangmeiheng/GEO/tpm.RData")

colnames(tcga_PAAD_01)<-paste("tumor",colnames(tcga_PAAD_01),sep="_")
colnames(normal_pancreas_data)<-paste("control",colnames(normal_pancreas_data),sep="_")
colnames(GSE171485_tpm)<-c(paste("tumor",colnames(GSE171485_tpm)[1:6],sep="_"),paste("control",colnames(GSE171485_tpm),sep="_")[7:12])


GSE151580_tpm <- data.frame(Gene=rownames(GSE151580_tpm), log2(GSE151580_tpm+0.001))
GSE171485_tpm <- data.frame(Gene=rownames(GSE171485_tpm), log2(GSE171485_tpm+0.001))
GSE114453_tpm <- data.frame(Gene=rownames(GSE114453_tpm), log2(GSE114453_tpm+0.001))

tcga_PAAD_01<- data.frame(Gene=rownames(tcga_PAAD_01), tcga_PAAD_01)
normal_pancreas_data <- data.frame(Gene=rownames(normal_pancreas_data),normal_pancreas_data)

list_of_dataframes <- list(GSE151580_tpm,GSE171485_tpm,GSE114453_tpm,tcga_PAAD_01,normal_pancreas_data)
#Merge all data frames using the Reduce function and full_join function
combined_exp <- Reduce(function(x, y) full_join(x, y, by = "Gene"),list_of_dataframes)

rownames(combined_exp)<-combined_exp$Gene
combined_exp<-na.omit(combined_exp)
#combined_exp[markerGenes,]
combined_exp <- combined_exp[,-1]
#Go to batch
#Create batch vector
batch <- c(
  rep("Batch1", ncol(GSE151580_tpm) - 1),
  rep("Batch2", ncol(GSE171485_tpm) - 1),  
  rep("Batch3", ncol(GSE114453_tpm) - 1), 
  rep("Batch4", ncol(tcga_PAAD_01)-1),
  rep("Batch5", ncol(normal_pancreas_data) - 1))
length(batch)
status <- sub('_.*$','',colnames(combined_exp))
print(batch)
print(status)

library(sva)
#Extract tumor data
tumor_data <- combined_exp[, status == "tumor"]
#Obtain the batch information corresponding to the tumour data
tumor_batch <- batch[status == "tumor"]
#Select Batch4 as the reference batch
ref_batch <- "Batch4"
#Using ComBat for batch effect correction
tumor_corrected <- ComBat(dat = as.matrix(tumor_data), batch = tumor_batch, mod = NULL, par.prior = TRUE, ref.batch = ref_batch)
#Merge the corrected tumour data back into the overall expression data matrix
combined_exp_corrected <- combined_exp
combined_exp_corrected[, status == "tumor"] <- tumor_corrected

#Extract control data
control_data <- combined_exp[, status == "control"]
#Obtain the batch information corresponding to the control data
control_batch <- batch[status == "control"]
#Chosen Batch5 as the reference batch
ref_batch <- "Batch5"
#Using ComBat for batch effect correction
control_corrected <- ComBat(dat = as.matrix(control_data), batch = control_batch, mod = NULL, par.prior = TRUE, ref.batch = ref_batch)
#Merge the corrected control data back into the overall expression data matrix
combined_exp_corrected[, status == "control"] <- control_corrected
save(combined_exp_corrected,status,file="bulk_combined_exp_corrected.RData")

###############################################################################################################################
###########################################Figure 5 data processing############################################################
###############################################################################################################################
###############################################spatial_A#######################################################################
options(stringsAsFactors=F)
#Import data
setwd('/400T/wangmeiheng/spatial_transcripts/A/outs')
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
expr.data <- Seurat::Read10X_h5(filename = 'filtered_feature_bc_matrix.h5')
spatial_A <- Seurat::CreateSeuratObject(counts = expr.data, project = 'spatial_A', assay = 'Spatial')
spatial_A$slice <- 1
spatial_A$region <- 'spatial_A'
img <- Seurat::Read10X_Image(image.dir = 'spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = spatial_A)]
spatial_A[['image']] <- img
save(spatial_A,file='/400T/wangmeiheng/spatial_transcripts/A/spatial_A.Rdata')
options(stringsAsFactors=F)
setwd('/400T/wangmeiheng/spatial_transcripts/A/')
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
load('spatial_A.Rdata')#36601 X4110
#UMI ratio
p1 <- VlnPlot(spatial_A, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(spatial_A, features = "nCount_Spatial") + theme(legend.position = "right")
#pdf('nCount_Spatial_volin.pdf')
#gene ratio
p3<-VlnPlot(spatial_A, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
p4<-SpatialFeaturePlot(spatial_A, features = "nFeature_Spatial") + theme(legend.position = "right")
#dev.off()
p1+p2
p3+p4
#Mitochondrial ratio
spatial_A[["percent.mt"]] <- PercentageFeatureSet(spatial_A, pattern = "^MT-")
p5 <- VlnPlot(spatial_A, features = "percent.mt", pt.size = 0.1) + NoLegend()
p6 <- SpatialFeaturePlot(spatial_A, features = "percent.mt") + theme(legend.position = "right")
#Data filtering
spatial_A <- subset(spatial_A, subset = nFeature_Spatial > 100 & nFeature_Spatial <7500 & nCount_Spatial > 200 & nCount_Spatial < 40000 & percent.mt < 10)#36601 X 3876
markerGenes  <- c(
  'KRT19','AMBP',"TM4SF1","S100A4",'MKI67',"MUC1",'EPCAM',#Ductal
  'DCN','LUM','COL1A1',#Fibroblast
  'CHGB',#Endocrine
  'PRSS1',#Pancreatic acinar cells
  'CD3D',#T
  'FCGR3A',#NK
  'CD79A',#B
  'AIF1','CD1C',#DC
  'CD68')#Macrophage
spatial_A<- SCTransform(spatial_A, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(spatial_A, features = c('CFTR'), pt.size.factor = 1)
SpatialPlot(spatial_A, features = "CFTR",image.alpha=0.2)
save(spatial_A,file='/400T/wangmeiheng/spatial_transcripts/A/spatial_A.Rdata')
spatial_A <- RunPCA(spatial_A, assay = "SCT", verbose = FALSE)
spatial_A <- FindNeighbors(spatial_A, reduction = "pca", dims = 1:30)
for(res in seq(0,1,0.05)){
  spatial_A <- FindClusters(spatial_A, verbose = FALSE,resolution = res)
}
spatial_A <- RunUMAP(spatial_A, reduction = "pca", dims = 1:30)
save(spatial_A,file='spatial_A_umap.Rdata')
#marker genes 
options(stringsAsFactors=F)
setwd('/400T/wangmeiheng/spatial_transcripts/A/')
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
load('spatial_A_umap.Rdata')
library(RColorBrewer)
library(clustree)
clustree(spatial_A, prefix = 'SCT_snn_res.') + coord_flip()
colors <- brewer.pal(n = 16, name = "Set3")
Idents(spatial_A) <- spatial_A$SCT_snn_res.0.3
temp_res='SCT_snn_res.0.3'
DoHeatmap(spatial_A, features = markerGenes, size = 3.5,group.colors = colors) + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill") + theme(axis.text.y = element_text(size = 9))
DotPlot(spatial_A, features = markerGenes) + scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+RotatedAxis()
VlnPlot(spatial_A,features = c('CFTR','KRT19'), group.by = temp_res,pt.size = 0)
#Manual annotation
spatial_A@meta.data$celltype<-as.character(spatial_A@meta.data$SCT_snn_res.0.35)
spatial_A@meta.data$celltype[spatial_A@meta.data$celltype=='0']<-'Stroma'
spatial_A@meta.data$celltype[spatial_A@meta.data$celltype=='1']<-'DUC'
spatial_A@meta.data$celltype[spatial_A@meta.data$celltype=='2']<-'Stroma'
spatial_A@meta.data$celltype[spatial_A@meta.data$celltype=='3']<-'Malignant DUC'
spatial_A@meta.data$celltype[spatial_A@meta.data$celltype=='4']<-'FIBROBLAST'
spatial_A@meta.data$celltype[spatial_A@meta.data$celltype=='5']<-'Normal Pancreatic Tissue'
spatial_A@meta.data$celltype[spatial_A@meta.data$celltype=='6']<-'Normal Pancreatic Tissue'
#save(spatial_A,file='spatial_A_umap.Rdata')
###############################################spatial_B#######################################################################
options(stringsAsFactors=F)
#Import data
setwd('/400T/wangmeiheng/spatial_transcripts/B/outs')
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
expr.data <- Seurat::Read10X_h5(filename = 'filtered_feature_bc_matrix.h5')
spatial_B <- Seurat::CreateSeuratObject(counts = expr.data, project = 'spatial_B', assay = 'Spatial')
spatial_B$slice <- 1
spatial_B$region <- 'spatial_B'
img <- Seurat::Read10X_Image(image.dir = 'spatial')
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = spatial_B)]
spatial_B[['image']] <- img
save(spatial_B,file='/400T/wangmeiheng/spatial_transcripts/B/spatial_B.Rdata')
options(stringsAsFactors=F)
setwd('/400T/wangmeiheng/spatial_transcripts/B/')
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
load('spatial_B.Rdata')#36601 X3899
#UMI
p1 <- VlnPlot(spatial_B, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(spatial_B, features = "nCount_Spatial") + theme(legend.position = "right")
#gene
p3<-VlnPlot(spatial_B, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
p4<-SpatialFeaturePlot(spatial_B, features = "nFeature_Spatial") + theme(legend.position = "right")
p1+p2
p3+p4
#Mitochondria
spatial_B[["percent.mt"]] <- PercentageFeatureSet(spatial_B, pattern = "^MT-")
p5 <- VlnPlot(spatial_B, features = "percent.mt", pt.size = 0.1) + NoLegend()
p6 <- SpatialFeaturePlot(spatial_B, features = "percent.mt") + theme(legend.position = "right")
plot_grid(p5, p6)
#Data filtering
spatial_B <- subset(spatial_B, subset = nFeature_Spatial > 100 & nFeature_Spatial <7000 & nCount_Spatial > 200 & nCount_Spatial < 30000 & percent.mt < 10)#36601 X 3701
markerGenes  <- c(
  'KRT19','AMBP',"TM4SF1","S100A4",'MKI67',"MUC1",'EPCAM',
  'DCN','LUM','COL1A1',
  'CHGB',
  'PRSS1',
  'CD3D',
  'FCGR3A',
  'CD79A',
  'AIF1','CD1C',
  'CD68')
spatial_B<- SCTransform(spatial_B, assay = "Spatial", verbose = FALSE)
SpatialPlot(spatial_B, features = markerGenes,image.alpha=0)
save(spatial_B,file='/400T/wangmeiheng/spatial_transcripts/spatial_B.Rdata')
spatial_B <- RunPCA(spatial_B, assay = "SCT", verbose = FALSE)
spatial_B <- FindNeighbors(spatial_B, reduction = "pca", dims = 1:30)
for(res in seq(0,1,0.05)){
  spatial_B <- FindClusters(spatial_B, verbose = FALSE,resolution = res)
}
spatial_B <- RunUMAP(spatial_B, reduction = "pca", dims = 1:30)
sample_cluster<-Idents(spatial_B)
DimPlot(spatial_B, reduction = "umap", label = TRUE)
SpatialDimPlot(spatial_B, label = TRUE, label.size = 3)
###marker genes
options(stringsAsFactors=F)
setwd('/400T/wangmeiheng/spatial_transcripts/B/')
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(RColorBrewer)
library(clustree)
clustree(spatial_B, prefix = 'SCT_snn_res.') + coord_flip()
colors <- brewer.pal(n = 16, name = "Set3")
Idents(spatial_B) <- spatial_B$SCT_snn_res.0.3
temp_res='SCT_snn_res.0.3'
SpatialPlot(spatial_B, group.by = temp_res,alpha = 0)
SpatialPlot(spatial_B, group.by = temp_res,image.alpha = 0)
SpatialFeaturePlot(spatial_B, features = c("CFTR",'KRT19'))
DoHeatmap(spatial_B, features = markerGenes, size = 3.5,group.colors = colors) + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill") + theme(axis.text.y = element_text(size = 9))
DotPlot(spatial_B, features = markerGenes) + scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+RotatedAxis()
#Manual annotation
spatial_B@meta.data$celltype<-as.character(spatial_B@meta.data$SCT_snn_res.0.3)
spatial_B@meta.data$celltype[spatial_B@meta.data$celltype=='0']<-'Stroma'
spatial_B@meta.data$celltype[spatial_B@meta.data$celltype=='1']<-'DUC'
spatial_B@meta.data$celltype[spatial_B@meta.data$celltype=='2']<-'Malignant DUC'
spatial_B@meta.data$celltype[spatial_B@meta.data$celltype=='3']<-'FIBROBLAST'
spatial_B@meta.data$celltype[spatial_B@meta.data$celltype=='4']<-'Normal Pancreatic Tissue'
spatial_B@meta.data$celltype[spatial_B@meta.data$celltype=='5']<-'Normal Pancreatic Tissue'
save(spatial_B,file='/400T/wangmeiheng/spatial_transcripts/spatial_B_umap.Rdata')
###############################################A_GSM3036911#######################################################################
library(png)
library(jsonlite)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(clustree)
library(patchwork)
rm(list = ls())
#Import data
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_A") 
files=list.files(pattern = '^GSM')
assay = "Spatial"
slice = "slice1"
temp_file<-files[3]
image <- readPNG(source = file.path(paste0(temp_file,"/spatial/example.png")))
scale.factors <- fromJSON(txt = file.path(paste0(temp_file,"/spatial/scalefactors_json.json")))
unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
spot.radius <- unnormalized.radius/max(dim(x = image))
image[,,1]<-1
image[,,2]<-1
image[,,3]<-1
file_count <- list.files(paste0("./",temp_file),pattern = ".txt|tsv")
slice_name <- unlist(strsplit(file_count,"_"))[[1]]
count_data<-read.delim(paste0(temp_file,'/',file_count),stringsAsFactors = F,check.names = F) %>% as.data.frame()  
#Data filtering
count_data<-count_data[!duplicated(count_data$Genes),] %>% as.data.frame()
rownames(count_data)<-count_data$Genes
count_data<-count_data[,-1]
for(i in 1:ncol(count_data)){count_data[,i] <- as.numeric(count_data[,i])}
tissue.positions <- data.frame(tissue=rep(1,ncol(count_data)),
                               imagerow=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[2]}))),
                               imagecol=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[1]}))))

colnames(count_data)<-paste0("cell_",1:ncol(count_data))
rownames(tissue.positions)<-colnames(count_data)
image_info<-new( Class = "VisiumV1", 
                 image = image, 
                 scale.factors = scalefactors(
                   spot = scale.factors$spot_diameter_fullres, 
                   fiducial = scale.factors$fiducial_diameter_fullres, 
                   hires = scale.factors$tissue_hires_scalef, 
                   scale.factors$tissue_lowres_scalef), 
                 coordinates = tissue.positions, 
                 spot.radius = spot.radius)
seurat_ST <- CreateSeuratObject(counts = count_data, assay = assay)
image_info <- image_info[Cells(x = seurat_ST)]
DefaultAssay(seurat_ST = image_info) <- assay
seurat_ST[[slice]] <- image_info

seurat_ST$orig.ident <- unlist(strsplit(file_count,split = "_"))[[1]]
seurat_ST <- SetIdent(seurat_ST, value = "orig.ident")
seurat_ST <- SCTransform(seurat_ST, assay = "Spatial", verbose = FALSE)
save(seurat_ST,file=paste0(temp_file,'/',temp_file,"_seurat_ST.RData"))
################################################A_GSM4100721 A_GSM4100722 are similar########################################################
rm(list = ls())
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_A") 
files=list.files(pattern ='^GSM')
load(paste0(files[1],'/',files[1],"_seurat_ST.RData"))
A_GSM3036911<-seurat_ST
load(paste0(files[2],'/',files[2],"_seurat_ST.RData"))
A_GSM4100721<-seurat_ST
load(paste0(files[3],'/',files[3],"_seurat_ST.RData"))
A_GSM4100722<-seurat_ST
##########################################A_GSM3036911################################################################################
#Dimensionality reduction, clustering
library(clustree)
A_GSM3036911 <- RunPCA(A_GSM3036911, verbose = FALSE)
ElbowPlot(A_GSM3036911)
A_GSM3036911 <- FindNeighbors(A_GSM3036911, dims = 1:30)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){A_GSM3036911<- FindClusters(A_GSM3036911, resolution = res)}
p1 <- clustree(A_GSM3036911, prefix = 'SCT_snn_res.') + coord_flip()
p1
#Cell annotation
markerGenes  <- c('CRISP3',#Ductal
                  'KRT19',"TM4SF1","S100A4",'MKI67','EPCAM',#Malignant
                  'DCN','LUM','COL1A1',#Fibroblast
                  'PRSS1'#Pancreatic acinar cells
)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(RColorBrewer)
library(clustree)
colors <- brewer.pal(n = 16, name = "Set3")
Idents(A_GSM3036911) <- A_GSM3036911$SCT_snn_res.1
temp_res='SCT_snn_res.1'
selected_genes <- union(VariableFeatures(A_GSM3036911), markerGenes)
A_GSM3036911 <- ScaleData(A_GSM3036911, assay = "SCT", features = selected_genes)
DoHeatmap(A_GSM3036911, features = markerGenes, size = 3.5,group.colors = colors) + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill") + theme(axis.text.y = element_text(size = 9))

DotPlot(A_GSM3036911, features = markerGenes) + scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+RotatedAxis()
VlnPlot(A_GSM3036911,features = c('CFTR','KRT19'), group.by = temp_res,pt.size = 0)
A_GSM3036911@meta.data$celltype<-as.character(A_GSM3036911@meta.data$SCT_snn_res.1)
A_GSM3036911@meta.data$celltype[A_GSM3036911@meta.data$celltype=='0']<-'Stroma'
A_GSM3036911@meta.data$celltype[A_GSM3036911@meta.data$celltype=='1']<-'Cancer region'
A_GSM3036911@meta.data$celltype[A_GSM3036911@meta.data$celltype=='2']<-'DUC'
A_GSM3036911@meta.data$celltype[A_GSM3036911@meta.data$celltype=='3']<-'Pancreatic tissue'
#save(A_GSM3036911,file='/400T/wangmeiheng/spatial_transcripts/PDAC_ST_A/A_GSM3036911_umap.Rdata')

###############################################B_GSM3405534#######################################################################
#Import data
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B") 
files=list.files(pattern = '^GSM')
assay = "Spatial"
slice = "slice1"
temp_file=files[3]
image <- readPNG(source = file.path(paste0(temp_file,"/spatial/example.png")))
scale.factors <- fromJSON(txt = file.path(paste0(temp_file,"/spatial/scalefactors_json.json")))
unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
spot.radius <- unnormalized.radius/max(dim(x = image))
image[,,1]<-1
image[,,2]<-1
image[,,3]<-1

file_count <- list.files(paste0("./",temp_file),pattern = ".txt|tsv")
slice_name <- unlist(strsplit(file_count,"_"))[[1]]
count_data<-read.delim(paste0(temp_file,'/',file_count),stringsAsFactors = F,check.names = F) %>% as.data.frame() 
#Data filtering
count_data<-count_data[!duplicated(count_data$Genes),] %>% as.data.frame()
rownames(count_data)<-count_data$Genes
count_data<-count_data[,-1]
for(i in 1:ncol(count_data)){count_data[,i] <- as.numeric(count_data[,i])}
count_data<-count_data[,which(apply(count_data,2,sum)!=0)]
tissue.positions <- data.frame(tissue=rep(1,ncol(count_data)),
                               imagerow=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[2]}))),
                               imagecol=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[1]}))))

colnames(count_data)<-paste0("cell_",1:ncol(count_data))
rownames(tissue.positions)<-colnames(count_data)
image_info<-new( Class = "VisiumV1", 
                 image = image, 
                 scale.factors = scalefactors(
                   spot = scale.factors$spot_diameter_fullres, 
                   fiducial = scale.factors$fiducial_diameter_fullres, 
                   hires = scale.factors$tissue_hires_scalef, 
                   scale.factors$tissue_lowres_scalef), 
                 coordinates = tissue.positions, 
                 spot.radius = spot.radius)
seurat_ST <- CreateSeuratObject(counts = count_data, assay = assay)
image_info <- image_info[Cells(x = seurat_ST)]
DefaultAssay(seurat_ST = image_info) <- assay
seurat_ST[[slice]] <- image_info

seurat_ST$orig.ident <- unlist(strsplit(file_count,split = "_"))[[1]]
seurat_ST <- SetIdent(seurat_ST, value = "orig.ident")
seurat_ST <- SCTransform(seurat_ST, assay = "Spatial", verbose = FALSE)
save(seurat_ST,file=paste0(temp_file,'/',temp_file,"_seurat_ST.RData"))
#############################B_GSM4100723 B_GSM4100724 are similar###########################################################
#Import data
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B")
files=list.files(pattern ='^GSM')
load(paste0(files[1],'/',files[1],"_seurat_ST.RData"))
B_GSM3405534<-seurat_ST
load(paste0(files[2],'/',files[2],"_seurat_ST.RData"))
B_GSM4100723<-seurat_ST
load(paste0(files[3],'/',files[3],"_seurat_ST.RData"))
B_GSM4100724<-seurat_ST
##########################################B_GSM3405534################################################################################
#Import data
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B")
files=list.files(pattern ='^GSM')
load(paste0(files[1],'/',files[1],"_seurat_ST.RData"))
B_GSM3405534<-seurat_ST
#Dimensionality reduction, clustering
library(clustree)
B_GSM3405534 <- RunPCA(B_GSM3405534, verbose = FALSE)
ElbowPlot(B_GSM3405534)
B_GSM3405534 <- FindNeighbors(B_GSM3405534, dims = 1:30)
seq <- seq(0.1, 1.5, by = 0.1)
for(res in seq){B_GSM3405534<- FindClusters(B_GSM3405534, resolution = res)}
p1 <- clustree(A_GSM3036911, prefix = 'SCT_snn_res.') + coord_flip()
p1
#Cell annotation
markerGenes  <- c('CRISP3',
                  'KRT19',"TM4SF1","S100A4",'MKI67','EPCAM',
                  'DCN','LUM','COL1A1',
                  'PRSS1')
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(RColorBrewer)
library(clustree)
colors <- brewer.pal(n = 16, name = "Set3")
Idents(B_GSM3405534) <- B_GSM3405534$SCT_snn_res.1.1
temp_res='SCT_snn_res.1.1'
SpatialPlot(B_GSM3405534, group.by = temp_res,alpha = 0)

SpatialPlot(B_GSM3405534, group.by = temp_res,image.alpha = 0)
SpatialFeaturePlot(B_GSM3405534, features = c("CFTR",'KRT19','COL1A1'))
selected_genes <- union(VariableFeatures(B_GSM3405534), markerGenes)
B_GSM3405534 <- ScaleData(B_GSM3405534, assay = "SCT", features = selected_genes)
DoHeatmap(B_GSM3405534, features = markerGenes, size = 3.5,group.colors = colors) + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')),
                        midpoint = 0, guide = "colourbar", aesthetics = "fill") + theme(axis.text.y = element_text(size = 9))
DotPlot(B_GSM3405534, features = markerGenes) + scale_colour_gradient2(low="steelblue", mid="lightgrey", high="darkgoldenrod1")+RotatedAxis()
VlnPlot(B_GSM3405534,features = c('CFTR','KRT19'), group.by = temp_res,pt.size = 0)
B_GSM3405534@meta.data$celltype<-as.character(B_GSM3405534@meta.data$SCT_snn_res.1.1)
B_GSM3405534@meta.data$celltype[B_GSM3405534@meta.data$celltype=='0']<-'Interstitium'#
B_GSM3405534@meta.data$celltype[B_GSM3405534@meta.data$celltype=='1']<-'Cancer region'
B_GSM3405534@meta.data$celltype[B_GSM3405534@meta.data$celltype=='2']<-'DUC'#
B_GSM3405534@meta.data$celltype[B_GSM3405534@meta.data$celltype=='3']<-'Interstitium'#
B_GSM3405534@meta.data$celltype[B_GSM3405534@meta.data$celltype=='4']<-'Interstitium'#
#save(B_GSM3405534,file='/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B/B_GSM3405534_umap.Rdata')

#############################################"D_GSM4100725"########################################################################################
#Import data
library(png)
library(jsonlite)
rm(list = ls())
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_D") 
files=list.files(pattern = '^GSM')
assay = "Spatial"
slice = "slice1"
image <- readPNG(source = file.path(paste0(files[1],"/spatial/example.png")))
scale.factors <- fromJSON(txt = file.path(paste0(files[1],"/spatial/scalefactors_json.json")))
unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
spot.radius <- unnormalized.radius/max(dim(x = image))
image[,,1]<-1
image[,,2]<-1
image[,,3]<-1
file_count <- list.files(paste0("./",files[1]),pattern = ".txt|tsv")
slice_name <- unlist(strsplit(file_count,"_"))[[1]]
count_data<-read.delim(paste0(files[1],'/',file_count),stringsAsFactors = F,check.names = F) %>% as.data.frame()
#Data filtering
count_data<-count_data[!duplicated(count_data$Genes),] %>% as.data.frame()
rownames(count_data)<-count_data$Genes
count_data<-count_data[,-1]
for(i in 1:ncol(count_data)){count_data[,i] <- as.numeric(count_data[,i])}
count_data<-count_data[,which(apply(count_data,2,sum)!=0)]
tissue.positions <- data.frame(tissue=rep(1,ncol(count_data)),
                               imagerow=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[2]}))),
                               imagecol=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[1]}))))

colnames(count_data)<-paste0("cell_",1:ncol(count_data))
rownames(tissue.positions)<-colnames(count_data)
image_info<-new( Class = "VisiumV1", 
                 image = image, 
                 scale.factors = scalefactors(
                   spot = scale.factors$spot_diameter_fullres, 
                   fiducial = scale.factors$fiducial_diameter_fullres, 
                   hires = scale.factors$tissue_hires_scalef, 
                   scale.factors$tissue_lowres_scalef), 
                 coordinates = tissue.positions, 
                 spot.radius = spot.radius)
seurat_ST <- CreateSeuratObject(counts = count_data, assay = assay)
image_info <- image_info[Cells(x = seurat_ST)]
DefaultAssay(seurat_ST = image_info) <- assay
seurat_ST[[slice]] <- image_info

seurat_ST$orig.ident <- unlist(strsplit(file_count,split = "_"))[[1]]
seurat_ST <- SetIdent(seurat_ST, value = "orig.ident")
seurat_ST <- SCTransform(seurat_ST, assay = "Spatial", verbose = FALSE)
#save(seurat_ST,file=paste0(files[1],'/',files[1],"_seurat_ST.RData"))
#######################################"E_GSM4100726"#################################################
#Import data
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_E") 
files=list.files(pattern = '^GSM')
assay = "Spatial"
slice = "slice1"
image <- readPNG(source = file.path(paste0(files[1],"/spatial/example.png")))
scale.factors <- fromJSON(txt = file.path(paste0(files[1],"/spatial/scalefactors_json.json")))
unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
spot.radius <- unnormalized.radius/max(dim(x = image))
image[,,1]<-1
image[,,2]<-1
image[,,3]<-1

file_count <- list.files(paste0("./",files[1]),pattern = ".txt|tsv")
slice_name <- unlist(strsplit(file_count,"_"))[[1]]
count_data<-read.delim(paste0(files[1],'/',file_count),stringsAsFactors = F,check.names = F) %>% as.data.frame() 
#Data filtering
count_data<-count_data[!duplicated(count_data$Genes),] %>% as.data.frame()
rownames(count_data)<-count_data$Genes
count_data<-count_data[,-1]
for(i in 1:ncol(count_data)){count_data[,i] <- as.numeric(count_data[,i])}
count_data <- count_data[, !grepl("\\.", colnames(count_data))]
count_data<-count_data[,which(apply(count_data,2,sum)!=0)]
tissue.positions <- data.frame(tissue=rep(1,ncol(count_data)),
                               imagerow=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[2]}))),
                               imagecol=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[1]}))))

colnames(count_data)<-paste0("cell_",1:ncol(count_data))
rownames(tissue.positions)<-colnames(count_data)
image_info<-new( Class = "VisiumV1", 
                 image = image, 
                 scale.factors = scalefactors(
                   spot = scale.factors$spot_diameter_fullres, 
                   fiducial = scale.factors$fiducial_diameter_fullres, 
                   hires = scale.factors$tissue_hires_scalef, 
                   scale.factors$tissue_lowres_scalef), 
                 coordinates = tissue.positions, 
                 spot.radius = spot.radius)
seurat_ST <- CreateSeuratObject(counts = count_data, assay = assay)
image_info <- image_info[Cells(x = seurat_ST)]
DefaultAssay(seurat_ST = image_info) <- assay
seurat_ST[[slice]] <- image_info

seurat_ST$orig.ident <- unlist(strsplit(file_count,split = "_"))[[1]]
seurat_ST <- SetIdent(seurat_ST, value = "orig.ident")
#SCTransform
seurat_ST <- SCTransform(seurat_ST, assay = "Spatial", verbose = FALSE)
#save(seurat_ST,file=paste0(files[1],'/',files[1],"_seurat_ST.RData"))
#############################################"GSM4100728"###########################################################
#Import data
setwd("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_G") 
files=list.files(pattern = '^GSM')
assay = "Spatial"
slice = "slice1"
image <- readPNG(source = file.path(paste0(files[1],"/spatial/example.png")))
scale.factors <- fromJSON(txt = file.path(paste0(files[1],"/spatial/scalefactors_json.json")))
unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
spot.radius <- unnormalized.radius/max(dim(x = image))
image[,,1]<-1
image[,,2]<-1
image[,,3]<-1
file_count <- list.files(paste0("./",files[1]),pattern = ".txt|tsv")
slice_name <- unlist(strsplit(file_count,"_"))[[1]]
count_data<-read.delim(paste0(files[1],'/',file_count),stringsAsFactors = F,check.names = F) %>% as.data.frame()  ## 最重要的输入数据
#Data filtering
count_data<-count_data[!duplicated(count_data$Genes),] %>% as.data.frame()
rownames(count_data)<-count_data$Genes
count_data<-count_data[,-1]
for(i in 1:ncol(count_data)){count_data[,i] <- as.numeric(count_data[,i])}
count_data <- count_data[, !grepl("\\.", colnames(count_data))]
count_data<-count_data[,which(apply(count_data,2,sum)!=0)]
tissue.positions <- data.frame(tissue=rep(1,ncol(count_data)),
                               imagerow=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[2]}))),
                               imagecol=as.numeric(unlist(lapply(strsplit(colnames(count_data),"x"),function(x){x[1]}))))

colnames(count_data)<-paste0("cell_",1:ncol(count_data))
rownames(tissue.positions)<-colnames(count_data)
image_info<-new( Class = "VisiumV1", 
                 image = image, 
                 scale.factors = scalefactors(
                   spot = scale.factors$spot_diameter_fullres, 
                   fiducial = scale.factors$fiducial_diameter_fullres, 
                   hires = scale.factors$tissue_hires_scalef, 
                   scale.factors$tissue_lowres_scalef), 
                 coordinates = tissue.positions, 
                 spot.radius = spot.radius)
seurat_ST <- CreateSeuratObject(counts = count_data, assay = assay)
image_info <- image_info[Cells(x = seurat_ST)]
DefaultAssay(seurat_ST = image_info) <- assay
seurat_ST[[slice]] <- image_info
seurat_ST$orig.ident <- unlist(strsplit(file_count,split = "_"))[[1]]
seurat_ST <- SetIdent(seurat_ST, value = "orig.ident")
#SCTransform
seurat_ST <- SCTransform(seurat_ST, assay = "Spatial", verbose = FALSE)
#save(seurat_ST,file=paste0(files[1],'/',files[1],"_seurat_ST.RData"))
###############################################################################################################################
###########################################Figure 4 data processing 2############################################################
###############################################################################################################################
load("/400T/wangmeiheng/spatial_transcripts/A/spatial_A.Rdata")
load("/400T/wangmeiheng/spatial_transcripts/B/spatial_B.Rdata")
#SCTransform
spatial_B <- SCTransform(spatial_B, assay = "Spatial", verbose = FALSE)
data_merge_new<-merge(spatial_A,spatial_B)
save(data_merge_new,file = "/400T/wangmeiheng/spatial_transcripts/data_merge_new.RData")

###############################################################################################################################
###########################################Figure 5 data processing############################################################
###############################################################################################################################
################################Enrichment Analysis#############################
#ST
##########################Enrichment Analysis spatial_A#########################
setwd("/400T/wangmeiheng/spatial_transcripts/")
load(paste0("/400T/wangmeiheng/spatial_transcripts/A/spatial_A_umap.Rdata"))
Idents(spatial_A) <- spatial_A@meta.data[["celltype"]]
# Conduct differential gene analysis
markers <- FindMarkers(
  object = spatial_A,
  ident.1 = "DUC", 
  ident.2 = "Malignant DUC", 
  test.use = "wilcox", 
  logfc.threshold = 0.25,
  min.pct = 0.1           
)

#Save the differential gene results as a CSV file
write.csv(markers, file = "spatial_A_markers.csv", row.names = TRUE)
library(org.Hs.eg.db)
library(clusterProfiler)
markers<-markers[markers$p_val_adj<0.05,]
#GO enrichment
geneList<-bitr(rownames(markers)[which(markers$avg_log2FC <0)],fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
GO_down<-enrichGO(geneList$ENTREZID,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "ENTREZID",
                  ont           = "BP",
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

temp<-GO_down@result[which(GO_down@result$pvalue<0.05),]
# Save the differential gene results as a txt file
write.table(temp,file="GO_tumor_spatial_A.txt",sep="\t")
##########################Enrichment Analysis spatial_B#########################
load(paste0("/400T/wangmeiheng/spatial_transcripts/B/spatial_B_umap.Rdata"))
Idents(spatial_B) <- spatial_B$celltype
# Conduct differential gene analysis
markers <- FindMarkers(
  object = spatial_B,
  ident.1 = "DUC" ,  
  ident.2 = "Malignant DUC",  
  test.use = "wilcox", 
  logfc.threshold = 0.25,  
  min.pct = 0.1           
)
# Save the differential gene results as a CSV file
write.csv(markers, file = "spatial_B_markers.csv", row.names = TRUE)
library(org.Hs.eg.db)
library(clusterProfiler)
markers<-markers[markers$p_val_adj<0.05,]
#GO enrichment
geneList<-bitr(rownames(markers)[which(markers$avg_log2FC <0)],fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
GO_down<-enrichGO(geneList$ENTREZID,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "ENTREZID",
                  ont           = "BP",
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

temp<-GO_down@result[which(GO_down@result$pvalue<0.05),]
write.table(temp,file="GO_tumor_spatial_B.txt",sep="\t")

##########################Enrichment analysis PDAC_ST_A-GSM3036911#########################
load(paste0("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_A/A_GSM3036911_umap.Rdata"))
Idents(A_GSM3036911) <- A_GSM3036911$celltype
# Conduct differential gene analysis
markers <- FindMarkers(
  object = A_GSM3036911,
  ident.1 = "DUC" , 
  ident.2 = "Cancer region",
  test.use = "wilcox", 
  logfc.threshold = 0.25,  
  min.pct = 0.1        
)
# Save the differential gene results as a CSV file
write.csv(markers, file = "GSM3036911_markers.csv", row.names = TRUE)

library(org.Hs.eg.db)
library(clusterProfiler)

markers<-markers[markers$p_val_adj<0.05,]
#GO enrichment
geneList<-bitr(rownames(markers)[which(markers$avg_log2FC <0)],fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
GO_down<-enrichGO(geneList$ENTREZID,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "ENTREZID",
                  ont           = "BP",#ALL,BP,CC,MF
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

temp<-GO_down@result[which(GO_down@result$pvalue<0.05),]
write.table(temp,file="GO_tumor_A_GSM3036911.txt",sep="\t")


#######################Enrichment analysis PDAC_ST_B_GSM3405534############################
load(paste0("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B/B_GSM3405534_umap.Rdata"))
Idents(B_GSM3405534) <- B_GSM3405534$celltype

#Conduct differential gene analysis
markers <- FindMarkers(
  object = B_GSM3405534,
  ident.1 = "DUC" , 
  ident.2 = "Cancer region", 
  test.use = "wilcox", 
  logfc.threshold = 0.25,  
  min.pct = 0.1           
)
#Save the differential gene results as a CSV file
write.csv(markers, file = "B_GSM3405534_markers.csv", row.names = TRUE)
library(org.Hs.eg.db)
library(clusterProfiler)
markers<-markers[markers$p_val_adj<0.05,]
#Go enrichment
geneList<-bitr(rownames(markers)[which(markers$avg_log2FC <0)],fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
GO_down<-enrichGO(geneList$ENTREZID,
                  OrgDb         = "org.Hs.eg.db",
                  keyType = "ENTREZID",
                  ont           = "BP",
                  pvalueCutoff  = 0.05,
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

temp<-GO_down@result[which(GO_down@result$pvalue<0.05),]
write.table(temp,file="GO_tumor_B_GSM3405534.txt",sep="\t")
#TCGA bulk
library(reshape2)
library(stringr)
rm(list = ls())
options(stringsAsFactors = F)
options(stringsAsFactors = FALSE) 
#Import data
load("/400T/wangmeiheng/TCGA/UCSC/TCGA.RData")
load("/400T/wangmeiheng/TCGA/UCSC/normal_pancreas_data.RData")
temp_n<-which(tcga_PAAD_clinical_01$histological_type=="Pancreas-Adenocarcinoma Ductal Type")
datExpr_cancer <-tcga_PAAD_01[,temp_n]
datExpr_normal <-normal_pancreas_data
scoredata <- ssgsea
ducgene <- c("CFTR")
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(rtracklayer)
# 1. Merge expression matrices
datExpr_all <- cbind(datExpr_cancer, datExpr_normal)
group <- c(rep("Tumor", ncol(datExpr_cancer)), rep("Normal", ncol(datExpr_normal)))
#Keep protein-coding genes
gtf <- import("/400T/wangmeiheng/TCGA/gencode_v43_GRCh38.p13_CHR_Comprehensive.gene.annotation.gtf")
head(gtf)
# Filter protein-coding genes
protein_coding_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
# Extract gene ID
protein_coding_gene_ids <- unique(protein_coding_genes$gene_name)
datExpr_all<-datExpr_all[protein_coding_gene_ids,]
datExpr_all<-na.omit(datExpr_all)
#save(datExpr_all,group ,file="/400T/wangmeiheng/FJ/bulk/datExpr.RData")
load("/400T/wangmeiheng/FJ/bulk/datExpr.RData")
# 2. Extract CFTR expression and group
cftr_exp <- as.numeric(datExpr_all["CFTR", ])
cftr_group <- ifelse(cftr_exp > median(cftr_exp), "CFTR_High", "CFTR_Low")
# 3. Create design matrix
design <- model.matrix(~0 + factor(cftr_group))
colnames(design) <- levels(factor(cftr_group))
# 4. Difference Analysis
fit <- lmFit(datExpr_all, design)
contrast.matrix <- makeContrasts(CFTR_High - CFTR_Low, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2,trend = TRUE)
deg_results <- topTable(fit2, adjust="BH", number=Inf)
#5. Extract differential genes
deg_sig <- deg_results[which(deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1), ]
gene_up <- rownames(deg_sig[deg_sig$logFC > 1, ])
gene_down <- rownames(deg_sig[deg_sig$logFC < -1, ])
#save(deg_results,file = "/400T/wangmeiheng/FJ/bulk/deg_sig.RData")
# 6. Enrichment Analysis
gene_up_id <- bitr(gene_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_down_id <- bitr(gene_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
#save(gene_up_id, gene_down_id, file = "/400T/wangmeiheng/FJ/bulk/CFTR_group_difgene.RData")
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
###############################################ST Deconvolution By Spotlight#######################################################################
# BiocManager::install("scran")
library(scran)
library(SPOTlight)
library(Seurat,lib.loc = "/home/wangmeiheng/R/library")
library(SeuratObject)
library(dplyr)
#Import data
options(stringsAsFactors=F)
load("/400T/wangmeiheng/scRNA/scRNA.RData")
load("/400T/wangmeiheng/spatial_transcripts/data_merge_new.RData")
seurat_obj<-data.merge
Idents(scRNA)
#filtering
scRNA$celltype <- droplevels(scRNA$celltype)
scRNA <- as.SingleCellExperiment(scRNA)
scRNA <- logNormCounts(scRNA) 
genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(scRNA)) # Remove ribosomal and mitochondrial genes
dec <- modelGeneVar(scRNA , subset.row = genes)
hvg <- getTopHVGs(dec, n = 2000)   
colLabels(scRNA) <- colData(scRNA)$celltype   
mgs <- scoreMarkers(scRNA, subset.row = genes)  # Compute marker genes
mgs_fil <- lapply(names(mgs), function(i) {   
  x <- mgs[[i]] %>% as.data.frame() 
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.5, ]
  if(nrow(x)<1){
    cat(i, "has no genes with AUC to meet the requirements\n")
    return(NULL) 
  }
  
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  return(x)
})
mgs_df <- do.call(rbind, mgs_fil)
table(mgs_df$cluster)

#spotlight deconvolution
res <- SPOTlight(
  x = scRNA, 
  y = seurat_obj@assays[["SCT"]]@data,
  n_top = 100,
  groups = as.character(scRNA$celltype), 
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

#Extract data from `SPOTlight`:
decon_mtrx <- res$mat
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, decon_mtrx)
head(seurat_obj)
#Extraction result
mat <-decon_mtrx
ct <- colnames(mat)
mat[mat < 0.1] <- 0
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
  "red", "#924900", "#db6d00", "#24ff24", "#ffff6d","cyan", "brown", "gray")
pal <- colorRampPalette(paletteMartin)(length(ct))

names(pal) <- ct
coords1 <- c(data_merge_new@images[["image"]]@coordinates[["row"]],data_merge_new@images[["image"]]@coordinates[["col"]])
coords2 <- c(data_merge_new@images[["image.spatial_B"]]@coordinates[["row"]],data_merge_new@images[["image.spatial_B"]]@coordinates[["col"]])

#Merge coordinates
coor <- rbind(coords1, coords2)
temp_data <- cbind(coor, mat[1:nrow(coor),])
temp_data <- as.data.frame(temp_data)
colnames(temp_data)[1:2]<-c("x","y")
#write.table(temp_data,"/400T/wangmeiheng/spatial_transcripts/result/SpatialAB_decon_mtrx.txt")
library(scater)
library(scran)
library(SPOTlight)
library(Seurat, lib.loc = "/home/wangmeiheng/R/library")
print(packageVersion("Seurat"))#‘4.3.0’
library(SeuratObject)
library(dplyr)
#Import data
options(stringsAsFactors=F)
load("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_A/GSM3036911/GSM3036911_seurat_ST.RData")
A_GSM3036911<-seurat_ST
load("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_B/GSM3405534/GSM3405534_seurat_ST.RData")
B_GSM3405534<-seurat_ST
load("/400T/wangmeiheng/spatial_transcripts/PDAC_ST_D/D_GSM4100725_umap.Rdata")
D_GSM4100725<-seurat_ST

library(Seurat)
#Merge data
data.merge1 <- merge(
  x = A_GSM3036911,
  y = list(B_GSM3405534, D_GSM4100725),
  add.cell.ids = c("A", "B", "D") 
)
#save(data.merge1,file="/400T/wangmeiheng/FJ/ST/data.mergeABD.rds")
load("/400T/wangmeiheng/scRNA/scRNA.RData")
seurat_obj=readRDS("/400T/wangmeiheng/FJ/ST/data.mergeABD.rds")
Idents(scRNA)
scRNA$celltype <- droplevels(scRNA$celltype)

scRNA <- as.SingleCellExperiment(scRNA)
scRNA <- logNormCounts(scRNA)  #    Feature selection
genes <- !grepl(pattern = "^RP[L|S]|MT", x = rownames(scRNA))   # Remove ribosomal and mitochondrial genes
dec <- modelGeneVar(scRNA , subset.row = genes)
hvg <- getTopHVGs(dec, n = 2000)   
colLabels(scRNA) <- colData(scRNA)$celltype   # Include cell annotation information
mgs <- scoreMarkers(scRNA, subset.row = genes)  # Compute marker genes
mgs_fil <- lapply(names(mgs), function(i) {   # Retain the most relevant marker genes
  x <- mgs[[i]] %>% as.data.frame() 
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x <- x[x$mean.AUC > 0.5, ]
  if(nrow(x)<1){
    cat(i, "has no genes with AUC to meet the requirements\n")
    return(NULL)  
  }
  # Sort the genes from highest to lowest weight
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the dataframe
  x$gene <- rownames(x)
  x$cluster <- i
  return(x)
})
mgs_df <- do.call(rbind, mgs_fil)
table(mgs_df$cluster)
#spotlight deconvolution
res <- SPOTlight(
  x = scRNA, 
  y = seurat_obj@assays[["SCT"]]@data,
  n_top = 100,
  groups = as.character(scRNA$celltype), 
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")
#Extract data from `SPOTlight`:
decon_mtrx <- res$mat
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, decon_mtrx)
head(seurat_obj)
mat <-decon_mtrx
#Pie chart plotting
ct <- colnames(mat)
mat[mat < 0.1] <- 0
#Custom colour
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff",
  "red", "#924900", "#db6d00", "#24ff24", "#ffff6d","cyan", "brown", "gray")
pal <- colorRampPalette(paletteMartin)(length(ct))

names(pal) <- ct
coor<-cbind(x=c(seurat_obj@images[["slice1"]]@coordinates[["imagerow"]],
                seurat_obj@images[["slice1.2"]]@coordinates[["imagerow"]],
                seurat_obj@images[["slice1.3"]]@coordinates[["imagerow"]]
),
y=c(seurat_obj@images[["slice1"]]@coordinates[["imagecol"]],
    seurat_obj@images[["slice1.2"]]@coordinates[["imagecol"]],
    seurat_obj@images[["slice1.3"]]@coordinates[["imagecol"]]
))

temp_data <- cbind(coor, mat[1:nrow(coor),])
temp_data <- as.data.frame(temp_data)
colnames(temp_data)[1:2]<-c("x","y")
#write.table(temp_data,"/400T/wangmeiheng/spatial_transcripts/result/SpatialGEO_decon_mtrx.txt")
###############################################################################################################################
###########################################Figure 7 data processing############################################################
###############################################################################################################################
library(reshape2)
library(stringr)
library(limma)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(rtracklayer)
rm(list = ls())
options(stringsAsFactors = F)
####################################################################################
# Import data
mouse_data<-read.table("/400T/wangmeiheng/bulk_knock/result/gene_count_0.txt",header = T,row.names = 1)
# Read GTF file
gtf <- import("/400T/reference_genome/gencode.vM25.annotation.gff3")
head(gtf)
# Filter non-protein-coding genes
protein_coding_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
library(stringr)
rownames(mouse_data) <- str_replace(rownames(mouse_data), "\\..*$", "")
protein_coding_genes$ID<- str_replace(protein_coding_genes$ID, "\\..*$", "")
# Extract gene ID or symbol
protein_coding_gene_ids <- unique(protein_coding_genes$ID)
mouse_data1<-mouse_data[protein_coding_gene_ids,]
mouse_data<-na.omit(mouse_data1)
#Gene ID conversion
id_to_gene <- setNames(protein_coding_genes$gene_name, protein_coding_genes$ID)
rownames(mouse_data) <- sub("\\..*$", "", rownames(mouse_data))
#Match ID with gene_name
id_to_gene <- setNames(protein_coding_genes$gene_name, protein_coding_genes$ID)
gene_names <- id_to_gene[rownames(mouse_data)]
head(gene_names)
mouse_data$gene_name<-gene_names
mouse_data<-mouse_data[!duplicated(mouse_data$gene_name),]
rownames(mouse_data)<-mouse_data$gene_name
mouse_data<-mouse_data[,-13]
#mouse_data<-mouse_data[,-c(5)]
#save(mouse_data,file="/400T/wangmeiheng/FJ/bulk/house/mouse_data.RData")
load("/400T/wangmeiheng/FJ/bulk/house/mouse_data.RData")
cftr_group <-c(rep("RNAsh",6),rep("control",6))
design <- model.matrix(~0 + factor(cftr_group))
colnames(design) <- levels(factor(cftr_group))
contrast.matrix <- makeContrasts(control-RNAsh, levels=design)
library(edgeR)

#Creates a DGEList
dge <- DGEList(counts = mouse_data,remove.zeros = T)
v <- voom(dge, design, plot=F) 
#Fit a linear model
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)

allDEG <- topTable(fit2, adjust="BH", number=Inf)
allDEG <- na.omit(allDEG)
padj = 0.05
foldChange= 1
diff_signif = allDEG[(allDEG$adj.P.Val < padj & abs(allDEG$logFC)>foldChange),]                    
diff_signif = diff_signif[order(diff_signif$logFC),]
library(org.Mm.eg.db) 
symbol_entrez <- bitr(rownames(diff_signif), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)

# Add ENTREZID to diff_signif
diff_signif$SYMBOL <- rownames(diff_signif)
diff_signif <- merge(diff_signif, symbol_entrez ,by="SYMBOL", all.x=TRUE)
# Extract differential genes
gene_up <- diff_signif[diff_signif$logFC > 1, "SYMBOL"]
gene_down <-diff_signif[diff_signif$logFC < -1, "SYMBOL"]
gene_up_id <- bitr(gene_up, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
gene_down_id <- bitr(gene_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
#save(allDEG,file = "/400T/wangmeiheng/FJ/bulk/house/house_allDEG.RData")
#save(diff_signif ,gene_up_id, gene_down_id, file = "/400T/wangmeiheng/FJ/bulk/house/mouse_difgene.RData")
#################################################################Go enrichment
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
ggplot(df1, aes(x = Gene_ratio, y = reorder(Description, Gene_ratio), size = count, color = pvalue)) +
  geom_point() +
  labs(title = "GO Enrichment", x = "Gene Ratio", y = "Description", color = "pvalue", size = "Count") +
  scale_color_gradientn(colors = c("#D73027", "#4575B4")) + 
  scale_x_continuous(labels = function(x) sprintf("%.3f", x)) +  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

#save(df,file = "/400T/wangmeiheng/FJ/bulk/house/df.RData")
