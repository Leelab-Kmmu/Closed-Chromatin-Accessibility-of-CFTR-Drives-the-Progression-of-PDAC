#######################################
#数据预处理
#########################################################################################################################
############################################1
#CRA001160
{
  library(dplyr)
  meta.data=CRA001160@meta.data
  # 假设你的数据框名为meta.data
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
  
  
  table(scRNA$sample)
  samples_name=c('CRA001160','GSE111672','GSE141017','GSE154778','GSE162708','GSE165399','GSE148673','GSE154673','GSE158356')
  scRNAlist <- list(CRA001160,GSE111672,GSE141017,GSE154778,GSE162708,GSE165399,GSE148673,GSE154673,GSE158356) 
  
  # 先为每个数据集添加样本名称标识
  for(i in 1:length(samples_name)){
    # 计算线粒体基因比例
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
    # 计算核糖体基因比例
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]") 
    # 计算红细胞基因比例
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")   
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))   
    scRNAlist[[i]][["percent.HB"]] <- PercentageFeatureSet(scRNAlist[[i]], features=HB.genes)
    
    # 为每个数据集添加样本名称标识
    scRNAlist[[i]][["sample"]] <- samples_name[i]
  }
  
  # 合并所有数据集
  scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)]) 
  
  # 查看结果（确认样本名称已添加）
  head(scRNA@meta.data$sample)
  
  
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
  save(scRNA,file = '/400T/ckn/胰腺癌/scRNA_temp.Rdata') 
  #scRNAlist <- SplitObject(scRNA, split.by = "orig.ident") #分割Seurat对象
  rm(CRA001160,GSE111672,GSE141017,GSE154778,GSE162708,GSE165399,GSE148673,GSE154673,GSE158356,scRNAlist)
  
  
}
####QC质控
{
  {
    ##计算质控指标
    #计算细胞中核糖体基因比例
    #test<-scRNA[["percent.mt"]]
    scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
    
    #计算红细胞比例
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    #test<-scRNA@assays$RNA
    #table(rownames(scRNA@assays$RNA))
    
    HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
    HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
    HB.genes <- HB.genes[!is.na(HB.genes)] 
    scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes)
    # 将每个细胞每个UMI的基因数目添加到元数据中
    scRNA@meta.data$log10GenesPerUMI <- log10(scRNA@meta.data$nFeature_RNA) / log10(scRNA@meta.data$nCount_RNA)
    
    ##数据质控
    minGene=500
    minRNA=500
    log10 =0.83
    scRNA <- subset(scRNA, subset = nFeature_RNA > minGene  & nCount_RNA>minRNA & log10GenesPerUMI>log10)
    dim(scRNA);
    
    # 数据标准化、寻找可变特征、数据缩放
    #scRNA = scRNA %>%
    #  NormalizeData() %>%
    #  FindVariableFeatures() %>%
    #  ScaleData()
    #寻找高变基因
    scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
    
    #数据中心化
    ##如果内存足够最好对所有基因进行中心化
    scale.genes <-  rownames(scRNA)
    scRNA <- ScaleData(scRNA, features = scale.genes)
    
    #PCA降维并提取主成分
    scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
    head(scRNA)
    
    ElbowPlot(scRNA, ndims=40, reduction="pca")
    pc.num=1:40
    
    ##细胞聚类
    
    scRNA <- FindNeighbors(scRNA, dims = pc.num) 
    scRNA <- FindClusters(scRNA, resolution = 0.6)
    table(scRNA@meta.data$seurat_clusters)
    ##非线性降维
    #tSNE
    
    scRNA = RunTSNE(scRNA, dims = pc.num)
    
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
  
  
  ###去批次效应
  {
    
    library(harmony)
    scRNA <- RunHarmony(scRNA,reduction = "pca",group.by.vars = "sample",reduction.save = "harmony")
    scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30,reduction.name = "umap")
    
  }
  
  ###细胞注释
  ##人工注释
  {
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
    library(stringr) 
    library(ggplot2)
    genes_to_check = unique(intersect(row.names(scRNA),markerGenes))
    genes_to_check
    P13 <- DotPlot(scRNA, features = genes_to_check,
                   assay='RNA'  )  + coord_flip()
    P14 <- VlnPlot(object = scRNA, features =genes_to_check,log =T )
    #P15 <- FeaturePlot(object = scRNA, features=genes_to_check )
    #ggsave(filename = "E:/胰腺癌/GSE154763/mk.pdf",P14)
    
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
    levels(scRNA@meta.data$celltype) <- new.cluster.ids#将celltype确定
    
    DotPlot(scRNA, features = unique(markerGenes),group.by = "celltype")+RotatedAxis()+
      scale_x_discrete("")+scale_y_discrete("")
    
    DimPlot(scRNA, reduction = "umap", group.by = "celltype",
            pt.size = 1.5,
            label = T,label.box = T
    )
    
  }
  
}
