###############################################################################################################################
###########################################Figure 1 A#########################################################################
###############################################################################################################################
library(maftools)
library(maftools)
library(mclust)
library(NMF)
library(pheatmap)
setwd("/400T/wangmeiheng/TCGA/PDAC_mutation")
files <- list.files(pattern = '*.gz',recursive = TRUE)
files=files[seq(1,length(files),2)]
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

gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))
gene_top <- gene_count$gene[1:20] 
library(RColorBrewer)
vc_cols <- brewer.pal(n = 8, name = "YlOrRd")  
vc_cols <- c(
  "Missense_Mutation" = "#FF2000",   
  "Nonsense_Mutation" = "#90FF00", 
  "Frame_Shift_Del" = "#800080",      
  "Frame_Shift_Ins" = "#0000FF",      
  "In_Frame_Ins" = "#FFA000",         
  "In_Frame_Del" = "#FFFF00",        
  "Splice_Site" = "#FF00FF",         
  "Multi_Hit" = "#A52A2A"             
)
print(vc_cols)
par(family = "Arial")  
oncoplot(
  maf = all_mut, 
  colors = vc_cols,
  top = 30,
  fontSize = 0.6, 
  showTumorSampleBarcodes = FALSE
)


###############################################################################################################################
###########################################Figure 1 B#########################################################################
###############################################################################################################################

plotmafSummary(maf = all_mut, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,color = vc_cols)
###############################################################################################################################
###########################################Figure 1 C#########################################################################
###############################################################################################################################

#Draw oncostrip
oncostrip(
  maf = all_mut,
  genes = c('KRAS', "TP53", "CFTR"),
  colors = vc_cols
)

# Get the total number of samples
total_samples <- length(unique(all_mut@data$Tumor_Sample_Barcode))
gene_summary <- getGeneSummary(all_mut)
gene_summary$Mutation_Rate <- round(gene_summary$MutatedSamples / total_samples * 100, 2)
# Sort by number of mutation samples in descending order
gene_summary_sorted <- gene_summary[order(gene_summary$MutatedSamples, decreasing = TRUE), ]
# Export to Excel
write.xlsx(gene_summary_sorted, 
           file = "2_TCGA_gene_mut.xlsx", 
           rowNames = FALSE)


###############################################################################################################################
###########################################Figure 1 D/E#########################################################################
###############################################################################################################################
library(BSgenome)
#BiocManager::install("MutationalPatterns", update = FALSE)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(BSgenome.Hsapiens.UCSC.hg19)
library(MutationalPatterns)
vcf_files <- list.files(
  path = "/400T/ckn/yxa/vcf",  
  pattern = ".vcf",  
  full.names = TRUE,
  recursive = TRUE  
)
sample_name=basename(vcf_files)
vcfs <- read_vcfs_as_granges(vcf_files, sample_name, ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
plot_96_profile(mut_mat, condensed = TRUE)

#Integrate mutation site information from 18 samples
mut_mat=as.data.frame(mut_mat)
mut_mat$PDAC=apply(mut_mat,1,sum)
mut_mat1=as.data.frame(mut_mat[,19])%>%as.matrix()
colnames(mut_mat1)='PDAC'
row.names(mut_mat1)=row.names(mut_mat)
plot_96_profile(mut_mat1, condensed = TRUE)


