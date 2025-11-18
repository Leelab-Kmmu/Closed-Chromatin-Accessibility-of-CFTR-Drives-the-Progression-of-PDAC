###############################################################################################################################
###########################################Figure 3 ###########################################################################
###############################################################################################################################
setwd("/400T/wangmeiheng/GEO")
###############################################################################################################################
###########################################Figure 3A ###########################################################################
###############################################################################################################################
load("array_combined_exp_corrected.RData")
# 查看结果
print(head(combined_exp_corrected))
# ##########################################################################################################箱视图
#############################################################################################CFTR单独
color = c("#35A1D3","#BE4E4D")
my_comparisons <- list(c("CFTR_control", "CFTR_tumor"))
markerGenes  <- c('CFTR')
marker_combined_exp_corrected<-combined_exp_corrected[markerGenes[markerGenes %in% rownames(combined_exp_corrected)],]
df<-as.data.frame(marker_combined_exp_corrected)
colnames(df)<-c("CFTR")
df$type<-status
df_melt<-reshape::melt(df,id.vars=c("type"))
df_melt$value<-log2(df_melt$value)
df_melt$combine_name<-factor(paste(df_melt$variable,df_melt$type,sep="_"),levels = c("CFTR_control", "CFTR_tumor"))
df_melt$type<-factor(df_melt$type,levels = c("control", "tumor"))
# p<- df_melt %>% ggplot(aes(x=combine_name,y=value,fill=type))+
#   geom_boxplot(width=0.5,position=position_dodge(width=1),outlier.colour = NA)+
#   theme_bw()+
#   theme(legend.position="none",
#         axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size=15), # 调整xlabel的倾斜角度,
#         axis.text.y = element_text(hjust = 0.5, vjust = 0.5,size=15),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text( size=rel(1)),
#         panel.border = element_blank(),
#         axis.line = element_line(colour = "black",size=1)
#   )+scale_fill_manual(values = color)
# p1<-p +stat_compare_means(comparisons =my_comparisons, label.y =5.8,method="wilcox.test",label= "p.signif" ,size = 3.5)#label= "p.signif"

p=ggboxplot(df_melt, x="combine_name", y="value", color = "type", 
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
                     label.y = c(5.8),
                     ,label= "p.signif",
                     method="wilcox.test")+
  ylim(0.001,6.2) # Add horizontal line at base mean
p
library(dplyr)

quartiles_df <- df_melt %>%
  dplyr::group_by(combine_name, type) %>%
  dplyr::summarise(
    Q1 = quantile(value, 0.25, na.rm = TRUE),
    Median = quantile(value, 0.50, na.rm = TRUE),
    Q3 = quantile(value, 0.75, na.rm = TRUE)
  )

# A tibble: 2 × 5
# Groups:   combine_name [2]
# combine_name type       Q1 Median    Q3
# <fct>        <fct>   <dbl>  <dbl> <dbl>
# CFTR_control control  4.72   5.02  5.10
# CFTR_tumor   tumor    4.35   4.63  4.93
table(df_melt$type)
# control   tumor   
# 270     499 
#p<2.22e-16
###################################################################################################over
###############################################################################################################################
###########################################Figure 3B ###########################################################################
###############################################################################################################################
setwd("/400T/wangmeiheng/TCGA/UCSC")
load("bulk_combined_exp_corrected.RData")
library(ggplot2)
library(ggpubr)

color = c("#35A1D3","#BE4E4D")


##########################################################################CFTR
detach("package:plyr", unload=TRUE)
my_comparisons <- list(c("CFTR_control", "CFTR_tumor"))
markerGenes  <- c('CFTR')
marker_combined_exp_corrected<-combined_exp_corrected[markerGenes,]
marker_combined_exp_corrected<-na.omit(marker_combined_exp_corrected)
df<-as.data.frame(t(marker_combined_exp_corrected))
colnames(df)<-'CFTR'
df$type<-status
df_melt<-reshape::melt(df,id.vars=c("type"))
df_melt$combine_name<-factor(paste(df_melt$variable,df_melt$type,sep="_"))
summary_df<-df_melt %>%group_by(variable,type) %>%summarize(n = n())
print(summary_df, n = Inf)

levels(df_melt$combine_name)
df_melt$combine_name<-factor(paste(df_melt$variable,df_melt$type,sep="_"),levels = c("CFTR_control", "CFTR_tumor"))
df_melt$type<-factor(df_melt$type,levels = c("control", "tumor"))
p=ggboxplot(df_melt, x="combine_name", y="value", color = "type", 
            ylab="Gene expression",
            xlab="",
            legend.title='type',
            palette =color,
            width=0.4, add = "",
            bxp.errorbar=T,#显示误差条
            bxp.errorbar.width=0.1, #误差条大小
            size=0.5, #箱型图边线的粗细
            outlier.shape=NA)+#不显示outlier+#rotate_x_text(20)+#jitter
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(8.5,9.5,8.5),
                     method="wilcox.test",
                     label="p.signif")+
  ylim(-10,11) # Add horizontal line at base mean
p
quartiles_df <- df_melt %>%
  group_by(combine_name, type) %>%
  summarise(
    Q1 = quantile(value, 0.25, na.rm = TRUE),
    Median = quantile(value, 0.50, na.rm = TRUE),
    Q3 = quantile(value, 0.75, na.rm = TRUE)
  )
# A tibble: 2 × 5
#combine_name type       Q1 Median    Q3
#CFTR_control control  5.16   5.49  5.96
#CFTR_tumor   tumor    3.53   4.93  5.89
table(df_melt$type)
# control   tumor    
# 188        204
#p=3.3e-07


###############################################################################################################################
###########################################Figure 3C ###########################################################################
###############################################################################################################################
my_comparisons=list(c("Normal", "Tumor"))
#蛋白质组数据展示绘制
load('/400T/wangmeiheng/data_result/protein/PDC_all.Rdata')

###总数据--401 Samples
table(df$tumor_stage)

ggboxplot(df, x="tumor_stage", y="CFTR", width = 0.6, 
          color = "black",#轮廓颜色
          fill="tumor_stage",#填充
          palette = "npg",
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=1, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right")+
  stat_compare_means(method = "kruskal.test",
                     label = "p.format",
                     label.x = 0.8, 
                     label.y = 2.3,
                     show.legend = F) +
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(2.5),
                     method="wilcox.test",
                     label="p.signif"
  )

##PDC00393
ggboxplot(df_393, x="tumor_stage", y="CFTR", width = 0.6, 
          color = "black",#轮廓颜色
          fill="tumor_stage",#填充
          palette = "npg",
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=1, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right")+
  stat_compare_means(method = "kruskal.test",
                     label = "p.format",
                     label.x = 0.8, 
                     label.y = 2.3,
                     show.legend = F) +
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(2.5),
                     method="wilcox.test",
                     label="p.signif"
  )

ggboxplot(df_248, x="tumor_stage", y="CFTR", width = 0.6, 
          color = "black",#轮廓颜色
          fill="tumor_stage",#填充
          palette = "npg",
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=1, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right")+
  stat_compare_means(method = "kruskal.test",
                     label = "p.format",
                     label.x = 0.8, 
                     label.y = 2.3,
                     show.legend = F) +
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(2.5),
                     method="wilcox.test",
                     label="p.signif"
  )

