library(pheatmap)
library(stringr)
library(reshape2)
###############################################################################################################################
###########################################Figure 6D###########################################################################
###############################################################################################################################
load("/400T/wangmeiheng/FJ/bulk/house/mouse_data.RData")
load("/400T/wangmeiheng/FJ/bulk/house/mouse_difgene.RData")
load("/400T/wangmeiheng/FJ/bulk/house/df.RData")
mat<-mouse_data

gene_list <- unlist(strsplit(df$geneID, "/"))
df_x<-data.frame(df$Description,df$geneID)

col_df <- data.frame(
  Group_col = c(rep("shRNA",6), rep("control",6))
)
rownames(col_df) <- colnames(mat)

row_df <- data.frame(
  Group = unlist(
    lapply(1:nrow(df_x), function(i) {
      rep(df_x[i,1], str_count(df_x[i,2], "/") + 1)
    })
  )
)
rownames(row_df) <- rownames(mat[gene_list ,])
library(RColorBrewer)
pathways <- unique(row_df$Group)

pathway_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(length(pathways)),
  pathways
)


col_colors <- c(
  shRNA = "#1f77b4",  
  control = "#ff7f0e"  
)

annotation_colors <- list(
  Group = pathway_colors,   
  Group_col = col_colors   
)
row_labels<- sub("\\..*$", "", rownames(mat[gene_list ,]))

p<-pheatmap(mat[gene_list ,],
            color =colorRampPalette(c("#113F8C","white","#AF1E23"))(50),
            scale = "row",         
            cluster_rows = F,    
            cluster_cols = F,   
            show_rownames = TRUE,   
            show_colnames = F,   
            fontsize = 10,          
            annotation_col =col_df,   
            annotation_row = row_df, 
            labels_row = row_labels,
            annotation_colors = annotation_colors,
            main = ""       
)

#################foldchange value##########################
library(tidyverse)
library(ggthemes)
library(ggprism)
diff <- diff_signif
rownames(diff)<-diff$SYMBOL
diff<-diff[gene_list,]
diff <- rbind(
  subset(diff, logFC > 0),
  subset(diff, logFC < 0)
)

# Load differential analysis results
dat_plot <- data.frame(id=row.names(diff), p=diff$P.Value,logFC= diff$logFC)
dat_plot$group<- ifelse(dat_plot$logFC>0 ,1,-1) 


dat_plot$g<-"Not"
dat_plot$g[ dat_plot$logFC>0 & dat_plot$p< 0.05 ] <-"Up"
dat_plot$g[ dat_plot$logFC<0 & dat_plot$p< 0.05 ] <-"Down"
table(dat_plot$g)
dat_plot$g<- factor(dat_plot$g, levels=c('Up','Down','Not'))

# Add label colour
dat_plot$color<- ifelse(dat_plot$g=="Not",'#cccccc',"black")

#dat_plot <- dat_plot[order(dat_plot$logFC,decreasing = T), ]
dat_plot$id<- factor(dat_plot$id,levels = rev(dat_plot$id))
dat_plot$lable_hjust<- ifelse(dat_plot$logFC>0, 1, 0)
dat_plot$lable_xloc<- ifelse(dat_plot$logFC>0, -0.05, 0.05)
fc_up <- min(dat_plot[dat_plot$g=="Up","logFC"]);
if(fc_up==Inf){fc_up=0}
fc_down <- min(dat_plot[dat_plot$g=="Down","logFC"]);


p1 <- ggplot(data = dat_plot,aes(x = id, y =logFC, fill = g)) +
  geom_col() +
  scale_y_continuous("",breaks = seq(-6,0),   expand = expansion(mult = c(2, 2)),
                     limits = c(fc_down, fc_up),position = "right")+
  scale_x_discrete("")+
  scale_fill_manual(values = c('Up'='#36648b','Not'='#cccccc','Down'='#7ccd7c')) +
  geom_hline(yintercept = c(fc_down,fc_up), color ='white', size = 0.5,lty='dashed') +
  ylab('t value of gsva score,ST11-K47 vs ST11-K64') +
  theme(
    panel.grid = element_blank(),
    legend.position = "none", 
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),  
    panel.background = element_rect(fill = "transparent", colour = NA),
    #axis.line.x = element_line(colour = "black")
    axis.ticks.length = unit(0, "cm")
  )+
  coord_flip()
p1
#####################FDR heatmap##############
mat1<-data.frame(diff$SYMBOL,diff$adj.P.Val)
rownames(mat1)<-rownames(diff)
p2<-pheatmap(as.matrix(mat1[,2]),
             color = colorRampPalette(c("#FF9900","#FFFFCC" ,"white" ))(50), 
             #scale = "row",         
             cluster_rows = F,   
             cluster_cols = F, 
             show_rownames = F,  
             show_colnames = F,  
             fontsize = 10,          
             #annotation_col =col_df,    
             #annotation_row = row_df, 
             labels_row = diff$SYMBOL,
             #annotation_colors = annotation_colors,
             main = ""     
)
p2




library(readxl)
library(survival)
library(survminer)
###############################################################################################################################
###########################################Figure 6H###########################################################################
###############################################################################################################################
df <- read_excel("/400T/wangmeiheng/TCGA/KM/house/house_os.xlsx")
unicox <-coxph(Surv(time = OS_time, event = OS_status) ~ group, data = df) 
unisum<- summary(unicox)   
pvalue <- round(unisum$coefficients[,5],3) 
hh<-cbind(
  HR=round(unisum$coefficients[,2],2),
  L95CI=round(unisum$conf.int[,3],2),
  H95CI=round(unisum$conf.int[,4],2),
  pvalue=round(unisum$coefficients[,5],3)
)
write.table(hh,"/400T/wangmeiheng/TCGA/KM/house/HR.txt",row.names = T,col.names = T)
###############################################################################################################################
###########################################Figure 6E ###########################################################################
###############################################################################################################################
#KM survival curve
fit<-survfit(Surv(time = OS_time, event = OS_status) ~ group,data=df)
survtest <- survdiff(Surv(time = OS_time, event = OS_status) ~ group,data=df)
ggsurvplot(fit, 
           data = df, 
           pval = TRUE, 
           pval.method.size = 3.5,         
           surv.median.line = "hv",         
           #conf.int = TRUE,                 
           #conf.int.style = "step",        
           conf.int.alpha = 0.1,
           risk.table.height = 0.1,        
           risk.table = TRUE, 
           risk.table.pos = 'in',          
           risk.table.title = "",          
           add.all = FALSE,                
           xlab = "Time (Days)", 
           xlim = c(0, max(df$OS_time)),
           
           palette = c("#3090a1", "grey", "#bc5148"),  
           legend.title = colnames(df)[4], 
           break.time.by = 20              
           
)
###############################################################################################################################
###########################################Figure 6F ###########################################################################
###############################################################################################################################
#KM survival curve
df_new1<-df[1:10,]
fit_new1<-survfit(Surv(time = OS_time, event = OS_status) ~ group,data=df_new1)
ggsurvplot(fit_new1, 
           data = df_new1, 
           pval = TRUE, 
           pval.method.size = 3.5,          
           surv.median.line = "hv",         
           #conf.int = TRUE,                
           #conf.int.style = "step",   
           conf.int.alpha = 0.1,
           risk.table.height = 0.1,        
           risk.table = TRUE, 
           risk.table.pos = 'in',          
           risk.table.title = "",          
           add.all = FALSE,                
           xlab = "Time (Days)", 
           xlim = c(0, max(df_new1$OS_time)),
           
           palette = c("#3090a1", "#bc5148"), 
           legend.title = colnames(df_new1)[4], 
           break.time.by = 20              
           
)
###############################################################################################################################
###########################################Figure 6G ###########################################################################
###############################################################################################################################
#KM survival curve
df_new2<-df[c(1:5,11:15),]
fit_new2<-survfit(Surv(time = OS_time, event = OS_status) ~ group,data=df_new2)
ggsurvplot(fit_new2, 
           data = df_new2, 
           pval = TRUE, 
           pval.method.size = 3.5,          
           surv.median.line = "hv",        
           #conf.int = TRUE,              
           #conf.int.style = "step",       
           conf.int.alpha = 0.1,
           risk.table.height = 0.1,         
           risk.table = TRUE, 
           risk.table.pos = 'in',        
           risk.table.title = "",         
           add.all = FALSE,                 
           xlab = "Time (Days)", 
           xlim = c(0, max(df_new2$OS_time)),
           
           palette = c("#3090a1", "#bc5148"),  
           legend.title = colnames(df_new2)[4], 
           break.time.by = 20              
           
)
