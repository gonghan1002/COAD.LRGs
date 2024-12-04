
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "4.Biology/"

# 读取TIDE数据
data_exprs <- read.csv("./TME/TIDE_res.TCGA-COAD.csv")

group <- read.csv("./res/3.signature/table/Risk_group.csv")
colnames(group)
group_sub <- group[,c("id","risk","riskScore")]
colnames(group_sub)[c(1,2,3)] <- c("id","group","risk_score")
colnames(data_exprs)
colnames(data_exprs)[1] <- "id"
data_expr_sub <- data_exprs[,c(1,3,4)]
data_expr_sub$id <- substring(data_expr_sub$id,1,16)
data_merge <- merge(data_expr_sub, group_sub, by = "id")

data_pre <- data_merge
table(data_pre$Responder)
data_pre$Responder <- ifelse(data_pre$Responder=="True","R","NR") 
cutpoint <- quantile(data_pre[ ,"risk_score"], probs=0.5);cutpoint
data_pre$group <- ifelse(data_pre$risk_score < cutpoint,"low","high")
data_plot <- data_pre[,c("Responder","TIDE","risk_score","group")]
chisq.test(data_plot$group,data_plot$Responder)
 
library(ggpubr)
data_plot$score <- scale(data_plot$risk_score)
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Barplot.","TIDE.risk",".pdf");outfile_tmp
pdf(file=outfile_tmp,onefile = FALSE, height = 4, width = 3)   
ggplot(data = data_plot,
       aes(y = TIDE,
           x = group))+
  geom_boxplot(alpha = 1,
               fill = c("#f89f68", "#4b84b3"))+
  stat_compare_means(method = "t.test",label = "p.format") + 
  theme_bw()+ theme_classic() +
  ylab('TIDE') +
  xlab('group')
dev.off()
