########################   GSVA   ######################## 
#### 

library(data.table)
library(tidyverse)
library(stringr)
library(openxlsx)
library(GSVA)
library(ggpubr)
library(enrichplot)
library(clusterProfiler)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "1.degs/"

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

# interest gene
interest_genes <- read.xlsx("./list/lactylation-related genes.xlsx")

##############################  富集分析 -- gsea  ##############################
# 转为矩阵
data_pre <- as.matrix(data_exprs)
# ssGSEA
genelist <- list(Lactylation = interest_genes$Gene)
data_ES <- GSVA::gsva(data_pre, genelist, method="gsva")
data_ES <- t(data_ES) %>% data.frame()

# 计算zscore
data_ES_zscore <- scale(data_ES) %>% data.frame()
## 保存数据
outfile_tmp <- paste0(res_home,proj_name,"/RData/","GSVA.","Lacty",".RData");outfile_tmp
save(data_ES,data_ES_zscore,file = outfile_tmp)

##############################  degs  ##############################
data_pre <- rownames_to_column(data_ES_zscore,"sample")
### 区分tumor和normal
table(substring(data_pre$sample,first = 13,last = 16))#观察样本类型
tumor_idx <- grep("01A",substring(data_pre$sample,first = 13,last = 16))
normol_idx <- grep("11A",substring(data_pre$sample,first = 13,last = 16))
data_pre$group <- "normal"
data_pre$group[tumor_idx] <- "tumor"
table(data_pre$group)

# 计算p值 ANOVA
# ref: https://www.omicsclass.com/article/743
# ref: https://www.omicsclass.com/article/733
shapiro.test(data_pre$Lactylation)
pval <- wilcox.test(Lactylation~group,data_pre)
pval$p.value
p = signif(pval$p.value,2);p 

## 绘图
str(data_pre)
data_plot <- data_pre
data_plot$group <- factor(data_plot$group,levels = c("tumor","normal"))

# color
my_palet <- c("#f89f68","#4b84b3")

# save
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","Boxplot.",
                      "tumor_normal.", "Lactylation",".pdf");outfile_tmp
pdf(outfile_tmp,width = 3,height = 4)
ggplot(data_plot, 
       aes(x = group, y = Lactylation)) +
  geom_boxplot(aes(fill = group), show.legend = T, width = 0.6) +
  scale_fill_manual(values = c('#f89f68','#4b84b3')) + # 箱线图的填充色
  scale_color_manual(values = c('#f89f68','#4b84b3')) + # 点的填充色
  geom_point(aes(color = group),size = 1,alpha=1) + #绘制样本点
  stat_compare_means(method = "wilcox.test",label = "p.format") +
  ##以下是ggplot2的主题设置，修改边框、背景、标题、字体等
  theme_bw()+ theme_classic() +
  ylab('Lactylation') +
  xlab('') +
  ggtitle('') + 
  theme(plot.title = element_text(color = 'black', size = 15, hjust = 0.5)) # 标题居中。
  # labs(subtitle = paste0("p-value = ", p))
dev.off()


##############################  degs  -- pair ##############################
data_pre <- rownames_to_column(data_ES_zscore,"sample")
### 区分tumor和normal
table(substring(data_pre$sample,first = 13,last = 16))#观察样本类型
data_pre$group <- "normal"
data_pre$group[str_sub(data_pre$sample,14,15) <= 10] <- "tumor"
table(data_pre$group)

# tumor
data_tumor <- data_pre[str_sub(data_pre$sample,14,15) <= 10,]
# 有一些样本会重复 delete
pres_idx = !duplicated(str_sub(data_tumor$sample,1,12));table(pres_idx)
data_tumor = data_tumor[pres_idx,]

# normal
data_normal <- data_pre[str_sub(data_pre$sample,14,15) > 10,]
pres_idx = !duplicated(str_sub(data_normal$sample,1,12));table(pres_idx)
data_normal = data_normal[pres_idx,]

# 根据正常样本 提取肿瘤样本
patient = str_sub(data_normal$sample,1,12)
pres_idx = str_sub(data_tumor$sample,1,12) %in% patient;table(pres_idx) 
data_tumor = data_tumor[pres_idx,]

# 保证配对顺序
data_normal <- data_normal[match(str_sub(data_tumor$sample,1,12),
                                 str_sub(data_normal$sample,1,12)),]
# merge
data_merge <- rbind(data_tumor, data_normal)
# add group2
data_merge$group2 <- str_sub(data_merge$sample,1,12)

###########################     配对箱线图    ###########################  
data_pre <- data_merge
# orde
table(data_pre$group)
data_pre$group <- factor(data_pre$group,levels = c("tumor","normal"))

# 去重
data_pre$sample <- substring(data_pre$sample,1,15)
data_pre <- data_pre[!duplicated(data_pre$sample),]

# t检验
t <- data_pre[which(data_pre$group == "tumor"),]
n <- data_pre[which(data_pre$group == "normal"),]
p.val.merge <- data.frame(p.val = 1,row.names = "Lactylation")
for (gene in "Lactylation"){
  p.val <-  t.test(t[,gene],n[,gene],paired = T)$p.value
  p.val = signif(p.val,3)
  print(p.val)
  p.val.merge[gene,] = p.val
}

#### plot 
data_plot <- data_pre
data_plot$group <- factor(data_plot$group,levels = c("tumor","normal"))

# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Boxplot.",
                      "Lactylation",".paired",".pdf");outfile_tmp
pdf(file = outfile_tmp,width = 3,height = 4,onefile = F)
ggplot(data_plot, 
       aes(x = group, y = Lactylation)) +
  geom_boxplot(aes(fill = group), show.legend = T, width = 0.6) +
  scale_fill_manual(values = c('#f89f68','#4b84b3')) + # 箱线图的填充色
  scale_color_manual(values = c('#f89f68','#4b84b3')) + # 点的填充色
  geom_point(aes(color = group),size = 1,alpha=1) + #绘制样本点
  geom_line(aes(group = group2), alpha=0.5, color = 'gray',lwd = 0.5) + #绘制配对样本间连线
  stat_compare_means(method = "wilcox.test",label = "p.format",paired = T) +  # 统计学检验
  ##以下是ggplot2的主题设置，修改边框、背景、标题、字体等
  theme_bw()+ theme_classic() +
  ylab('Lactylation') +
  xlab('') +
  ggtitle('') + 
  theme(plot.title = element_text(color = 'black', size = 15, hjust = 0.5)) # 标题居中。
# labs(subtitle = paste0("p.value = ", p.val))
dev.off()

