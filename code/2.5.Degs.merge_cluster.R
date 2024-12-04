########################  差异基因分析  ######################## 
#### 

library(data.table)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrastr)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "2.Cluster/"

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.Tumor.RDS")
head(data_exprs[,1:4])

# 读取分组数据
group <- read.csv("./res/2.Cluster/table/Cluster.Consensus.k3.csv")
# 合并cluster2 3
colnames(group) <- c("sample","group")
group$group <- paste0("Cluster", group$group)
table(group$group)

group$group[group$group %in% c("Cluster2","Cluster3")] <- "Cluster2"
table(group$group)

###########################  degs  --limma  ########################### 
### 设置注释矩阵
data_anno <- data.frame(group = group[,c("group")],
                        row.names = group$sample)
# match
data_pre <- data_exprs[,match(rownames(data_anno),colnames(data_exprs))]

# # 过滤
# pres = apply(data_pre>1,1,sum) 
# to_keep = pres > 0.5 * ncol(data_pre)
# table(to_keep) # 看一下有多少满足条件
# data_pre = data_pre[to_keep,]

# degs
degs_res <- diff_limma(inputData = data_pre, colData = data_anno,
                       contrast_name = c("Cluster1", "Cluster2"),
                       log2trans = T, adjust.method = "BH", 
                       pval_cutoff = NULL,log2fc_cutoff = NULL)
# 加标签列
degs_res$threshold = factor(ifelse(degs_res$adj.P.Val < 0.05 & abs(degs_res$logFC) >= 0.585, 
                                   ifelse(degs_res$logFC >= 0.585,'Up','Down'),'NoSig'),
                            levels=c('Up','Down','NoSig'))
table(degs_res$threshold != "NoSig" )
# save
outfile_tmp <- paste0(res_home, proj_name,"/table/","degs.",
                      "Cluster1",".vs.","merge_Cluster2",".merge",".fc1.5", ".csv");outfile_tmp
write.csv(degs_res, file = outfile_tmp)

#################     火山图  volcano     #####################
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","degs.",
                      "Cluster1",".vs.","merge_Cluster2",".fc1.5.2", ".pdf")
pdf(file = outfile_tmp,onefile = F,width = 6,height = 6)
ggplot(degs_res,aes(x = logFC, y = -log10(adj.P.Val)))+
  geom_point_rast(size = 2, aes(color = threshold))+
  scale_color_manual(values = c("#f89f68","#4b84b3","gray50"),
                     labels = c("Up","Down","NoSig"))+ #确定点的颜色
  geom_vline(xintercept=c(-0.585,0.585),linetype ="longdash",col="black",lwd=0.5) + 
  geom_hline(yintercept = -log10(0.05),linetype ="longdash",col="black",lwd=0.5)+ 
  theme_bw()+ #修改图片背景
  theme(
    legend.title = element_blank(), # 不显示图例标题
    panel.grid=element_blank(),
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.key.size = unit(24, "pt")
  )
dev.off()

