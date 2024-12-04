########################  差异基因分析  ######################## 
#### 

library(data.table)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ggthemes)

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

# read group
group <- read.csv("./res/2.Cluster/table/Cluster.Consensus.k3.csv")
colnames(group) <- c("sample","group")
group$group <- paste0("Cluster", group$group)

### 选取两两分组
group_sub <- group[group$group %in% c("Cluster1","Cluster2"),]

###########################  degs  --limma  ########################### 
### 设置注释矩阵
data_anno <- data.frame(group = group_sub[,c("group")],
                        row.names = group_sub$sample)
# match
data_pre <- data_exprs[,match(rownames(data_anno),colnames(data_exprs))]

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
                      "Cluster1","_vs_","Cluster2",".fc1.5", ".csv");outfile_tmp
write.csv(degs_res, file = outfile_tmp)

#################     火山图  volcano     #####################
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","degs.",
                      "Cluster1","_vs_","Cluster2",".fc1.5", ".pdf")
volcano_res <- volcano_plot(input=degs_res, output=outfile_tmp,
                            log2fc_cutoff=0.585, pval_cutoff=0.05,
                            title=NULL,padj_cutoff=0.05, genes_highlight=NULL, 
                            color_map=c(Down="#4b84b3", Not="grey", Up="#f89f68"), 
                            width=7, height=7, dpi=300, device="pdf")

