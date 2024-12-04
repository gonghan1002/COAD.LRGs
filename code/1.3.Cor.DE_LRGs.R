########################   验证 --   GSE157011  ######################## 
#### median fail
#### maxstat 
#### ROC fail 

library(data.table)
library(stringr)
library(tidyverse)
library(openxlsx)
library(igraph)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "1.degs/"

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.Tumor.RDS")
head(data_exprs[,1:4])
range(data_exprs)

# 读取显著基因
interest_genes <- read.csv("./res/1.degs/table/degs.LRGs.fc1.csv")
interest_genes <- interest_genes[interest_genes$threshold != "NoSig",][,1]
interest_genes

### data handle
data_sub <- data_exprs[interest_genes,]

#########################   cor   #########################
# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Cor.","DE_LRGs",".pdf");outfile_tmp
pdf(file = outfile_tmp,width = 6,height = 6)
t(data_sub) %>% corrr::correlate(method = "spearman") %>%
  corrr::network_plot(colours = c("#4b84b3","white", "#f89f68"),
                      repel = T,min_cor = .3)
dev.off()


