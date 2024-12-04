########################     ESTIMATE     ######################## 
# install.packages("estimate", repos="http://R-Forge.R-project.org")

library(estimate)
library(tidyverse)
library(data.table)

rm(list = ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "4.Biology/"

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

##############################  ESTIMATE  ##############################
### 保存exprs data为 txt文件
outfile_tmp <- paste0(res_home,proj_name, "/table/","TCGA-COAD.for_ESTIMATE", ".txt");outfile_tmp
write.table(data_exprs, file = outfile_tmp, sep = '\t',quote = F)

### 表达矩阵与作者gene dataset 交集
input_file <- paste0(res_home,proj_name, "/table/","TCGA-COAD.for_ESTIMATE", ".txt")
outfile_tmp <- paste0(res_home,proj_name, "/table/","TCGA-COAD", ".estimate_gene.gct")
filterCommonGenes(input.f = input_file,
                  output.f = outfile_tmp ,
                  id="GeneSymbol")
### estimate Score
input_file <- paste0(res_home,proj_name, "/table/","TCGA-COAD", ".estimate_gene.gct")
outfile_tmp <- paste0(res_home,proj_name, "/table/","TCGA-COAD",".ESTIMATE_score", ".gct")
estimateScore(input.ds = input_file,
              output.ds = outfile_tmp,
              platform = "illumina") ## 注意platform

##############################  读取ESTIMATE  ##############################
outfile_tmp <- paste0(res_home,proj_name, "/table/","TCGA-COAD",".ESTIMATE_score", ".gct")
estimate_res <- read.table(outfile_tmp, skip = 2, header = T)
rownames(estimate_res) <- estimate_res[,1]
estimate_res <- t(estimate_res[,3:ncol(estimate_res)])
# 计算肿瘤纯度
TumorPurity = cos(0.6049872018+0.0001467884 * estimate_res[,3])
# 合并
ESTIMATE_scores <- cbind(estimate_res,TumorPurity)
ESTIMATE_scores <- as.data.frame(ESTIMATE_scores)

# save
write.csv(ESTIMATE_scores,file = "./TME/ESTIMATE_res.TCGA-COAD.csv")


