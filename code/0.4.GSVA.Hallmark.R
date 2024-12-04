########################   ######################## 
#### 

library(openxlsx)
library(tidyverse)
library(openxlsx)
library(GSVA)
library(data.table)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
# data_exprs <- log2(data_exprs + 1)
head(data_exprs[,1:4])

### ssGSEA
data_input <- data_exprs
gmtfile <- paste0()
gsva_res <- gsva_score(input_exprs=data_input,
                       gene_set_list=NULL,
                       gmtfile="E:/BaiduSyncdisk/r_project/data/Msigdb/h.all.v2022.1.Hs.symbols.gmt", 
                       method="gsva", kcdf="Gaussian", ssgsea.norm=TRUE, threads=8)
# save
write.csv(gsva_res, file = "./GSVA/GSVA.hallmark.csv")

