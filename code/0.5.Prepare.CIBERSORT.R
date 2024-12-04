########################  CIBERSORT   ######################## 

library(data.table)
library(tidyverse)
library(devtools)
library(e1071)
library(preprocessCore)
library(parallel)
library(CIBERSORT)
# devtools::install_github("Moonerss/CIBERSORT")

rm(list = ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

#### CIBERSORT
# cibersort
results <- cibersort(sig_matrix = LM22, mixture_file = data_exprs,
                     perm = 1000,QN = F) # 芯片 T, 测序 F

# save
write.csv(results,file = "./TME/CIBERSORT_res.TCGA-COAD.csv")

