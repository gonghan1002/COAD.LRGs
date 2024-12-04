########################   clinical feature   ######################## 

library(data.table)
library(tidyverse)
library(survival)
library(openxlsx)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "2.Prognosis/"

# read data
data_clinical <- fread("./data/TCGA-COAD.GDC_phenotype.tsv.gz",data.table = F)

### extract
colnames(data_clinical)
table(data_clinical$sample_type.samples)
data_clinical_sub <- data_clinical[,c(1,2,73,97,56,28,29,93,41,40,39)]

### 处理临床数据
meta <- data_clinical_sub
colnames(meta)
colnames(meta) <- c("sample","age","gender","BMI","venous_invasion","lymphatic_invasion",
                    "MSI","stage","stage_T","stage_N","stage_M")
meta$sample <- gsub("-",".",meta$sample)
str(meta)

### venous invasion
idx_tmp <- meta$venous_invasion
table(idx_tmp)
idx_tmp[idx_tmp == ""] = NA
table(idx_tmp)
meta$venous_invasion <- idx_tmp

### lymphatic invasion
idx_tmp <- meta$lymphatic_invasion
table(idx_tmp)
idx_tmp[idx_tmp == ""] = NA
table(idx_tmp)
meta$lymphatic_invasion <- idx_tmp

### MSI
idx_tmp <- meta$MSI
table(idx_tmp)
idx_tmp[idx_tmp == ""] = NA
table(idx_tmp)
meta$MSI <- idx_tmp

### stage
idx_tmp <- meta$stage
table(idx_tmp)
idx_tmp <- gsub("stage iv[abc]","Stage IV",
                gsub("stage iii[abc]", "Stage III", 
                     gsub("stage ii[abc]", "Stage II",
                          gsub("stage i[abc]", "Stage I", idx_tmp))))
idx_tmp[idx_tmp == "not reported"] = NA
idx_tmp[idx_tmp == ""] = NA
idx_tmp[idx_tmp == "stage iv"] = "Stage IV"
idx_tmp[idx_tmp == "stage iii"] = "Stage III"
idx_tmp[idx_tmp == "stage ii"] = "Stage II"
idx_tmp[idx_tmp == "stage i"] = "Stage I"
table(idx_tmp)
meta$stage <- idx_tmp

### T
idx_tmp <- meta$stage_T
table(idx_tmp)
idx_tmp <- gsub("T4[a-z]", "T4", idx_tmp)
idx_tmp <- gsub("Tis", NA, idx_tmp)
idx_tmp[idx_tmp == ""] = NA
table(idx_tmp)
meta$stage_T <- idx_tmp

### N
idx_tmp <- meta$stage_N
table(idx_tmp)
idx_tmp <- gsub("N1[a-z]", "N1", idx_tmp)
idx_tmp <- gsub("N2[a-z]", "N2", idx_tmp)
idx_tmp <- gsub("NX", NA, idx_tmp)
idx_tmp[idx_tmp == ""] = NA
table(idx_tmp)
meta$stage_N <- idx_tmp

### M
idx_tmp <- meta$stage_M
table(idx_tmp)
idx_tmp <- gsub("M1[a-z]", "M1", idx_tmp)
idx_tmp <- gsub("MX", NA, idx_tmp)
idx_tmp[idx_tmp == ""] = NA
table(idx_tmp)
meta$stage_M <- idx_tmp

# save
write.csv(meta,"./data/TCGA-COAD.clinical_info.csv",row.names = F)


