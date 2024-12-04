
library(tidyverse)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(rstatix)
library(stringr)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "5.Drug/"

# load data
data_drug <- read.csv("./calcPhenotype_Output/DrugPredictions.csv")
head(data_drug[,1:4])
colnames(data_drug)[1] <- "sample"

group <- read.csv("./res/3.signature/table/Risk_group.csv")
colnames(group)[1] <- "sample"
group_sub <- group[,c("sample","riskScore","risk")]
data_pre <- data_drug
data_pre$sample <- substring(data_pre$sample,1,16)
data_merge <- merge(group_sub, data_pre, by = "sample")
data_merge <- data_merge[!duplicated(data_merge$sample),]
data_merge <- data.frame(data_merge[,-c(1,3)], row.names = data_merge[,1])
data_pre <- data_merge
colnames(data_pre)[1] <- "RiskScore"
colnames(data_pre) <- str_split(colnames(data_pre), "_",simplify = T)[,1]

drug <- colnames(data_pre)[-1]
# 计算
res_cor <- cor_asymmetry(input1=data_pre[, "RiskScore",drop = F],
                         input2=data_pre[, drug],
                         cor_method="spearman")
res_cor <- lapply(res_cor, function(x){
  res <- data.frame(geneset=rownames(x), x)
  res <- convArrType(inputArr=res, keyCols="geneset", valueCols=colnames(res)[-1], wide2long=TRUE)
  return(res)
})
cor_res <- merge(res_cor$correlation, res_cor$pvalue, by=c("geneset", "variable"))
colnames(cor_res) <- c("var1", "var2", "corr", "pval")
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Cor.","risk_","drug",".csv");outfile_tmp
write.csv(cor_res,outfile_tmp,row.names = F)
