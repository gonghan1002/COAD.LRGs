
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

data_drug <- read.csv("./calcPhenotype_Output/DrugPredictions.csv")
head(data_drug[,1:4])
colnames(data_drug)[1] <- "sample"
group <- read.csv("./res/3.signature/table/Risk_group.csv")
colnames(group)[1] <- "sample"
group_sub <- group[,c("sample","riskScore","risk")]
interest_drug <- c("Dihydrorotenone","TAF1","Ulixertinib","Lapatinib","VX.11e",
                   "BMS.754807","NU7441","AZD8055","JQ1","AZD8186")
data_pre <- data_drug
data_pre$sample <- substring(data_pre$sample,1,16)
data_merge <- merge(group_sub, data_pre, by = "sample")
data_merge <- data_merge[!duplicated(data_merge$sample),]
data_merge <- data.frame(data_merge[,-c(1,3)], row.names = data_merge[,1])

data_pre <- data_merge
colnames(data_pre)[1] <- "RiskScore"
colnames(data_pre) <- str_split(colnames(data_pre), "_",simplify = T)[,1]
data_pre <- data_pre[,c("RiskScore",interest_drug)]
data_pre <- as.data.frame(scale(data_pre))
# calculate
cor_res_list <- list()
for(drug in interest_drug){
  cor_res_list[[drug]] <- scatter_ggplot(inputArr=data_pre, outfile=NULL,
                                         xcol="RiskScore", ycol= drug,
                                         title = NULL, color="#4b84b3",
                                         color_manual=NULL, col_pal_name="lancet",
                                         facet.by=NULL, Marginal_Density=FALSE)
}

outfile_tmp <- paste0(res_home, proj_name,"/pdf/","Cor.","risk_drug" , ".pdf")
multiplot_ggarrange(ggobj_list=cor_res_list, outfile=outfile_tmp, labels=NULL, 
                    ncol= 5, nrow = 2,
                    legend="bottom", common.legend=TRUE, width=16, height=8)
