
library(data.table)
library(tidyverse)
library(maftools)
library(openxlsx)
library(survival)
library(survminer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "4.Biology/"

# load data
load("E:/BaiduSyncdisk/r_project/data/TCGA_data/SNV/TCGA-COAD.SNV.Rdata")
head(snv[,1:4])
snv <- snv[,-1]

group <- read.csv("./res/3.Signature/table/Risk_group.csv")
colnames(group)
group <- group[,c("id","risk")]
colnames(group)[c(1,2)] <- c("id","group")
group_sub <- group[!duplicated(group$id),]
group_sub$id <- substring(group_sub$id,1,12)
group_sub <- data.frame(group_sub,row.names = group[,1])

maf_data <- read.maf(maf = snv)
maf.sample <- as.data.frame(maf_data@clinical.data)
maf.sample$id <- substring(maf.sample$Tumor_Sample_Barcode,1,12)
maf.sample$id <- gsub("-",".",maf.sample$id)
maf.sample <- merge(maf.sample,group_sub,"id")

pres_idx1 <- maf.sample$Tumor_Sample_Barcode[maf.sample$group %in% "high"]
pres_idx2 <- maf.sample$Tumor_Sample_Barcode[maf.sample$group %in% "low"]
maf_c1 <-  subsetMaf(maf = maf_data, tsb = pres_idx1)
maf_c2 <- subsetMaf(maf = maf_data, tsb = pres_idx2) 
risk_color <- list(group = c(high = "#f89f68",low = "#4b84b3"))

# 运行时切换
outfile_tmp <- paste0(res_home,proj_name,"pdf/","SNV.","risk_","high",".pdf");outfile_tmp
pdf(file = outfile_tmp, onefile = F, height = 6, width = 7)
oncoplot(maf = maf_c1, bgCol = "#F4F4F4")
dev.off()