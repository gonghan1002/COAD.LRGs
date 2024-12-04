########################  snv  ######################## 

library(data.table)
library(tidyverse)
library(maftools)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "1.degs/"

# load data
load("E:/BaiduSyncdisk/r_project/data/TCGA_data/SNV/TCGA-COAD.SNV.Rdata")
head(snv[,1:4])
snv <- snv[,-1]

# read genes
interested_genes <- read.csv("./res/1.degs/table/degs.LRGs.fc1.csv")
interested_genes <- interested_genes[interested_genes$threshold != "NoSig",]
interested_genes <- interested_genes[,1];interested_genes

###############  读入maf文件   ########################
# 读取maf文件
maf_data <- read.maf(maf = snv)

#### plot 
outfile_tmp <- paste0(res_home,proj_name,"pdf/","SNV.","LRGs",".pdf");outfile_tmp
# save
pdf(file = outfile_tmp, onefile = F, height = 7, width = 6)
oncoplot(maf = maf_data,genes = interested_genes,
         fontSize = 0.4,
         draw_titv = T,bgCol = "#F4F4F4"
)
dev.off()

