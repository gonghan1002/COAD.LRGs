
library(FactoMineR)
library(factoextra)
library(missMDA)
library(RColorBrewer)
library(data.table)
library(tidyverse)
library(openxlsx)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "3.Signature/"

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.Tumor.RDS")
head(data_exprs[,1:4])

# 读取分组数据
group <- read.csv("./res/3.Signature/table/Risk_group.csv")
colnames(group)[1] <- "sample"
group_sub <- group[,c(1,8)]
group_sub <- group_sub[!duplicated(group_sub$sample),]

interest_genes <- read.csv("./res/3.Signature/table/LASSO.gene_coef.csv")
interest_genes <- interest_genes[,1][-1]
interest_genes

data_pre <- data_exprs[interest_genes,]
colnames(data_pre) <- substring(colnames(data_pre),1,16)
group_sub <- group_sub[group_sub$sample %in% colnames(data_pre),]
data_pre <- data_pre[,group_sub$sample]

## PCA
pca_data <- as.data.frame(t(data_pre))
pca_res <- PCA(pca_data, graph = FALSE) #imputePCA()
# plot
my_palet <- c("#f89f68", "#4b84b3")
# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","PCA.","risk_group",".pdf");outfile_tmp
pdf(file=outfile_tmp, onefile = F,width = 5,height =5)
fviz_pca_ind(pca_res,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group_sub$risk, # color by groups
             palette = my_palet,
             addEllipses = F, # Concentration ellipses
             legend.title = "group",
             title = ""
)
dev.off()
