
library(tidyverse)
library(data.table)
library(ggthemes)
library(ggpubr)
library(openxlsx)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "4.Biology/"

# 读入estimate res
ESTIMATE_res <- read.csv("./TME/ESTIMATE_res.TCGA-COAD.csv")

group <- read.csv("./res/3.Signature/table/Risk_group.csv")
colnames(group)
group <- group[,c("id","riskScore")]
colnames(group)[c(1,2)] <- c("sample","risk_score")
data_pre <- ESTIMATE_res
colnames(data_pre)[1] <- "sample"
data_pre$sample <- substring(data_pre$sample,1,16)
data_merge <- merge(group, data_pre, by="sample")
data_merge <- data_merge[,-6]

scatter_ggplot_list <- list()
for(target_col in colnames(data_pre)[-c(1,5)]){
  scatter_ggplot_list[[target_col]] <- scatter_ggplot(inputArr=data_merge, 
                                                      outfile=NULL, title = NULL, 
                                                      xcol="risk_score", ycol=target_col, 
                                                      color="#4b84b3",
                                                      color_manual=NULL, 
                                                      col_pal_name="lancet",
                                                      facet.by=NULL, Marginal_Density=FALSE)
}
outfile_tmp <- paste0(res_home,proj_name, "/pdf/","Cor.","risk_TME", ".pdf")
multiplot_ggarrange(ggobj_list=scatter_ggplot_list, outfile=outfile_tmp, labels="AUTO", ncol=3, nrow=1, 
                    legend="bottom", common.legend=TRUE, width=12, height=4)

