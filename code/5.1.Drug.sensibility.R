
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

checkDir(paste0(res_home,proj_name,"/pdf/"))
checkDir(paste0(res_home,proj_name,"/RData/"))
checkDir(paste0(res_home,proj_name,"/table/"))

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
data_merge <- data.frame(data_merge[,-c(1,2)], row.names = data_merge[,1])

res_test <- melt(data_merge,id.vars=c("risk"),variable.name = "drug") %>%
  group_by(drug) %>% t_test(value ~ risk) %>% 
  adjust_pvalue(method = "fdr") %>% add_significance("p.adj")
res_test$threshold = factor(ifelse(res_test$p.adj.signif < 0.05 & abs(res_test$statistic)>=0, 
                                   ifelse(res_test$statistic>=0,'Up','Down'),'NoSig'),
                            levels=c('Up','Down','NoSig'))
outfile_tmp <-  paste0(res_home,proj_name,"table/","Drug_auc.t_test",".csv");outfile_tmp
write.csv(res_test,file = outfile_tmp,row.names = F)
interest_drug <- c("Dihydrorotenone","TAF1","Ulixertinib","Lapatinib","VX.11e",
                   "BMS.754807","NU7441","AZD8055","JQ1","AZD8186")
res_test.order <- res_test[order(res_test$statistic),]
res_test.order <- as.data.frame(res_test.order)

data_boxplot <- data_merge
colnames(data_boxplot) <- str_split(colnames(data_boxplot), "_",simplify = T)[,1]
common_col <- intersect(interest_drug, colnames(data_boxplot));common_col
data_boxplot <- data_boxplot[,c("risk",common_col)]
head(data_boxplot[,1:4])
data_boxplot <- as.data.frame(data_boxplot)

# 每次运行时切换
data_pre <- data_boxplot
colnames(data_pre)
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/",
                      "Boxplot.","drug_auc.","BMS.754807",".pdf");outfile_tmp
pdf(outfile_tmp,width = 3,height = 4,onefile = F)
ggplot(data_pre,aes(x=risk,y= BMS.754807 )) + 
  geom_boxplot(aes(fill=risk), fill = c("#f89f68", "#4b84b3")) + 
  stat_compare_means(method = "t.test",label = "p.format") +
  xlab('risk')+
  ylab('BMS.754807 (IC50)') + 
  theme_classic()
dev.off()