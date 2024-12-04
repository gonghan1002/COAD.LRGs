########################   cluster  ######################## 
#### 不行 换一致性聚类

library(data.table)
library(tidyverse)
library(openxlsx)
library(ConsensusClusterPlus)
library(survival)
library(survminer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "2.Cluster/"

checkDir(paste0(res_home,proj_name,"/pdf/"))
checkDir(paste0(res_home,proj_name,"/RData/"))
checkDir(paste0(res_home,proj_name,"/table/"))

# read data
group <- read.csv("./res/2.Cluster/table/Cluster.Consensus.k3.csv")

# surv data
data_surv <- fread("./data/TCGA-COAD.survival.tsv", data.table = F)

# 合并cluster2 3
colnames(group) <- c("sample","group")
group$group <- paste0("Cluster", group$group)
table(group$group)

group$group[group$group %in% c("Cluster2","Cluster3")] <- "Cluster2"
table(group$group)

##############################  surv ##############################
#### 合并生存数据
group$sample <- gsub("\\.","-",group$sample)
group$sample <- substring(group$sample,1,16)
data_merge <- merge(data_surv, group,by = "sample")

# 过滤 删除生存时间少于30天的患者
data_merge$months <- data_merge$OS.time/30
data_merge <- data_merge[data_merge$months >= 1,]
str(data_merge)

### plot
fit <- survfit(Surv(months, OS) ~ group, data = data_merge)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Surv.","consensus.k2",".pdf");outfile_tmp
pdf(file=outfile_tmp, onefile = F,width = 5,height =6)
ggsurvplot(fit, 
           data=data_merge,
           conf.int=F,
           pval = T,
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("Cluster1", "Cluster2"), #, "Cluster4"
           legend.title="Cluster",
           xlab="Time(months)",
           ylab="Overall Survival",
           break.time.by = 24,
           risk.table.title="",
           palette=c("#f89f68", "#4b84b3"), #, "darkblue"
           risk.table.height=.3)
dev.off()





