

library(openxlsx)
library(tidyverse)
library(data.table)
library(rrtable)
library(ggplot2)
library(survival)
library(survminer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "3.Signature/"

# 读入临床信息
data_clinical <- fread("./data/TCGA-COAD.clinical_info.csv",data.table = F)
colnames(data_clinical)
data_surv <- fread("./data/TCGA-COAD.survival.tsv", data.table = F)
colnames(data_surv)
group <- read.csv("./res/3.Signature/table/Risk_group.csv")
group$id <- gsub("[.]","-",group$id)
data_surv <- fread("./data/TCGA-COAD.survival.tsv",data.table = F)
colnames(data_surv)
colnames(group)
group_sub <- group[,c("id","risk")]
colnames(group_sub)[1] <- "sample"
data_surv <- data_surv[,-3]
data_merge <- merge(data_surv, group_sub, "sample")
colnames(data_merge)[3] <- 'months'
data_merge$months <- data_merge$months / 30
data_merge <- data_merge[data_merge$months >= 1,]
data_merge$sample <- gsub("-",".",data_merge$sample)
colnames(data_clinical)
data_clinical_sub <- data_clinical[,c(1,2,3,5,6,8:11)]
colnames(group_sub)[1] <- "sample"
data_merge <- merge(data_merge,data_clinical_sub,by="sample")
colnames(data_merge) 
str(data_merge)
data_merge$age <- ifelse(data_merge$age > 60,"> 60","<= 60")

### 每次运行时切换临床参数
data_clinical <- data_merge
str(data_clinical)

pres_idx1 <- which(data_clinical$age %in% "> 60")
pres_idx2 <- which(data_clinical$age %in% "<= 60")
data_clinical <- data_clinical[c(pres_idx1),]

### plot
fit <- survfit(Surv(months, OS) ~ risk, data = data_clinical)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Surv.","age.","than_60",".pdf");outfile_tmp
pdf(file=outfile_tmp, onefile = F,width = 5,height =6)   
ggsurvplot(fit, 
           data=data_clinical,
           conf.int=F,
           pval=T,
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("High", "Low"), #, "Cluster4"
           legend.title="Risk group",
           xlab="Time(months)",
           ylab="Overall Survival",
           break.time.by = 36,
           risk.table.title="",
           palette=c("#f89f68", "#4b84b3"), #, "darkblue"
           risk.table.height=.3)
dev.off()

