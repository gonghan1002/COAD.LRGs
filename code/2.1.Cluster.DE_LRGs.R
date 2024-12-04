########################   cluster  ######################## 

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
data_exprs <- readRDS("./data/TCGA-COAD.TPM.Tumor.RDS")
head(data_exprs[,1:4])
range(data_exprs)

# surv data
data_surv <- fread("./data/TCGA-COAD.survival.tsv", data.table = F)

# 读取目的基因
interest_genes <- read.csv("./res/1.degs/table/degs.LRGs.fc1.csv")
interest_genes <- interest_genes[interest_genes$threshold != "NoSig",][,1]
interest_genes

#### cluster
data_sub <- data_exprs[interest_genes,]

# Consensus Cluster
data_consens <- as.matrix(data_sub)
results = ConsensusClusterPlus(data_consens, maxK=6, reps=100,
                               pItem=0.8,pFeature=1,
                               title= "./res/2.Cluster/pdf/",
                               clusterAlg="pam",distance="spearman",
                               seed=1262118388.71279,plot="pdf")

##############################  不同cluster与生存的关系  ##############################
cluster_res <- results[[3]]
cluster_class <- cluster_res$consensusClass
outfile_tmp <- paste0(res_home,proj_name,"/table/","Cluster.","Consensus.k3",".csv");outfile_tmp
write.csv(cluster_class, file = outfile_tmp, row.names = T)

# 提取聚类样本
cluster_class <- as.data.frame(cluster_class)
table(cluster_class[,1])
cluster_class <- rownames_to_column(cluster_class,var = "sample")

#### 合并生存数据
cluster_class$sample <- gsub("\\.","-",cluster_class$sample)
cluster_class$sample <- substring(cluster_class$sample,1,16)
data_merge <- merge(data_surv, cluster_class,by = "sample")

data_merge$months <- data_merge$OS.time/30
data_merge <- data_merge[data_merge$months >= 1,]
str(data_merge)
data_merge$group <- paste0("Cluster",data_merge$cluster_class)

### survival
diff = survdiff(Surv(months, OS) ~ group,data = data_merge)
pValue = 1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,2)
pValue

##############################  plot  ##############################
fit <- survfit(Surv(months, OS) ~ group, data = data_merge)
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Surv.","consensus.k3",".pdf");outfile_tmp
pdf(file=outfile_tmp, onefile = F,width = 5,height =6)
ggsurvplot(fit, 
           data=data_merge,
           conf.int=F,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("Cluster1", "Cluster2","Cluster3"), #, "Cluster4"
           legend.title="Cluster",
           xlab="Time(months)",
           ylab="Overall Survival",
           break.time.by = 24,
           risk.table.title="",
           palette=c("#f89f68", "#4b84b3","#00C8B5"), #, "darkblue"
           risk.table.height=.3)
dev.off()





