########################    准备TIDE数据    ######################## 

library(data.table)
library(tidyverse)
library(stringr)

rm(list = ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

# 去掉在一半样品中 表达量 < 1 的基因
pres = apply(data_exprs>1, 1, sum) 
to_keep = pres > 0.5 * ncol(data_exprs)
table(to_keep) # 看一下有多少满足条件
data_pre = data_exprs[to_keep,]

### 针对每一个基因 与 平均表达值进行比较
data_norm <- apply(data_pre,1,function(x){x - mean(x)})
data_norm <- as.data.frame(t(data_norm))

# NOTE: 需要保存为tab分隔符文件
# data_new <- rownames_to_column(data_norm, var = "SYMBOL")

# save
write.table(data_norm, "./TME/TCGA-COAD.TIDE.txt",sep = "\t")


