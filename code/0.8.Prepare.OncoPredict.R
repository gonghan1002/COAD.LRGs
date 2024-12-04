########################  药物敏感性分析 ######################## 

library(data.table)
library(tidyverse)
library(openxlsx)
library(oncoPredict)
library(gtools)
library(reshape2)
library(ggpubr)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

### 读取下载的GDSC2数据
GDSC2_Expr <-  readRDS("E:/r_project/data/oncoPredict/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res <- readRDS("E:/r_project/data/oncoPredict/DataFiles/Training Data/GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) 

# 读取数据
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

##############################  药物预测  ############################## 
### oncoPredict
data_pre <- data_exprs
data_pre <- as.matrix(data_pre)

### 计算IC50 不需要过滤数据 包有这个流程
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = data_pre,
              batchCorrect = 'eb',
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,  
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )



