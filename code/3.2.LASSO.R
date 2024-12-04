

library(openxlsx)
library(tidyverse)
library(data.table)
library(glmnet)
library(survival)
library(stringr)
library(broom)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "3.Signature/"

# 读取数据
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

# 导入生存信息
survival_data <- fread("./data/TCGA-COAD.survival.tsv",data.table = F)

# 读取list
interest_genes <- read.csv("./res/3.Signature/table/Cox_univ.hub_genes.csv")
interest_genes <- interest_genes[interest_genes$pvalue <= 0.05,]
interest_genes <- interest_genes[,1];interest_genes
interest_genes <- gsub("[.]","-",interest_genes)

##############################   合并临床数据  ##############################
data_pre <- data_exprs[interest_genes,]
data_pre <- as.data.frame(t(data_pre))
data_surv <- survival_data[,c(1,2,4)]
data_pre <- data.frame(sample = rownames(data_pre),data_pre)
data_pre$sample <- gsub("[.]","-",data_pre$sample)
data_pre$sample <- substring(data_pre$sample,1,16)
data_exprs_clinical <- merge(data_surv, data_pre, by="sample")
data_exprs_clinical <- data.frame(months = data_exprs_clinical$OS.time / 30,
                                  data_exprs_clinical)
data_exprs_clinical <- data_exprs_clinical[data_exprs_clinical$months >= 1,]

data_pre <- data_exprs_clinical[!duplicated(data_exprs_clinical$sample),]
data_pre <- data.frame(data_pre, row.names = data_pre[,2])
data_pre <- na.omit(data_pre)
data_pre <- data_pre[,-2]
colnames(data_pre)
colnames(data_pre)
x <- as.matrix(data_pre[,c(4:ncol(data_pre))])
y <- data.matrix(Surv(data_pre$months,data_pre$OS))
fit <- glmnet(x, y, family = "cox")
cvfit <- cv.glmnet(x, y, family="cox")

tidy_df <- broom::tidy(fit)
tidy_cvdf <- broom::tidy(cvfit)
mypalette <- c(brewer.pal(11,"BrBG"),brewer.pal(11,"Spectral"),brewer.pal(11,"PiYG"),
               brewer.pal(11,"RdBu"),brewer.pal(8,"Accent"))
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","LASSO.","lambda_coef",".pdf");outfile_tmp
pdf(file = outfile_tmp,width = 8,height = 6)
ggplot(tidy_df, aes(lambda, estimate, group = term, color = term)) +
  geom_line(size=1.2)+
  geom_hline(yintercept = 0)+
  scale_x_log10(name = "Log Lambda")+
  ylab("Coefficients")+
  scale_color_manual(name="variable",values = mypalette)+
  theme_bw()
dev.off()

outfile_tmp <- paste0(res_home,proj_name,"/pdf/","LASSO.","lambda_cvfit.2",".pdf");outfile_tmp
pdf(file = outfile_tmp,width = 8,height = 8)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
  
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
coef.df <- as.data.frame(geneCoef)
coef.df$Coef <- as.numeric(coef.df$Coef)
coef.df[,2] <- signif(coef.df[,2],3)

geneCoef[,1] <- gsub("[.]","-",geneCoef[,1])
outfile_tmp <- paste0(res_home,proj_name,"/table/","LASSO.","gene_coef",".csv");outfile_tmp
write.csv(geneCoef,file=outfile_tmp,row.names = F)

# 计算风险得分
riskScore=predict(cvfit, newx = x, s = "lambda.min",type="response")
outCol=c("months","OS",lassoGene)
risk = as.vector(ifelse(riskScore>median(riskScore),"high","low"))
outTab <- cbind(riskScore,risk)
colnames(outTab)[1] <- "risk_score"
outTab <- rownames_to_column(as.data.frame(outTab),"sample")
data_merge <- merge(data_surv, outTab, by="sample")
# save
outfile_tmp <- paste0(res_home,proj_name,"/table/","LASSO.","risk_group",".csv");outfile_tmp
write.csv(data_merge,outfile_tmp,row.names = F)

