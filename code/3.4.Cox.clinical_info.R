

library(data.table)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(forestplot)

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
group <- read.csv("./res/3.Signature/table/Risk_group.csv")
data_surv <- fread("./data/TCGA-COAD.survival.tsv",data.table = F)
colnames(data_surv)
colnames(group)
group_sub <- group[,c("id","riskScore")]
colnames(group_sub)[1] <- "sample"
data_surv <- data_surv[,-3]
data_surv$sample <- gsub("-",".",data_surv$sample)
data_merge <- merge(data_surv, group_sub, "sample")
colnames(data_merge)[3] <- 'months'
data_merge$months <- data_merge$months / 30
data_merge <- data_merge[data_merge$months >= 1,]
colnames(data_clinical)
data_clinical_sub <- data_clinical[,c(1,2,3,5,6,8:11)]

# merge
colnames(group_sub)[1] <- "sample"
data_merge <- merge(data_merge,data_clinical_sub,by="sample")
colnames(data_merge) 
colnames(data_merge)[4] <- c("risk_score")
str(data_merge)

# age
data_merge$age <- ifelse(data_merge$age >= 60,">= 60","< 60")

data_clinical <- data_merge
{
# age
tmp_idx <- data_clinical$age
table(tmp_idx)
tmp_idx[which(tmp_idx=="< 60")] <- 0
tmp_idx[which(tmp_idx==">= 60")] <- 1
table(tmp_idx)
data_clinical$age <- as.numeric(tmp_idx)

# gender
tmp_idx <- data_clinical$gender
table(tmp_idx)
tmp_idx[which(tmp_idx=="female")] <- 0
tmp_idx[which(tmp_idx=="male")] <- 1
table(tmp_idx)
data_clinical$gender <- as.numeric(tmp_idx)

# stage
# 一定要注意 遵守iv iii ii i顺序
tmp_idx <- data_clinical$stage
table(tmp_idx)
tmp_idx[which(tmp_idx=="Stage I")] <- 1
tmp_idx[which(tmp_idx=="Stage II")] <- 1
tmp_idx[which(tmp_idx=="Stage III")] <- 2
tmp_idx[which(tmp_idx=="Stage IV")] <- 2
table(tmp_idx)
data_clinical$stage <- as.numeric(tmp_idx)

# T
tmp_idx <- data_clinical$stage_T
table(tmp_idx)
tmp_idx[which(tmp_idx=="T1")] <- 1
tmp_idx[which(tmp_idx=="T2")] <- 1
tmp_idx[which(tmp_idx=="T3")] <- 2
tmp_idx[which(tmp_idx=="T4")] <- 2
table(tmp_idx)
data_clinical$stage_T <- as.numeric(tmp_idx)

# N
tmp_idx <- data_clinical$stage_N
table(tmp_idx)
tmp_idx[which(tmp_idx=="N0")] <- 0
tmp_idx[which(tmp_idx=="N1")] <- 1
tmp_idx[which(tmp_idx=="N2")] <- 1
# tmp_idx[which(tmp_idx=="N3")] <- 3
table(tmp_idx)
data_clinical$stage_N <- as.numeric(tmp_idx)

# N
tmp_idx <- data_clinical$stage_M
table(tmp_idx)
tmp_idx[which(tmp_idx=="M0")] <- 0
tmp_idx[which(tmp_idx=="M1")] <- 1
table(tmp_idx)
data_clinical$stage_M <- as.numeric(tmp_idx)

# venous_invasion
tmp_idx <- data_clinical$venous_invasion
table(tmp_idx)
tmp_idx[which(tmp_idx=="NO")] <- 0
tmp_idx[which(tmp_idx=="YES")] <- 1
table(tmp_idx)
data_clinical$venous_invasion <- as.numeric(tmp_idx)

# lymphatic_invasion
tmp_idx <- data_clinical$lymphatic_invasion
table(tmp_idx)
tmp_idx[which(tmp_idx=="NO")] <- 0
tmp_idx[which(tmp_idx=="YES")] <- 1
table(tmp_idx)
data_clinical$lymphatic_invasion <- as.numeric(tmp_idx)
}

str(data_clinical)

HR_res_list <- cox_regression(input=data_clinical,
                              var_cols=c("stage",
                                         "stage_T","stage_N","stage_M",
                                         "risk_score"),
                              time_col="months", status_col="OS")
hr_res_univ <- HR_res_list$univ_res
hr_res_multiv <- HR_res_list$multiv_res

sub_inf <- function(x, sub_method=max){
  sub_value <- sub_method(x[which(is.finite(x))], na.rm=TRUE)
  inf_idx <- which(is.infinite(x))
  if(length(inf_idx)>0)
    x[inf_idx] <- sub_value
  return(x)
}

HR_df = hr_res_multiv
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Cox.","multiv.","clinical_factor",".pdf");outfile_tmp
zero_value=1; auto_col=TRUE; width = 6; height = 4; outfile <- outfile_tmp

HR_df[, 3] <- sub_inf(as.numeric(HR_df[, 3]), sub_method=min)
HR_df[, 4] <- sub_inf(as.numeric(HR_df[, 4]), sub_method=max)
HR_df[, 5] <- round(HR_df[, 5], digits=4)
HR_df_table <- as.matrix(HR_df[,c(1, 6, 5)])
HR_df_table <- rbind(colnames(HR_df_table), HR_df_table)
HR_df2 <- rbind(NA, HR_df[, c(2:4)])
HR_df2 <- data.frame(HR_df2, box_col='#f89f68') # 第一次univ #8DD3C7 第二次 multiv #E41A1C
if(auto_col){
  box_col <- as.character(HR_df2[, "box_col"])
  box_col[which(HR_df2[, "HR"]<1)] <- "blue"
}
clip_tmp_min <- min(c(0, HR_df2[, 2]), na.rm=TRUE)
clip_tmp_max <- max(c(2, HR_df2[, 3]), na.rm=TRUE)
clip_tmp <- c(clip_tmp_min, clip_tmp_max)
pdf(outfile, onefile=FALSE, width = width, height = height)
forestplot(labeltext = HR_df_table, 
           mean = as.numeric(HR_df2[, 1]), lower = as.numeric(HR_df2[, 2]), upper = as.numeric(HR_df2[, 3]),
           is.summary = c(TRUE, rep(FALSE, nrow(HR_df2)-1)), # 第一行字体加粗
           align=c("l", "c", "c"),
           clip=clip_tmp,  
           xlog=FALSE, 
           zero = zero_value,
           boxsize = 0.2, 
           new_page = TRUE,
           txt_gp =fpTxtGp(ticks=gpar(cex=0.8)), # x轴字体
           lwd.zero = 1,
           lwd.ci = 1,
           col=fpColors(box=as.character(HR_df2[,"box_col"]), 
                        summary= "#E41A1C",lines = 'lightgray',zero = 'lightgray'))
dev.off()
