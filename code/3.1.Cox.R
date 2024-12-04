########################   Cox  ########################   

library(data.table)
library(tidyverse)
library(stringr)
library(openxlsx)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "3.Signature/"

checkDir(paste0(res_home,proj_name,"/pdf/"))
checkDir(paste0(res_home,proj_name,"/RData/"))
checkDir(paste0(res_home,proj_name,"/table/"))

# 读取数据
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

# 导入生存信息
survival_data <- fread("./data/TCGA-COAD.survival.tsv",data.table = F)

# 读取list
interest_genes <- read.csv("./res/2.Cluster/table/degs.Cluster1.vs.merge_Cluster2.merge.fc1.5.csv")
interest_genes <- interest_genes[interest_genes$threshold != "NoSig",][,1]
interest_genes

##############################   合并临床数据  ##############################
data_pre <- data_exprs[interest_genes,]
data_pre <- as.data.frame(t(data_pre))
colnames(survival_data)
data_pre <- data.frame(sample = rownames(data_pre),data_pre)
data_pre$sample <- substring(data_pre$sample,1,16)

# 修改sample格式 方便合并
data_pre$sample <- gsub("[.]","-",data_pre$sample)
# merge
data_merge <- merge(survival_data, data_pre, by="sample")

##############################   cox回归  ##############################
target_genes <- gsub("-",".",interest_genes);target_genes

# cox回归
HR_res_list_ICB <- cox_regression(input=data_merge,var_cols=target_genes, 
                                  time_col="OS.time", status_col="OS")
HR_res_univ <- HR_res_list_ICB$univ_res # 提取单因素回归结果
colnames(HR_res_univ)[1] <- "genes"
# 选择预后基因
HR_res_univ_sig <- HR_res_univ[HR_res_univ$pvalue < 0.05,]
# save
outfile_tmp <- paste0(res_home,proj_name,"/table/","Cox_univ.","hub_genes",".csv");outfile_tmp
write.csv(HR_res_univ,file = outfile_tmp,row.names = F)

## plot
data_plot <- HR_res_univ_sig
data_plot$risk_froup <- ifelse(data_plot$HR>1,"risk","protect")
data_plot <- data_plot[order(data_plot$HR,decreasing = F),]
data_plot$p_signif <- ifelse(data_plot$pvalue<0.001,"***",
                              ifelse(data_plot$pvalue<0.01,"**",
                                     ifelse(data_plot$pvalue<0.05,"*","ns")))
data_plot$genes <- factor(data_plot$genes,levels = data_plot$genes)

# save
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Cox_univ.","hub_genes.2",".pdf")
pdf(file = outfile_tmp,width = 4,height = 6)
ggplot(data_plot, aes(HR, genes, col=risk_froup))+ # 不同形状shape= Factor
  geom_point(size=3) +
  geom_errorbarh(aes(xmax =HR_CI_975, xmin = HR_CI_025), height = 0.4) +
  #scale_x_log10()+
  scale_x_continuous(limits= c(0.9, 2.1)) + #, breaks= seq(0, 2.5, 0.5)
  scale_color_manual(limits=c("risk","protect"),values=c("#f89f68","#4b84b3"),name="Risk Factor") + 
  scale_y_discrete(limits=data_plot$genes)+
  geom_vline(aes(xintercept = 1)) +
  geom_hline(aes(yintercept = 0.5))+
  # geom_text(data=data_plot,aes(x= 2.1,y=genes,label=p_signif),color="black")+
  xlab('Hazard Ratio ') + ylab(' ')+theme_bw()
dev.off()


