
library(openxlsx)
library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(stringr)
library(timeROC)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "3.Signature/"

# read data
group <- read.csv("./res/3.Signature/table/Risk_group.csv")

# data handle
colnames(group)
melcli_group <- group[,-c(4,5,6)]
colnames(melcli_group)
melcli_group <- melcli_group[melcli_group$months >= 1,]

diff = survdiff(Surv(months, event) ~ risk,data = melcli_group)
pValue = 1-pchisq(diff$chisq,df=1)
pValue = signif(pValue,2);pValue


fit <- survfit(Surv(months, event) ~ risk, data = melcli_group)
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","Surv.","risk_group",".pdf")
pdf(file=outfile_tmp,onefile=F,height=6,width=5)
ggsurvplot(fit, data = melcli_group, 
           conf.int=F, 
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=T,   # table
           risk.table.height=.3,
           legend.labs = c("high risk", "low risk"),
           legend.title="risk",
           xlab="Time(months)", ylab="Overall Survival",
           break.time.by = 36, 
           risk.table.title="",
           palette = c("#f89f68", "#4b84b3"))
dev.off()

melcli_group <- group[,-c(4,5,6)]
colnames(melcli_group)
melcli_group <- melcli_group[melcli_group$months >= 1,]

ROC<- timeROC(T = melcli_group$months, 
              delta = melcli_group$event, 
              marker = melcli_group$riskScore,
              cause=1, 
              weighting="marginal", 
              times=c(1*12,3*12,5*12),
              ROC = TRUE,
              iid = TRUE)
auc_1 = ROC$AUC[[1]]; auc_2 = ROC$AUC[[2]]; auc_3 = ROC$AUC[[3]]
auc_1; auc_2; auc_3

dat = data.frame(tpr_1 = ROC$TP[,1],fpr_1 = ROC$FP[,1],tpr_2 = ROC$TP[,2],
                 fpr_2 = ROC$FP[,2],tpr_3 = ROC$TP[,3],fpr_3 = ROC$FP[,3])

my_palet <- c("#f89f68", "#4b84b3","#00C8B5")
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","ROC.","risk_group.","years.","1_3_5",".pdf");outfile_tmp
pdf(file = outfile_tmp, width=4.5, height=4)
ggplot() + 
  geom_smooth(data = dat,aes(x = fpr_1, y = tpr_1),color = my_palet[1]) + 
  geom_smooth(data = dat,aes(x = fpr_2, y = tpr_2),color = my_palet[2])+
  geom_smooth(data = dat,aes(x = fpr_3, y = tpr_3),color = my_palet[3])+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme(panel.background=element_rect(colour=NA,fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"))+
  annotate("text",x = .75, y = .25,label = paste("AUC of 1 years = ",round(auc_1,2)),color = my_palet[1])+
  annotate("text",x = .75, y = .15,label = paste("AUC of 3 years = ",round(auc_2,2)),color = my_palet[2])+
  annotate("text",x = .75, y = .05,label = paste("AUC of 5 years = ",round(auc_3,2)),color = my_palet[3])+
  scale_x_continuous(name  = "1-Specificity")+
  scale_y_continuous(name = "Specificity")
dev.off()
