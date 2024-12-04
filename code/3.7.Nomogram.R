
library(rms)
library(foreign)
library(survival)
library(regplot)
library(mstate)
library(timeROC)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "3.Signature/"

#### 读入临床与生存信息
# 读入临床信息
data_clinical <- fread("./data/TCGA-COAD.clinical_info.csv",data.table = F)
colnames(data_clinical)

# 读取分组数据
group <- read.csv("./res/3.Signature/table/Risk_group.csv")
# group$sample <- gsub("[.]","-",group$id)
data_surv <- fread("./data/TCGA-COAD.survival.tsv",data.table = F)
colnames(data_surv)

### 合并
colnames(group)
group_sub <- group[,c("id","risk")]
colnames(group_sub)[1] <- "sample"
data_surv <- data_surv[,-3]
data_surv$sample <- gsub("-",".",data_surv$sample)
data_merge <- merge(data_surv, group_sub, "sample")
colnames(data_merge)[3] <- 'months'
data_merge$months <- data_merge$months / 30
data_merge <- data_merge[data_merge$months >= 1,]
colnames(data_clinical)
data_clinical_sub <- data_clinical[,c(1,8,9,11)]

# merge
colnames(group_sub)[1] <- "sample"
data_merge <- merge(data_merge,data_clinical_sub,by="sample")
colnames(data_merge) 
colnames(data_merge)[4] <- c("risk")
str(data_merge)

data_clinical <- data_merge

{
  # stage
  # 一定要注意 遵守iv iii ii i顺序
  tmp_idx <- data_clinical$stage
  table(tmp_idx)
  tmp_idx[which(tmp_idx=="Stage I")] <- "Stage I-II"
  tmp_idx[which(tmp_idx=="Stage II")] <- "Stage I-II"
  tmp_idx[which(tmp_idx=="Stage III")] <- "Stage III-IV"
  tmp_idx[which(tmp_idx=="Stage IV")] <- "Stage III-IV"
  table(tmp_idx)
  data_clinical$stage <- tmp_idx
  
  # T
  tmp_idx <- data_clinical$stage_T
  table(tmp_idx)
  tmp_idx[which(tmp_idx=="T1")] <- "T1-2"
  tmp_idx[which(tmp_idx=="T2")] <- "T1-2"
  tmp_idx[which(tmp_idx=="T3")] <- "T3-4"
  tmp_idx[which(tmp_idx=="T4")] <- "T3-4"
  table(tmp_idx)
  data_clinical$stage_T <- tmp_idx
  
  # M
  tmp_idx <- data_clinical$stage_M
  table(tmp_idx)
  tmp_idx[which(tmp_idx=="M0")] <- "M0"
  tmp_idx[which(tmp_idx=="M1")] <- "M1"
  table(tmp_idx)
  data_clinical$stage_M <- tmp_idx
}

ddist <- datadist(data_clinical)
options(datadist='ddist')
rt <- data_clinical
colnames(rt)
cox <- cph(Surv(months,OS) ~ stage+stage_T+stage_M+risk,surv=T,x=T, y=T,data=rt) 

surv <- Survival(cox)
sur_1_year<-function(x)surv(1*12,lp=x)
sur_3_year<-function(x)surv(1*12*3,lp=x)
sur_5_year<-function(x)surv(1*12*5,lp=x)
nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),
                    lp= F,funlabel=c('1-Year survival','3-Year survival','5-Year survival'),
                    maxscale=100,
                    fun.at= c('1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'))
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Nomogram.","classical",".pdf");outfile_tmp
pdf(outfile_tmp, height = 6, width = 8, onefile = F)
plot(nom_sur,xfrac=0.4)
dev.off()


###########   Calibration   ###########
### 1-year
cox_1 <- cph(Surv(months,OS) ~ stage+stage_T+stage_M+risk,surv=T,x=T, y=T,time.inc = 1*12,data=rt) 
cal_1 <- calibrate(cox_1, cmethod="KM", method="boot", u=1*12, m= 100, B=1000)
### 2-year
cox_2 <- cph(Surv(months,OS) ~ stage+stage_T+stage_M+risk,surv=T,x=T, y=T,time.inc = 3*12,data=rt) 
cal_2 <- calibrate(cox_2, cmethod="KM", method="boot", u=3*12, m= 100, B=1000)
### 3-year
cox_3 <- cph(Surv(months,OS) ~ stage+stage_T+stage_M+risk,surv=T,x=T, y=T,time.inc = 5*12,data=rt) 
cal_3 <- calibrate(cox_3, cmethod="KM", method="boot", u=5*12, m= 100, B=1000)


### plot
my_color <- RColorBrewer::brewer.pal(3,"Set1")
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Nomogram.","calibration",".pdf");outfile_tmp
pdf(outfile_tmp,onefile = F,width = 8,height = 8)
plot(cal_1, lwd=3, lty=2, errbar.col="black",col= my_color[1],
     xlim = c(0,1.0), ylim = c(0,1.0),
     xlab ="Nomogram-predicted OS",ylab="Observed OS")
par(new = T)
plot(cal_2, lwd=3, lty=2, errbar.col="black",col= my_color[2],
     xlim = c(0,1.0), ylim = c(0,1.0),
     xlab ="",ylab="")
par(new = T)
plot(cal_3, lwd=3, lty=2, errbar.col="black",col= my_color[3],
     xlim = c(0,1.0), ylim = c(0,1.0),
     xlab ="",ylab="")
lines(cal_1[,c('mean.predicted','KM')],type = 'l',lwd = 2,col =my_color[1] ,pch = 16)
lines(cal_2[,c('mean.predicted','KM')],type = 'l',lwd = 2,col =my_color[2] ,pch = 16)
lines(cal_3[,c('mean.predicted','KM')],type = 'l',lwd = 2,col =my_color[3] ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 2,lwd = 2,col = "black")
legend(0.8,0.2,
       c("1-year","3-year","5-year"), 
       lty = c(1,1,1), lwd = c(3,3,3), 
       col = my_color, 
       bty = "n")
dev.off()

f<-coxph(Surv(months,OS==1)~stage+stage_T+stage_M+risk,data=rt)
sum.surv<-summary(f)
c_index<-sum.surv$concordance
c_index  ##


### ROC
pret <- predict(cox,data_clinical,type = "lp")
ROC <- data.frame(months = data_clinical[,"months"],status = data_clinical[,"OS"],
                  score = pret)
ROC <- na.omit(ROC)
res_ROC <- timeROC(T = ROC$months,
                   delta = ROC$status,
                   marker = ROC$score, 
                   cause=1, 
                   weighting="marginal", 
                   times=c(1*12,2*12,3*12), 
                   ROC = TRUE,
                   iid = TRUE)
ROC <- res_ROC
auc_1 = ROC$AUC[[1]]; auc_2 = ROC$AUC[[2]]; auc_3 = ROC$AUC[[3]]
auc_1; auc_2; auc_3

dat = data.frame(tpr_1 = ROC$TP[,1],fpr_1 = ROC$FP[,1],tpr_2 = ROC$TP[,2],
                 fpr_2 = ROC$FP[,2],tpr_3 = ROC$TP[,3],fpr_3 = ROC$FP[,3])

my_palet <- c("#f89f68", "#4b84b3","#00C8B5")
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","ROC.","nomogram.","years.","1_2_3",".pdf")
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
  annotate("text",x = .75, y = .15,label = paste("AUC of 2 years = ",round(auc_2,2)),color = my_palet[2])+
  annotate("text",x = .75, y = .05,label = paste("AUC of 3 years = ",round(auc_3,2)),color = my_palet[3])+
  scale_x_continuous(name  = "1-Specificity")+
  scale_y_continuous(name = "Specificity")
dev.off()

