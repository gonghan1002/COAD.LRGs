
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
proj_name <- "4.Biology/"

checkDir(paste0(res_home,proj_name,"/pdf/"))
checkDir(paste0(res_home,proj_name,"/RData/"))
checkDir(paste0(res_home,proj_name,"/table/"))

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)
group <- read.csv("./res/3.Signature/table/Risk_group.csv")
colnames(group)[1] <- "sample"
group_sub <- group[,c(1,8)]
group_sub <- group_sub[!duplicated(group_sub$sample),]
colnames(data_exprs) <- substring(colnames(data_exprs),1,16)
group_sub <- group_sub[group_sub$sample %in% colnames(data_exprs),]
data_exprs <- data_exprs[,group_sub$sample]
data_anno <- data.frame(group = group[,c("risk")],row.names = group$sample)

# match
data_exprs <- data_exprs[,rownames(data_anno)]

# degs
degs_res <- diff_limma(inputData = data_exprs, colData = data_anno,
                       contrast_name = c("high", "low"),
                       log2trans = F, adjust.method = "BH", 
                       pval_cutoff = NULL,log2fc_cutoff = NULL)
degs_res$threshold = factor(ifelse(degs_res$adj.P.Val < 0.05 & abs(degs_res$logFC) >= 1, 
                                   ifelse(degs_res$logFC >= 1,'Up','Down'),'NoSig'),
                            levels=c('Up','Down','NoSig'))
table(degs_res$threshold != "NoSig" )
outfile_tmp <- paste0(res_home, proj_name,"/table/","degs.","limma.","fc2", ".csv");outfile_tmp
write.csv(degs_res, file = outfile_tmp)
data_plot <- degs_res
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","Volcano.","degs.", "limma.","fc2", ".pdf");outfile_tmp
pdf(file = outfile_tmp,onefile = F,width = 6,height = 6)
ggplot(data_plot,aes(x = logFC, y = -log10(adj.P.Val)))+
  geom_point_rast(size = 1.5, aes(color = threshold))+
  scale_color_manual(values = c("#f89f68","#4b84b3","gray50"),
                     labels = c("Up","Down","NoSig"))+ #确定点的颜色
  geom_vline(xintercept=c(-1,1),linetype ="longdash",col="black",lwd=0.5) + 
  geom_hline(yintercept = -log10(0.05),linetype ="longdash",col="black",lwd=0.5)+ 
  theme_bw()+ #修改图片背景
  theme(
    legend.title = element_blank(), # 不显示图例标题
    panel.grid=element_blank(),
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.key.size = unit(24, "pt")
  )
dev.off()
