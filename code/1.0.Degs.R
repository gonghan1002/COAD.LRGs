########################   差异分析   ######################## 

library(data.table)
library(tidyverse)
library(stringr)
library(openxlsx)
library(ggrastr)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "1.degs/"

checkDir(paste0(res_home,proj_name,"/pdf/"))
checkDir(paste0(res_home,proj_name,"/RData/"))
checkDir(paste0(res_home,proj_name,"/table/"))

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

# interest gene
interest_genes <- read.xlsx("./list/lactylation-related genes.xlsx")
interest_genes <- c(unlist(interest_genes),use.names = F)
interest_genes <- interest_genes[!duplicated(interest_genes)]

### 区分tumor normal
### group 
# 观察样本类型
table(substring(colnames(data_exprs),first = 13,last = 16))
# 不要01B B是石蜡样本
tumor_idx <- grep("01A",substring(colnames(data_exprs),first = 13,last = 16))
normol_idx <- grep("11A",substring(colnames(data_exprs),first = 13,last = 16))
# tumor 在前， normal 在后
data_exprs <- data_exprs[,c(tumor_idx, normol_idx)]

### 构建注释组
table(substring(colnames(data_exprs),first = 13,last = 16))
condition <- c(rep("tumor", 462),rep("normal", 41))
data_anno <- data.frame(condition, row.names = colnames(data_exprs))

### degs
# data_pre <- log2(data_exprs + 1)
res.degs <- diff_limma(inputData = data_exprs,colData = data_anno,
                       contrast_name = c("tumor","normal"),
                       log2trans = F,adjust.method = "BH",
                       pval_cutoff = NULL,log2fc_cutoff = NULL)

# 加标签列
res.degs$threshold = factor(ifelse(res.degs$adj.P.Val < 0.05 & abs(res.degs$logFC) >= 1, 
                                   ifelse(res.degs$logFC >= 1,'Up','Down'),'NoSig'),
                            levels=c('Up','Down','NoSig'))
table(res.degs$threshold)

# save
outfile_tmp <- paste0(res_home,proj_name, "/table/","degs.","tumor_normal.", "fc1", ".csv");outfile_tmp
write.csv(res.degs,file = outfile_tmp)

# sub
data_sub <- res.degs[rownames(res.degs) %in% interest_genes,]
table(data_sub$threshold)
# save
outfile_tmp <- paste0(res_home,proj_name, "/table/","degs.LRGs.", "fc1", ".csv");outfile_tmp
write.csv(data_sub,file = outfile_tmp)

#################     火山图  volcano     #####################
data_plot <- res.degs
colnames(data_plot)

# save
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","Volcano.","tumor_normal.", "fc2", ".pdf");outfile_tmp
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

