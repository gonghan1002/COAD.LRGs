
library(tidyverse)
library(data.table)
library(ggthemes)
library(ggpubr)
library(RColorBrewer)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "4.Biology/"

# 读取数据
data_exprs <- readRDS("./data/TCGA-COAD.TPM.RDS")
head(data_exprs[,1:4])
range(data_exprs)

# group
group <- read.csv("./res/3.Signature/table/Risk_group.csv")
colnames(group)
group <- group[,c("id","risk")]
colnames(group)[c(1,2)] <- c("sample","group")
interest_genes <- c("CD274","PDCD1","CTLA4","LAG3","TNFRSF18","HAVCR2","CD3D","CD3E",
                    "TIGIT","CD80","CD276","CD40","CD4","CD200",
                    "CD160","TNFSF4","TNFRSF8","CD27","CD86")
data_exprs_sub <- data_exprs[interest_genes,]
data_exprs_sub <- as.data.frame(t(data_exprs_sub))
data_pre <- rownames_to_column(data_exprs_sub,var="sample")
data_pre$sample <- substring(data_pre$sample,1,16)
data_pre_merge <- merge(group, data_pre, by="sample")

TME.cells <- colnames(data_pre_merge)[3:ncol(data_pre_merge)];TME.cells
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(data_pre_merge) == TME.cells[i])
  sub <- data.frame(
    PATIENT_ID = data_pre_merge$sample,
    CellType = TME.cells[i],
    group = data_pre_merge$group,
    Composition = data_pre_merge[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
plot.info$group <- factor(plot.info$group,levels = c("high","low"))

outfile_tmp <- paste0(res_home,proj_name, "/pdf/","Boxplot.","diff_ICBs", ".pdf")
pdf(file = outfile_tmp,width = 10,height = 5)
ggboxplot(data = plot.info,x = "CellType",y = "Composition",
          palette = c("#f89f68", "#4b84b3"), fill = "group",
          xlab = "",ylab = "Expression levels") +
  stat_compare_means(label = "p.signif",method = "wilcox.test",
                     aes(group=group),# Pairwise comparison against all
                     hide.ns = F) + 
  theme_base() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()