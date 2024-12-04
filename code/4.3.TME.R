
library(tidyverse)
library(data.table)
library(ggthemes)
library(ggpubr)
library(openxlsx)

rm(list = ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "4.Biology/"

# 读入cibersort res
CIBERSORT_res <- read.csv("./TME/CIBERSORT_res.TCGA-COAD.csv")

# 读入estimate res
ESTIMATE_res <- read.csv("./TME/ESTIMATE_res.TCGA-COAD.csv")

# group
group <- read.csv("./res/3.Signature/table/Risk_group.csv")
colnames(group)
group <- group[,c("id","risk")]
colnames(group)[c(1,2)] <- c("sample","group")

################################   比较ESTIMATE  ################################ 
data_pre <- ESTIMATE_res
colnames(data_pre)[1] <- "sample"

### merge
data_pre$sample <- substring(data_pre$sample,1,16)
data_merge <- merge(group, data_pre, by="sample")

### 转化数据类型
colnames(data_merge)
TME.cells <- colnames(data_merge)[3:ncol(data_merge)][-4];TME.cells
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(data_merge) == TME.cells[i])
  sub <- data.frame(
    PATIENT_ID = data_merge$sample,
    CellType = TME.cells[i],
    group = data_merge$group,
    Composition = data_merge[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
plot.info$group <- factor(plot.info$group,levels = c("high","low"))
outfile_tmp <- paste0(res_home,proj_name, "/pdf/","Boxplot.","diff_ESTIMATE", ".pdf")
pdf(file = outfile_tmp,width = 5,height = 4)
ggboxplot(data = plot.info,x = "CellType",y = "Composition",
          palette = c("#f89f68", "#4b84b3"), fill = "group",
          xlab = "",ylab = "scores") +
  stat_compare_means(label = "p.signif",method = "t.test",
                     aes(group=group),
                     hide.ns = F) + 
  theme_base() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()


################################   比较CIBERSORT  ################################ 
data_pre <- CIBERSORT_res
colnames(data_pre)
colnames(data_pre)[1] <- "sample"
data_pre <- data_pre[data_pre$P.value < 0.05,]
data_pre <- data_pre[,-c(24,25,26)]

### merge
data_pre$sample <- substring(data_pre$sample,1,16)
data_pre_merge <- merge(group, data_pre, by="sample")
colnames(data_pre_merge)

TME.cells <- colnames(data_pre_merge)[3:24];TME.cells
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
outfile_tmp <- paste0(res_home,proj_name, "/pdf/","Boxplot.","diff_CIBERSORT", ".pdf")
pdf(file = outfile_tmp,width = 12,height = 6)
ggboxplot(data = plot.info,x = "CellType",y = "Composition",
          palette = c("#f89f68", "#4b84b3"), fill = "group",
          xlab = "",ylab = "Immune cell composition") +
  stat_compare_means(label = "p.signif",method = "wilcox.test",
                     aes(group=group),# Pairwise comparison against all
                     hide.ns = F) + 
  theme_base() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()
