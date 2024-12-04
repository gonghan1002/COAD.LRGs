######################## 

library(FactoMineR)
library(factoextra)
library(missMDA)
library(RColorBrewer)
library(data.table)
library(tidyverse)
library(openxlsx)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "2.Cluster/"

# read data
data_exprs <- readRDS("./data/TCGA-COAD.TPM.Tumor.RDS")
head(data_exprs[,1:4])

# 读取分组数据
group <- read.csv("./res/2.Cluster/table/Cluster.Consensus.k3.csv")
# 合并cluster2 3
colnames(group) <- c("sample","group")
group$group <- paste0("Cluster", group$group)
table(group$group)

group$group[group$group %in% c("Cluster2","Cluster3")] <- "Cluster2"
table(group$group)

# 读取list 这里是提取预后基因
interest_genes <- read.csv("./res/1.degs/table/degs.LRGs.fc1.csv")
interest_genes <- interest_genes[interest_genes$threshold != "NoSig",][,1]
interest_genes


##############################    t-SNE    ##############################
tsne_out <- Rtsne::Rtsne(t(data_pre), pca = T)
data_tsne <- data.frame(tsne_out$Y,group$group)
colnames(data_tsne) <- c("tSNE_1","tSNE_2","group")

# plot
outfile_tmp <- paste0(res_home,proj_name,"/pdf/","tSNE.Consensus.","cluster.k2",".pdf");outfile_tmp
pdf(file=outfile_tmp, onefile = F,width = 6,height =5)
ggplot(data_tsne,aes(tSNE_1,tSNE_2,fill=group))+
  geom_point(size=3,colour="black",alpha=0.6,shape=21)+
  scale_fill_manual(values=my_palet)+ theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
dev.off()

