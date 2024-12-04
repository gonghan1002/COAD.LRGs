#############       富集分析    #############
#### 

library(data.table)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

rm(list=ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "2.Cluster/"

# 读取差异基因
degs <- read.csv("./res/2.cluster/table/degs.Cluster1.vs.merge_Cluster2.merge.fc1.5.csv")
degs_list <- degs[degs$threshold != "NoSig",]

# id convert
gene_id <- id_convert(degs_list[,1],"SYMBOL","ENTREZID")

#######################  富集分析 -- GO  #######################
kk <- enrichGO(gene = gene_id$ENTREZID,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, qvalueCutoff = 0.05,
               ont="BP", readable =T)
# 去重
# kk <- clusterProfiler::simplify(kk, cutoff=0.7, by="p.adjust", select_fun=min)
kk.df <- as.data.frame(kk)
# save
outfile_tmp <- paste0(res_home, proj_name,"/table/","GO_BP.","Cluster.","degs",".csv");outfile_tmp
write.csv(kk.df,file=outfile_tmp)

#### plot -1
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","GO.","cluster",".pdf");outfile_tmp
pdf(file = outfile_tmp,onefile = F,width = 8,height = 8)
barplot(kk,showCategory = 10)
dev.off()

#### plot -2
echo <- pairwise_termsim(kk)
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","GO.","emapplot.2",".pdf");outfile_tmp
pdf(file = outfile_tmp,onefile = F,width = 8,height = 5)
emapplot_cluster(echo)
dev.off()


#######################  富集分析 -- KEGG  #######################
kk <- enrichKEGG(gene = gene_id$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 0.25)
kk.df <- as.data.frame(kk)

# save
outfile_tmp <- paste0(res_home, proj_name,"/table/","KEGG.","Cluster.","degs",".csv");outfile_tmp
write.csv(kk.df,file=outfile_tmp)

# plot
# 不够十条 画个简单的气泡图
outfile_tmp <- paste0(res_home,proj_name,"pdf/","Kegg.","Cluster.",".pdf");outfile_tmp
pdf(file = outfile_tmp,onefile = F,width=6,height = 4)
barplot(kk,showCategory = 10)
dev.off()

#### plot -2
echo <- pairwise_termsim(kk)
# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","KEGG.","emapplot",".pdf");outfile_tmp
pdf(file = outfile_tmp,onefile = F,width = 8,height = 5)
emapplot(echo)
dev.off()
