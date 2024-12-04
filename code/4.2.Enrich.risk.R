
library(data.table)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

rm(list = ls())
# 导入函数
libSources <- list.files("E:/BaiduSyncdisk/r_project/CodeLib/", recursive=TRUE, full.names=TRUE, pattern="\\.R$")
for(i in 1:length(libSources))
  source(libSources[i])

# 设定输出目录与项目名字
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/res/"
proj_name <- "4.Biology/"

# 读取差异基因
degs <- read.csv("./res/4.Biology/table/degs.limma.fc_2.csv")
degs_list <- degs[degs$threshold != "NoSig",]
gene_id <- id_convert(degs_list[,1],"SYMBOL","ENTREZID")

### go
kk <- enrichGO(gene = gene_id$ENTREZID,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, qvalueCutoff = 0.05,
               ont="BP", readable =T)
kk.df <- as.data.frame(kk)
outfile_tmp <- paste0(res_home, proj_name,"/table/","GO_BP.","risk.","degs",".csv");outfile_tmp
write.csv(kk.df,file=outfile_tmp)

#### plot
colnames(kk.df)
data_pre <- kk.df[,c(2,6)]
data_pre$score <- -log10(data_pre$p.adjust)
colnames(data_pre)[2] <- c("term")
data_pre$term <- factor(data_pre$term,levels=rev(data_pre$term))

# plot
my_plat <- colorRampPalette(RColorBrewer::brewer.pal(5,"Reds"),space="rgb")(50)

# save
outfile_tmp <- paste0(res_home,proj_name,"pdf/","GO-BP.","risk.","top10",".pdf");outfile_tmp
pdf(file = outfile_tmp,onefile = F,width= 8,height = 6)
ggplot(data = data_pre, aes(x = score, y = term)) +
  geom_bar(stat="identity",position = "stack",width = 0.8,
           fill = "#F26553") +
  #  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x=element_text(color="black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(color="black",size = 12),
        axis.title.y = element_blank(), 
        panel.background = element_rect(color="black",fill=NA),
        panel.grid=element_blank(),
        legend.title=element_text(size=8), 
        legend.text = element_text(size=10)) +
  labs(title = " ") +
  xlab(expression(paste(-Log[10],"(p.adjust)",sep = "")))+
  scale_x_continuous(expand=c(0,0),limits = c(0,23))+
  geom_vline(xintercept = -log10(0.05),lty=2,col="#D2B509",lwd=1) # 添加0.05显著线
dev.off()

# kegg
kk <- enrichKEGG(gene = gene_id$ENTREZID,
                 # organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 0.05)
kk.df <- as.data.frame(kk)

# save
outfile_tmp <- paste0(res_home, proj_name,"/table/","KEGG.","risk.","degs",".csv");outfile_tmp
write.csv(kk.df,file=outfile_tmp)

#### 数据处理
colnames(kk.df)
data_pre <- kk.df[,c(2,6)]
data_pre$score <- -log10(data_pre$p.adjust)
colnames(data_pre)[1] <- c("term")

# order
data_pre$term <- factor(data_pre$term,levels=rev(data_pre$term))

my_plat <- colorRampPalette(RColorBrewer::brewer.pal(5,"Blues"),space="rgb")(50)
outfile_tmp <- paste0(res_home,proj_name,"pdf/","KEGG.","risk.","top10",".pdf");outfile_tmp
pdf(file = outfile_tmp,onefile = F,width= 8,height = 6)
ggplot(data = data_pre, aes(x = score, y = term)) +
  geom_bar(stat="identity",position = "stack",width = 0.8,
           fill = my_plat[38]) +  # 代码 #3081BC
  #  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x=element_text(color="black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 12), # 调整标签
        axis.text.y = element_text(color="black",size = 12), # 调整坐标刻度
        axis.title.y = element_blank(), # 调整标签
        panel.background = element_rect(color="black",fill=NA),
        panel.grid=element_blank(),
        legend.title=element_text(size=8), 
        legend.text = element_text(size=10)) +
  labs(title = " ") +
  xlab(expression(paste(-Log[10],"(p value)",sep = "")))+
  scale_x_continuous(expand=c(0,0),limits = c(0,8))+
  geom_vline(xintercept = -log10(0.05),lty=2,col="#D2B509",lwd=1) # 添加0.05显著线
dev.off()

#### GSEA
gene_id <- id_convert(input = degs[,1],from = "SYMBOL",to = "ENTREZID")
data_pre <- merge(gene_id, degs, by.x = "SYMBOL", by.y = "X")
colnames(data_pre)
logFC.list <- data_pre$logFC
names(logFC.list) <- data_pre$ENTREZID
logFC.list = sort(logFC.list,decreasing = T)

TERM2GENE_set <- read.gmt("E:/r_project/data/Msigdb/h.all.v2022.1.Hs.entrez.gmt")
colnames(TERM2GENE_set) <- c("ont", "gene")
kk <- GSEA(logFC.list, TERM2GENE = TERM2GENE_set, 
           exponent=1, minGSSize=5,  
           pvalueCutoff=0.05, pAdjustMethod="BH",
           seed=F, by="fgsea")
kk.df <- as.data.frame(kk)

library(GseaVis)
setid <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_INFLAMMATORY_RESPONSE",
           "HALLMARK_INTERFERON_GAMMA_RESPONSE")
outfile_tmp <- paste0(res_home, proj_name,"/pdf/","GSEA.","junjun",".pdf");outfile_tmp
# save
pdf(file = outfile_tmp, onefile = F,width = 6,height = 6)
gseaNb(object = kk,
       geneSetID = setid,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.8,
       pvalY = 0.5,
       pFill = "white")
dev.off()
