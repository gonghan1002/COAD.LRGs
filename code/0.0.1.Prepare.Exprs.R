

library(data.table)
library(tidyverse)

# 设定输出目录与项目名字
getwd()
res_home <- "E:/BaiduSyncdisk/r_project/COAD.LRGs/"

checkDir(paste0(res_home,"/code/"))
checkDir(paste0(res_home,"/data/"))
checkDir(paste0(res_home,"/list/"))
checkDir(paste0(res_home,"/res/"))

# read data
data_exprs <- readRDS("../data/TCGA_data/TPM/TCGA-COAD.TPM.rds")
head(data_exprs[,1:4])
range(data_exprs)

# data anno
data_gene <- fread("../data/TCGA_data/TCGA-mRNA.gene_anno.csv",data.table = F)
head(data_gene[,1:4])

### sub
data_gene_pre <- data_gene[data_gene$gene_type %in% c("protein_coding",
                                                      "nonsense_mediated_decay",
                                                      "non_stop_decay",
                                                      "IG_C_gene","IG_D_gene",
                                                      "IG_J_gene","IG_V_gene",
                                                      "polymorphic_pseudogene",
                                                      "protein_coding_LoF",
                                                      "TR_C_gene","TR_D_gene",
                                                      "TR_J_gene","TR_V_gene"),]
data_gene_pre <- separate(data_gene_pre,V1,into = c("V1"),sep="\\.")

data_exprs <- as.data.frame(data_exprs)
data_exprs <- rownames_to_column(data_exprs,"Ensembl_ID")
data_exprs <- separate(data_exprs,Ensembl_ID,into = c("Ensembl_ID"),sep="\\.") 
data_exprs <- data_exprs[data_exprs[,1] %in% data_gene_pre$V1,]

### id convert
colnames(data_exprs)[1] <- "ENSEMBL"
gene_id <- id_convert(data_exprs[,1],"ENSEMBL","SYMBOL")
# merge 
data_exprs <- merge(gene_id, data_exprs, "ENSEMBL")
# duplicated
data_exprs <- data_exprs[!duplicated(data_exprs[,2]),]
data_exprs <- data.frame(data_exprs[,-c(1,2)], row.names = data_exprs[,2])

range(data_exprs)
data_exprs <- log2(data_exprs + 1)

#### save
table(substring(colnames(data_exprs),first = 13,last = 16))#观察样本类型
tumor_idx <- grep("01A",substring(colnames(data_exprs),first = 13,last = 16))
normol_idx <- grep("11A",substring(colnames(data_exprs),first = 13,last = 16))
data_pre <- data_exprs[,c(tumor_idx, normol_idx)]
head(data_pre[,1:4])
# save
saveRDS(data_pre, file = "./data/TCGA-COAD.TPM.RDS")

# screen tumor
data_pre <- data_exprs[,tumor_idx]
saveRDS(data_pre, file = "./data/TCGA-COAD.TPM.Tumor.RDS")

