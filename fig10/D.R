setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/GSVA_hallmark/cor_drug_hallmark")

sample_name <- "TCGA_ICGC"

library(data.table)
dat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/GSVA_hallmark/t_results_hallmark.csv",sep=''))
dat <- as.data.frame(dat)
dat <- dat[which(dat$t_value >= 2.58),]
pathway <- dat$V1
dat <- readRDS(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/GSVA_hallmark/hallmark.rds",sep=''))
dat <- as.data.frame(dat)
dat <- dat[match(pathway,rownames(dat)),]

dat_gene <- t(dat)
colnames(dat_gene) <- sub("HALLMARK_", "", colnames(dat_gene))


dat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/out_put_",sample_name,"_DEGs_down_zong.csv",sep=''))
dat <- as.data.frame(dat)

rownames(dat) <- dat[,1]
dat <- dat[,-1]

datar <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/zong/cor_gene/data.r2.csv",sep=''))
drug <- colnames(datar)[-1]

dat_im <- dat[,match(drug,colnames(dat))]

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene) ))==2))


dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]

# for(i in 1:length(colnames(dat_gene)))
# {
#   dat_gene[,i] = as.numeric(unlist(dat_gene[,i]))
# }
# i=1
# nrow(dat_gene)
# ncol(dat_gene)
# 

colSums(dat_im)


library(psych)
data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


write.csv(data.r,file= paste("data.r_all.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_all.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(0.9/paletteLength, 0.9, length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_drug_hallmark.pdf",sep=''),width =12,height = 16)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()








#####################################
Sel <- 'hot-low'

rf_risk <- read.csv('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/clusterhc_riskscore/KM4/sur_clusterhc_riskscore_KM4.csv')
rf_risk <- rf_risk[,c(1,7)]
colnames(rf_risk) <- c("SampleName","risk")
hot_low <- rf_risk[which(rf_risk$risk==Sel),]$SampleName

Sel <- 'cold-high'

rf_risk <- read.csv('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/clusterhc_riskscore/KM4/sur_clusterhc_riskscore_KM4.csv')
rf_risk <- rf_risk[,c(1,7)]
colnames(rf_risk) <- c("SampleName","risk")
cold_high <- rf_risk[which(rf_risk$risk==Sel),]$SampleName





################################################
dat_genehot_low <- dat_gene[match(hot_low,rownames(dat_gene)),]
dat_imhot_low <- dat_im[match(hot_low,rownames(dat_im)),]


library(psych)
data.corr <- corr.test(dat_genehot_low,dat_imhot_low,  method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


write.csv(data.r,file= paste("data.r_hot_low.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_hot_low.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(0.9/paletteLength, 0.9, length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_drug_hallmark_hot_low.pdf",sep=''),width =12,height = 16)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()





#######################################
dat_genecold_high <- dat_gene[match(cold_high,rownames(dat_gene)),]
dat_imcold_high <- dat_im[match(cold_high,rownames(dat_im)),]


library(psych)
data.corr <- corr.test( dat_genecold_high, dat_imcold_high, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


write.csv(data.r,file= paste("data.r_cold_high.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_cold_high.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(0.9/paletteLength, 0.9, length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_drug_hallmark_cold_high.pdf",sep=''),width =12,height = 16)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()










###########################################################
setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/GSVA_hallmark/cor_drug_hallmark/10")

sample_name <- "TCGA_ICGC"

library(data.table)
dat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/GSVA_hallmark/t_results_hallmark.csv",sep=''))
dat <- as.data.frame(dat)
dat <- dat[which(dat$t_value >= 2.58),]
pathway <- dat$V1
dat <- readRDS(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/GSVA_hallmark/hallmark.rds",sep=''))
dat <- as.data.frame(dat)
dat <- dat[match(pathway,rownames(dat)),]

dat_gene <- t(dat)
colnames(dat_gene) <- sub("HALLMARK_", "", colnames(dat_gene))


dat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/out_put_",sample_name,"_DEGs_down_zong.csv",sep=''))
dat <- as.data.frame(dat)

rownames(dat) <- dat[,1]
dat <- dat[,-1]

datar <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/zong/cor_gene/data.r2.csv",sep=''))
drug <- colnames(datar)[-1]

dat_im <- dat[,match(drug,colnames(dat))]

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene) ))==2))


dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]

# for(i in 1:length(colnames(dat_gene)))
# {
#   dat_gene[,i] = as.numeric(unlist(dat_gene[,i]))
# }
# i=1
# nrow(dat_gene)
# ncol(dat_gene)
# 

colSums(dat_im)


library(psych)
data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

data.r <- data.r[1:10,]
data.p <- data.p[1:10,]

write.csv(data.r,file= paste("data.r_all.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_all.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(0.9/paletteLength, 0.9, length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_drug_hallmark.pdf",sep=''),width =8,height = 5)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()








#####################################
Sel <- 'hot-low'

rf_risk <- read.csv('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/clusterhc_riskscore/KM4/sur_clusterhc_riskscore_KM4.csv')
rf_risk <- rf_risk[,c(1,7)]
colnames(rf_risk) <- c("SampleName","risk")
hot_low <- rf_risk[which(rf_risk$risk==Sel),]$SampleName

Sel <- 'cold-high'

rf_risk <- read.csv('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/clusterhc_riskscore/KM4/sur_clusterhc_riskscore_KM4.csv')
rf_risk <- rf_risk[,c(1,7)]
colnames(rf_risk) <- c("SampleName","risk")
cold_high <- rf_risk[which(rf_risk$risk==Sel),]$SampleName





################################################
dat_genehot_low <- dat_gene[match(hot_low,rownames(dat_gene)),]
dat_imhot_low <- dat_im[match(hot_low,rownames(dat_im)),]


library(psych)
data.corr <- corr.test(dat_genehot_low,dat_imhot_low,  method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

data.r <- data.r[1:10,]
data.p <- data.p[1:10,]

write.csv(data.r,file= paste("data.r_hot_low.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_hot_low.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(0.9/paletteLength, 0.9, length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_drug_hallmark_hot_low.pdf",sep=''),width =8,height = 5)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()





#######################################
dat_genecold_high <- dat_gene[match(cold_high,rownames(dat_gene)),]
dat_imcold_high <- dat_im[match(cold_high,rownames(dat_im)),]


library(psych)
data.corr <- corr.test( dat_genecold_high, dat_imcold_high, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


data.r <- data.r[1:10,]
data.p <- data.p[1:10,]

write.csv(data.r,file= paste("data.r_cold_high.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_cold_high.csv",sep=''),quote=F)


library(pheatmap)
getSig <- function(dc) {
  print(dc)
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''
  }
  return(sc)
}
sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
str(sig.mat)

paletteLength <- 1000
myColor <- colorRampPalette(c("#36648b", "white", "#e94644"))(paletteLength)

test <- data.r
myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(0.9/paletteLength, 0.9, length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_drug_hallmark_cold_high.pdf",sep=''),width =8,height = 5)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()

