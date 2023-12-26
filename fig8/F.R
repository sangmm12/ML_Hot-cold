setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/cor_gene_cell")

library(data.table)
CIBER <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/CIBERSORT/3data.csv",sep=''),header = T)
CIBER <- t(CIBER)
CIBER <- as.data.frame(CIBER)
colnames(CIBER) <- CIBER[1,]
CIBER <- CIBER[-1,]
dat_im <- CIBER
dat_im[] <- lapply(dat_im, as.numeric)

sample_name <- "TCGA_ICGC"
data<- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name,'_HRp.csv',sep=''), header = T)
data <- as.data.frame(data)
rownames(data) <- data[,1]
dat_gene <- data[,-1]

all_name <- names(which(table(c(rownames(dat_im),rownames(dat_gene) ))==2))


dat_gene <- dat_gene[match(all_name,rownames(dat_gene)),]
dat_im <- dat_im[match(all_name,rownames(dat_im)),]


colSums(dat_im)


library(psych)
data.corr <- corr.test( dat_gene, dat_im,method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


write.csv(data.r,file= paste("data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p.csv",sep=''),quote=F)


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
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.5, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.6/paletteLength, 0.6, length.out=floor(paletteLength/2)))



#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_genedown_cell.pdf",sep=''),width =7,height = 22)
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
data.corr <- corr.test( dat_genehot_low, dat_imhot_low,method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


write.csv(data.r,file= paste("data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p.csv",sep=''),quote=F)


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
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.5, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.6/paletteLength, 0.6, length.out=floor(paletteLength/2)))



#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_genedown_hot_low_cell.pdf",sep=''),width =7,height = 22)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()






######################################
dat_genecold_high <- dat_gene[match(cold_high,rownames(dat_gene)),]
dat_imcold_high <- dat_im[match(cold_high,rownames(dat_im)),]


library(psych)
data.corr <- corr.test( dat_genecold_high, dat_imcold_high,method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


write.csv(data.r,file= paste("data.r.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p.csv",sep=''),quote=F)


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
# myBreaks <- c(
#   seq(-0.7, 0, length.out = floor(paletteLength/2) + 1)
# )
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
myBreaks <- c(seq(-0.5, 0, length.out=ceiling(paletteLength/2) + 1),
seq(0.6/paletteLength, 0.6, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_genedown_cold_high_cell.pdf",sep=''),width =7,height = 22)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()
