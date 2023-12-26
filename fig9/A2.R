setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/mirna/cor_gene_mirna")

load("mirna_gene_all.rda")
load("mirna_gene_hot_low.rda")
load("mirna_gene_cold_high.rda")

gene0 <- c(unique(gene_all),unique(gene_hot_low),unique(gene_cold_high))
gene <- unique(gene0)

mirna0 <- c(unique(mirna_all),unique(mirna_hot_low),unique(mirna_cold_high))
mirna <- unique(mirna0)

load("data.r_p.rda")

data.r <- data.r[,match(unique(mirna),colnames(data.r))]
data.p <- data.p[,match(unique(mirna),colnames(data.p))]

data.r <- data.r[match(unique(gene),rownames(data.r)),]
data.p <- data.p[match(unique(gene),rownames(data.p)),]

data.r <- t(data.r)
data.p <- t(data.p)

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
#myBreaks <- c(seq(-0.9, 0, length.out = floor(paletteLength/2) + 1))
myBreaks <- c(seq(-0.7, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(1/paletteLength,1, length.out=floor(paletteLength/2)))
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_gene_mirna_0.7_3.pdf",sep=''),width =5,height = 7)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()




load("data.r_p_hot_low.rda")

data.r <- data.r[,match(unique(mirna),colnames(data.r))]
data.p <- data.p[,match(unique(mirna),colnames(data.p))]

data.r <- data.r[match(unique(gene),rownames(data.r)),]
data.p <- data.p[match(unique(gene),rownames(data.p)),]

data.r <- t(data.r)
data.p <- t(data.p)

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
#myBreaks <- c(seq(-0.9, 0, length.out = floor(paletteLength/2) + 1))
myBreaks <- c(seq(-0.7, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(1/paletteLength,1, length.out=floor(paletteLength/2)))
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_gene_mirna_0.7_3_hot_low.pdf",sep=''),width =5,height = 7)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()





load("data.r_p_cold_high.rda")

data.r <- data.r[,match(unique(mirna),colnames(data.r))]
data.p <- data.p[,match(unique(mirna),colnames(data.p))]

data.r <- data.r[match(unique(gene),rownames(data.r)),]
data.p <- data.p[match(unique(gene),rownames(data.p)),]

data.r <- t(data.r)
data.p <- t(data.p)

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
#myBreaks <- c(seq(-0.9, 0, length.out = floor(paletteLength/2) + 1))
myBreaks <- c(seq(-0.7, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(1/paletteLength,1, length.out=floor(paletteLength/2)))
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_gene_mirna_0.7_3_cold_high.pdf",sep=''),width =5,height = 7)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()


