setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/zong/cor_gene")

sample_name <- "TCGA_ICGC"

library(data.table)
dat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/out_put_",sample_name,"_DEGs_down_zong.csv",sep=''))
dat <- as.data.frame(dat)


rownames(dat) <- dat[,1]
dat_gene <- dat[,-1]


is.numeric(dat_gene)

genedat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/KM/single_coxTCGA_ICGC_down_p_HR.csv",sep=''),header = T)
genedown <- genedat$genename
genedat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/KM/single_coxTCGA_ICGC_up_p_HR.csv",sep=''),header = T)
geneup <- genedat$genename
geneupdown <- c(geneup,genedown)

sample_name <- "TCGA_ICGC"

data <- fread(paste('D:/R/ML_PAAD/data1/alldata/exp_',sample_name,'.csv',sep=''), header = T)
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- data[,-1]

dat_im <- data[match(genedown,rownames(data)),]
dat_im <- t(dat_im)

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
data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
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
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_",sample_name,"_DEGs_down.pdf",sep=''),width =100,height = 32)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()


drug <- colnames(data.r)[which(colSums(data.r < 0)==length(rownames(data.r)))]


# drug <- c("Daporinad","Vinblastine","Vinblastine","Vinorelbine","Eg5","Dinaciclib","Paclitaxel","CDK9","Vincristine","Podophyllotoxin bromide",
#           "PD0325901","MK-1775","Trametinib","Dihydrorotenone","BMS-754807","AZD5153","Osimertinib","Afatinib","Wee1 Inhibito","Taselisib","AZD5438",
#           "Ulixertinib","Entinostat","YK-4-279","ULK1","AZD5582","VSP34","PAK","OTX015","Afuresertib","Erlotinib","ERK","SCH772984","AZD3759","Sorafenib",
#           "IWP-2","NU7441","JQ1","GSK343","VX-11e","Fulvestran","GSK269962A","Lapatinib","ZM447439","MK-2206","Crizotinib","RO-3306","Gefitinib","AZD6482",
#           "PRT062607","I-BET-762","BMS-345541","VE-822","Elephantin","Ipatasertib","Sinularin","Alpelisib","Tamoxifen","Nilotinib","Palbociclib",
#           "WIKI4","GSK2606414","Oxaliplatin","Linsitinib","TAF1","PF−4708671","Ribociclib","Sapitinib","LGK974","OF-1","VE821","AGI6780",
#           "Wnt-C59","Selumetinib","GSK1904529A","AZD5991","KRAS(G12C)Inhibitor−12","ML323","")


data.r <- data.r[,match(drug,colnames(data.r))]
data.p <- data.p[,match(drug,colnames(data.p))]

write.csv(data.r,file= paste("data.r1.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p1.csv",sep=''),quote=F)


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
myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_",sample_name,"_DEGs_down1.pdf",sep=''),width =35,height = 32)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()















drug <- colnames(data.p)[which(colSums(data.p < 0.0001)==length(rownames(data.p)))]

data.r <- data.r[,match(drug,colnames(data.r))]
data.p <- data.p[,match(drug,colnames(data.p))]

write.csv(data.r,file= paste("data.r2.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p2.csv",sep=''),quote=F)


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
myColor <- colorRampPalette(c("#36648b", "white"))(paletteLength)

test <- data.r
myBreaks <- c(seq(-1, 0, length.out = floor(paletteLength/2) + 1))
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1), 
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_",sample_name,"_DEGs_down2.pdf",sep=''),width =10,height = 32)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()


