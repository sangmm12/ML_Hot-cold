setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/mirna/cor_gene_mirna")

library(data.table)
dat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/mirna/Merge_miRNA_pre_Count.txt",sep=''))
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]
dat <- dat[,-1]
dat_gene <- t(dat)

mirna <- colnames(dat_gene)[which(colSums(dat_gene ==0) < nrow(dat_gene)*0.1)]
dat_gene <- dat_gene[,match(mirna,colnames(dat_gene))]

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
data.corr <- corr.test(log(dat_im+1), log(dat_gene+1), method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值

write.csv(data.r,file= paste("data.r_all0.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_all0.csv",sep=''),quote=F)

# 计算每列中NA值的数量
na_countsr <- colSums(is.na(data.r))
na_countsp <- colSums(is.na(data.p))
# 选择不全部是NA的列
data.r <- data.r[, na_countsr != nrow(data.r)]
data.p <- data.p[, na_countsp != nrow(data.p)]

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
myBreaks <- c(seq(-0.9, 0, length.out = floor(paletteLength/2) + 1))
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_gene_mirna_all.pdf",sep=''),width =100,height = 24)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()

#save(data.r,data.p, file = "data.r_p.rda")
# load("data.r_p.rda")



dat_r <- data.r
number <- rep(0, times=length(colnames(dat_r)))

for (i in 1:length(colnames(dat_r))) {
  aaa <- subset(dat_r, abs(dat_r[,i]) > 0.2)
  number[i] <- length(rownames(aaa))
  
}

numberr <- number

dat_r <- data.r
number <- rep(0, times=length(colnames(dat_r)))

for (i in 1:length(colnames(dat_r))) {
  aaa <- subset(dat_r, abs(dat_r[,i]) > 0.3)
  number[i] <- length(rownames(aaa))
  
}

numberr <- rbind(numberr,number)

dat_r <- data.r
number <- rep(0, times=length(colnames(dat_r)))

for (i in 1:length(colnames(dat_r))) {
  aaa <- subset(dat_r, abs(dat_r[,i]) > 0.4)
  number[i] <- length(rownames(aaa))
  
}

numberr <- rbind(numberr,number)



dat_r <- data.r
number <- rep(0, times=length(colnames(dat_r)))

for (i in 1:length(colnames(dat_r))) {
  aaa <- subset(dat_r, abs(dat_r[,i]) > 0.5)
  number[i] <- length(rownames(aaa))
  
}

numberr <- rbind(numberr,number)

dat_r <- data.r
number <- rep(0, times=length(colnames(dat_r)))

for (i in 1:length(colnames(dat_r))) {
  aaa <- subset(dat_r, abs(dat_r[,i]) > 0.6)
  number[i] <- length(rownames(aaa))
  
}

numberr <- rbind(numberr,number)
colnames(numberr) <- colnames(data.r)
rownames(numberr) <- c("0.2","0.3","0.4","0.5","0.6")

#numberr$Total <-rowSums(numberr[,c(1:length(colnames(numberr)))])

write.csv(numberr,file=paste("numberr.csv",sep=''),row.names = T)

out_put_list <- c()
mirna <- c()
gene <- c()
for (i in 1:length(rownames(data.r))) {
  for (j in 1:length(colnames(data.r))) {
    if(abs(data.r[i,j]) > 0.5){
      
      
      # which(colnames(data.r)=="hsa-mir-5583-2")
      # [1] 1017
      # which(rownames(data.r)=="AC018797.2")
      # [1] 65
      mirna <- c(mirna,colnames(data.r)[j])
      gene <- c(gene,rownames(data.r)[i])
      
      #colnames(dat_gene)
      dat_1 <- dat_im[,rownames(data.r)[i]]
      dat_2 <- dat_gene[, colnames(data.r)[j]]
      plot(dat_1,dat_2)
      
      
      library(ggpubr)
      library(reshape2)
      library(ggsci)
      d1 <- data.frame(NT5E=log(dat_1+1) ,VEGFA=log(dat_2+1))
      
      library(Hmisc)
      
      #ggscatterstats(data=d1,x=log(NT5E+1),y=log(VEGFA+1),bf.message = F)
      
      
      pdf( paste("0.5_log/cor_",rownames(data.r)[i],"_",colnames(data.r)[j],".pdf",sep=''),height=5,width=5)
      # p <- ggplot(d1, aes( x = NT5E, y =VEGFA)) + geom_point( size=3) + geom_smooth(method = 'lm', formula = y ~ x, se = T)+
      #   theme(axis.text = element_text(size = 25),axis.title = element_text(size = 25))+
      #   theme_classic() +labs(x="log(NT5E+1)",y="log(CD276+1)") + stat_cor(method = "pearson",size=8,,label.x =min(dat_1), label.y =max(dat_2))
      
      p <- ggplot(d1,aes(x=NT5E,y=VEGFA)) + geom_point(shape=19) +
        xlab( paste(rownames(data.r)[i],sep='')) + ylab(paste(colnames(data.r)[j],sep=''))+geom_smooth(method = lm)+theme_classic()+
        theme(axis.text = element_text(size = 15),axis.title = element_text(size = 25))+
        stat_cor(method = "pearson",size=8,label.x =min(log(dat_1+1)), label.y =max(log(dat_2+1)))
      
      print(p)
      eval(parse(text = paste('list',i,'_',j,' <- p',sep='')))
      out_put_list <- eval(parse(text = paste("append(out_put_list,'","list",i,"_",j,"')",sep=''))) 
      dev.off()
      
    }else {
      aa <- 1
    }
    
  }
}

pdf(paste("0.6_log/cor_miRNA_0.6.pdf",sep=''),width = 40/2,height = 10*length(out_put_list)/8)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
#eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))

dev.off()

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
myBreaks <- c(seq(-0.6, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.7/paletteLength, 0.7, length.out=floor(paletteLength/2)))
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_gene_mirna_0.5.pdf",sep=''),width =3.6,height = 7)
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
data.corr <- corr.test(dat_imhot_low, dat_genehot_low, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


write.csv(data.r,file= paste("data.r_hot_low.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_hot_low.csv",sep=''),quote=F)

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
myBreaks <- c(seq(-0.6, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.7/paletteLength, 0.7, length.out=floor(paletteLength/2)))
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_gene_mirna_0.5_hot_low.pdf",sep=''),width =3.6,height = 7)
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
data.corr <- corr.test(dat_imcold_high, dat_genecold_high, method="pearson", adjust="fdr")
data.r <- data.corr$r  # 相关系数
data.p <- data.corr$p  # p值


write.csv(data.r,file= paste("data.r_cold_high.csv",sep=''),quote=F)
write.csv(data.p,file= paste("data.p_cold_high.csv",sep=''),quote=F)

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
myBreaks <- c(seq(-0.6, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(0.7/paletteLength, 0.7, length.out=floor(paletteLength/2)))
# myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
#               seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
#myBreaks <- c(seq(-0.2, 0, length.out=ceiling(paletteLength/2) + 1), 
#seq(0.4/paletteLength, 0.4, length.out=floor(paletteLength/2)))


#chr [1:6, 1:5] "*" "***" "" "***" "***" "***" "***" "" "***" "**" ...

pdf(paste("cor_gene_mirna_0.5_cold_high.pdf",sep=''),width =3.6,height = 7)
pheatmap(data.r, 
         color=myColor,
         breaks=myBreaks,
         clustering_method="average", cluster_rows=F, cluster_cols=F,display_numbers=sig.mat
)

dev.off()



