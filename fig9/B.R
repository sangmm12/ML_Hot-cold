setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/mirna/cluster_mirna")

library(data.table)
rf_risk <- read.csv('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/clusterhc_riskscore/KM4/sur_clusterhc_riskscore_KM4.csv')
rf_risk <- rf_risk[,c(1,7)]
cluster <- subset(rf_risk, group %in% c("hot-low", "cold-high"))
cluster <- as.data.frame(cluster)
rownames(cluster) <- cluster[,1]
colnames(cluster) <- c("V1","V2")   
    

library(data.table)
dat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/mirna/Merge_miRNA_pre_Count.txt",sep=''))
dat <- as.data.frame(dat)
rownames(dat) <- dat[,1]
dat <- dat[,-1]
dat <- t(dat)

#run一下DEGs_downRMA_zong_cor_gene.R
# drug <- c("BI-2536","Trametinib","Dasatinib","SB505124","SCH772984","ERK","NU7441","JQ1","RO-3306","ZM447439","Axitinib","Tozasertib",
#           "Selumetinib","KU-55933","SB216763")

load("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/mirna/cor_gene_mirna/mirna_gene_hot_low.rda")
load("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/mirna/cor_gene_mirna/mirna_gene_cold_high.rda")

mirna <- c(mirna_hot_low,mirna_cold_high)
mirna <- unique(mirna)

dat <- dat[,match(mirna,colnames(dat))]


dat_immu <- dat


all_name <- names(which(table(c(rownames(dat_immu),cluster[,1] ))==2))


dat_cluster <- cluster[match(all_name,cluster[,1]),]

dat_im <- dat_immu[match(all_name,rownames(dat_immu)),]
dat_im <- as.data.frame(dat_im)

dat_cluster <- dat_cluster[c('V2')]

dat <- data.frame()
for(coln in colnames(dat_im))
{
  for(i in all_name)
  {
    dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],dat_im[c(coln)][match(i,rownames(dat_im)),]))
  }
}

dat[,3] = as.numeric(dat[,3])

colnames(dat) <- c("Gene","Group","value")


library(ggpubr)
compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")


pdf(paste("cluster_mirna.pdf",sep=''),width=40,height = 8)
p <- ggboxplot(dat, x = "Gene", y = "value",
               color = "Group", palette = c("#2ecc71","#e74c3c"), 
               add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
p <- p+xlab("Drug")+ylab("Expression Value")
p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
p
#print(p)
dev.off()

#}

out_put_list <- c()

for(cell in unique(dat$Gene)){
  
  dat1 <- subset(dat,dat$Gene==cell)
  
  dat1$Group <- factor(dat1$Group,levels=c("hot-low","cold-high"))
  
  cell1 <- gsub("-", "_", cell)
  dat1$Gene <- gsub("-", "_", dat1$Gene)
  
  library(ggpubr)
  compare_means(value ~ Group, data = dat1, group.by = "Gene",method = "anova")
  
  #  pdf(paste("DEGs_down_cluster_",cell1,".pdf",sep=''),width=2.7,height = 6)
  p <- ggboxplot(dat1, x = "Group", y = "value",
                 color = "Group", palette = c("FA7F6F"),
                 add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab(" ")+ylab(paste(cell,sep=''))
  p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  p <- p + theme(axis.text = element_text(size = 15),axis.title.y = element_text(size = 15))
  # p <- p + scale_x_discrete(limits=c("hot","cold")
  # p
  # 
  # print(p)
  eval(parse(text = paste(cell1,' <- p',sep='')))
  out_put_list <- append(out_put_list,cell1)
  
  
  #dev.off()
  
}


pdf(paste("cluster_mirna_fen.pdf",sep=''),width = 40/2,height = 4*length(out_put_list)/7)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=7))',sep='')))
#eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))

dev.off()


