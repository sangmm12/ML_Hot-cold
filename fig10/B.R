setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/zong/cluster_drug")



library(data.table)
rf_risk <- read.csv('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/clusterhc_riskscore/KM4/sur_clusterhc_riskscore_KM4.csv')
rf_risk <- rf_risk[,c(1,7)]
cluster <- subset(rf_risk, group %in% c("hot-low", "cold-high"))
cluster <- as.data.frame(cluster)
rownames(cluster) <- cluster[,1]
colnames(cluster) <- c("V1","V2") 


sample_name <- "TCGA_ICGC"
library(data.table)
dat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/out_put_",sample_name,"_DEGs_down_zong.csv",sep=''))
dat <- as.data.frame(dat)


rownames(dat) <- dat[,1]
dat<- dat[,-1]

#run一下DEGs_downRMA_zong_cor_gene.R
# drug <- c("BI-2536","Trametinib","Dasatinib","SB505124","SCH772984","ERK","NU7441","JQ1","RO-3306","ZM447439","Axitinib","Tozasertib",
#           "Selumetinib","KU-55933","SB216763")

datar <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/zong/cor_gene/data.r2.csv",sep=''))
drug <- colnames(datar)[-1]

dat <- dat[,match(drug,colnames(dat))]


dat_immu <- dat


all_name <- names(which(table(c(rownames(dat_immu),cluster[,1] ))==2))


dat_cluster <- cluster[match(all_name,cluster[,1]),]

dat_imm <- dat_immu[match(all_name,rownames(dat_immu)),]
cell <- colnames(dat_imm)
#cell <- c("B_cells_naive","B_cells_memory","Plasma_cells","T_cells_CD8","T_cells_CD4_naive","T_cells_CD4_memory_resting","T_cells_CD4_memory_activated","T_cells_follicular_helper","T_cells_regulatory_.Tregs.","T_cells_gamma_delta","NK_cells_resting","NK_cells_activated","Monocytes","Macrophages_M0","Macrophages_M1","Macrophages_M2","Dendritic_cells_resting","Dendritic_cells_activated","Mast_cells_resting","Mast_cells_activated","Eosinophils","Neutrophils")
dat_im <- dat_imm[cell]

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


pdf(paste("PAAD_DEGs_down_cluster_Drug2_zong_hlch.pdf",sep=''),width=40,height = 8)
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
  
  
  if(cell=="968"){
    cell1 <- "drug_968"
    
  }else {
    cell1 <- gsub("-", "_", cell)
    cell1 <- gsub("[^A-Za-z0-9_]", "_", cell1)
  }
  
  library(ggpubr)
  compare_means(value ~ Group, data = dat1, group.by = "Gene",method = "anova")
  
  #  pdf(paste("DEGs_down_cluster_",cell1,".pdf",sep=''),width=2.7,height = 6)
  p <- ggboxplot(dat1, x = "Group", y = "value",
                 color = "Group", palette = c("FA7F6F"),
                 add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+xlab(" ")+ylab(paste(cell," (IC50)",sep=''))
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


pdf(paste("drug_cluster_down2_hlch.pdf",sep=''),width = 3*6,height = 5*length(out_put_list)/6)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=6))',sep='')))
#eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))

dev.off()



