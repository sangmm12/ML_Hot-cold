setwd("D:/R/ML_PAAD/data1/AUCAtcga/ESTIMATE/cluster_ESTIMATE")

library(data.table)
cluster <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/CIBERSORT/clusterresult/clusterresult.csv",sep=''),header =F)
cluster <- as.data.frame(cluster)
rownames(cluster) <- cluster[,1]


data <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/ESTIMATE/3data_immueScore.csv",sep=''),header = T)
data <- as.data.frame(data)

#library(dplyr)
#data<- data %>% distinct(ID,B_cells_naive, .keep_all = T)

rownames(data) <- data[,1]
dat_immu <- data[,-1]


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
pdf(paste("3datacluster_ESTIMATE.pdf",sep=''),width=4,height = 6)
p <- ggboxplot(dat, x = "Gene", y = "value",
               color = "Group", palette = c("#f8ac8c","#2878b5"), 
               add = "jitter",x.text.angle=60) # palette可以按照期刊选择相应的配色，如"npg"等
p <- p+xlab("Gene")+ylab("Expression Value")
p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
p
#print(p)
dev.off()



cells <- c("T_cells_CD8","NK_cells_activated","Macrophages_M2")

out_put_list <- c()

for(cell in unique(dat$Gene)){
  
  dat1 <- subset(dat,dat$Gene==cell)
  
  dat1$Group <- factor(dat1$Group,levels=c("C1","C2"))
  library(ggpubr)
  compare_means(value ~ Group, data = dat1, group.by = "Gene",method = "anova")
  
  
  pdf(paste("3datacluster_",cell,".pdf",sep=''),width=2.4,height = 6)
  p <- ggboxplot(dat1, x = "Group",y = "value",
                 color = "Group", palette = c("#f8ac8c","#2878b5"), 
                 add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
  p <- p+ylab(cell)
  #p <- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
  #p <- p+ theme(axis.text = element_text(size = 15),axis.title=element_text(size=30))
  print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
  #print(p)
  
  eval(parse(text = paste(cell,'<- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")')))
  out_put_list <- append(out_put_list,cell)
  
  
  dev.off()
}

pdf(paste("3datacluster_immuzong.pdf",sep=''),width = 8,height = 8*length(out_put_list)/4)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))
#eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))

dev.off()

