setwd('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/pici/start')
library(stringr)

final_data = data.frame()

sample_name <- 'TCGA_ICGC'
my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name,'.csv',sep=''), header = T)
gene_exp <- as.data.frame(my_data)
gene_exp$SampleName = rep(str_to_upper(sample_name),length(gene_exp$SampleName))
gene_exp = gene_exp[,order(colnames(gene_exp))]
final_data = rbind(final_data,gene_exp)


sample_name <- 'GSE62452'
for(sample_name in c("GSE85916","GSE28735","GSE62452","GSE78229"))
{
  my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name,'.csv',sep=''), header = T)
  gene_exp <- as.data.frame(my_data)
  gene_exp$SampleName = rep(str_to_upper(sample_name),length(gene_exp$SampleName))
  gene_exp = gene_exp[,order(colnames(gene_exp))]
  all_name <- names(which(table(c(colnames(gene_exp),colnames(final_data)))==2))
  gene_exp <- gene_exp[,match(all_name,colnames(gene_exp))]
  final_data <- final_data[,match(all_name,colnames(final_data))]
  final_data = rbind(final_data,gene_exp)
}


final_data0 <- data.frame(lapply(final_data[,which(colnames(final_data)!='SampleName')], as.numeric))
com1 <- prcomp(final_data0, center = TRUE,scale. = TRUE)
summary(com1)
df1<-com1$x
head(df1)

df1<-data.frame(df1,final_data$SampleName)
head(df1)
library(ggplot2)
  summ<-summary(com1)
  xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
  ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
  pdf('out.pdf',width=8,height=7)
  p2<-ggplot(data = df1,aes(x=PC1,y=PC2,color=final_data.SampleName))+ 
    stat_ellipse(aes(fill=final_data.SampleName), type = "norm", geom ="polygon",alpha=0.2,color=NA)+ 
    geom_point()+labs(x=xlab,y=ylab,color="")+guides(fill=F)+
    scale_fill_manual(values = c('#ff4d4d','#d01257','#3bb4c1','#0f1021','#ffc300'))+ scale_colour_manual(values = c('#ff4d4d','#d01257','#3bb4c1','#0f1021','#ffc300'))+
    theme_bw()+
    guides(color=guide_legend(override.aes = list(size=5,alpha=1)))+
    theme(axis.title = element_text(size = 15,color='black'),axis.text = element_text(size = 20,color='black'))
  print(p2)
  dev.off()
