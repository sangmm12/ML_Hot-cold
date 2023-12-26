setwd('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/pici/highlow')
library(stringr)
library(ggplot2)
library(ggside)
library(data.table)

sample_name <- 'TCGA_ICGC'
sample_name <- 'GSE62452'
for(sample_name in c('TCGA_ICGC',"GSE85916","GSE28735","GSE62452","GSE78229"))
{
  my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name,'.csv',sep=''), header = T)
  gene_exp <- as.data.frame(my_data)
  riskscore <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/queue/risk/",sample_name,"/risk.csv",sep=''),header = T)
  risk <- as.data.frame(riskscore)
  colnames(risk)[1] = 'SampleName'
  final_data = merge(risk,gene_exp,by='SampleName')[,-c(1,2,4,5)]
  colnames(final_data)[1] = 'SampleName'
  
  #prcomp() 函数来执行主成分分析
  final_data0 <- data.frame(lapply(final_data[,which(colnames(final_data)!='SampleName')], as.numeric))
  com1 <- prcomp(final_data0, center = TRUE,scale. = TRUE)
  summary(com1)
  df1<-com1$x
  head(df1)
  
  df1<-data.frame(df1,final_data$SampleName)
  head(df1)
  #c('#ff4d4d','#d01257','#3bb4c1','#0f1021','#ffc300')
  summ<-summary(com1)
  xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
  ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
  pdf(paste0(sample_name,'.pdf'),width=7,height=7)
  p2<-ggplot(data = df1,aes(x=PC1,y=PC2,color=final_data.SampleName))+ 
    stat_ellipse(aes(fill=final_data.SampleName), type = "norm", geom ="polygon",alpha=0.2,color=NA)+ 
    geom_point()+labs(x=xlab,y=ylab,color="")+guides(fill = "none")+
    scale_fill_manual(values = c('darkred','darkgreen'))+ scale_colour_manual(values = c('darkred','darkgreen'))+
    geom_xsidedensity(aes(fill=final_data.SampleName),show.legend = FALSE)+
    geom_ysidedensity(aes(fill=final_data.SampleName),show.legend = FALSE)+
    theme_bw()+
    guides(color=guide_legend(override.aes = list(size=5,alpha=1)))+
    theme(legend.position = c(0.07,0.833),axis.title = element_text(size = 25,color='black'),axis.text = element_text(size = 25,color='black'),ggside.panel.scale=0.1, ggside.axis.text=element_blank(),ggside.panel.background=element_blank())
  print(p2)
  dev.off()
}
