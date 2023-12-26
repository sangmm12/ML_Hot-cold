library(pROC)
library(stringr)
library(data.table)
setwd('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/pROC/KM4')
# 假设你有两个向量，分别存储了正常组和疾病组中该基因的表达量

#gene_name = 'GLP2R'

sample_name <- "TCGA_ICGC"

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


genedat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/KM/single_coxTCGA_ICGC_down_p_HR.csv",sep=''),header = T)
genedown <- genedat$genename
genedat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/KM/single_coxTCGA_ICGC_down_p_HR.csv",sep=''),header = T)
geneup <- genedat$genename
geneupdown <- c(geneup,genedown)

genedown <- gsub("-", "_", genedown)
gene_names <- genedown
#gene_names <- c("CYP2U1","KCTD12","AP003486.1","ZBED3","IGIP","GGTA1P","SPARCL1","WAC-AS1","LINC00847","CCR2")


auc_normal_cold_high <- c()
out_put_list <- c()

gene_name <- "AP003486.1"
for(gene_name in gene_names){
  
  dat <- fread(paste('D:/R/ML_PAAD/data1/alldata/exp_',sample_name,'.csv',sep=''), header = T)
  dat$V1 <- gsub("-", "_", dat$V1)
  temp_exp <- as.data.frame(dat)
  temp_exp <- temp_exp[which(temp_exp[,c(1)]==gene_name),]
  temp_exp <- temp_exp[,-1]
  expc0 <- temp_exp[,match(hot_low,colnames(temp_exp))]
  exp_hot_low <- as.numeric(as.vector(t(expc0)))
  
  expc0 <- temp_exp[,match(cold_high,colnames(temp_exp))]
  exp_cold_high <- as.numeric(as.vector(t(expc0)))
  
  normal_expression <- exp_cold_high
  disease_expression <- exp_hot_low
  
  # 合并两组表达量数据，并创建对应的类别标签 (0表示正常组，1表示疾病组)
  all_expression <- c(normal_expression, disease_expression)
  labels <- factor(c(rep(0, length(normal_expression)), rep(1, length(disease_expression))))
  
  # 计算 ROC 曲线的真阳性率和假阳性率
  roc_curve <- roc(labels, all_expression)
  
  gene_name1 <- gsub("[\r\n]", "_", gene_name)
  gene_name1 <- gsub("\\.", "_", gene_name1)
  
  # 绘制 ROC 曲线
  pdf(paste('cold_high/normal_cold_high_',gene_name,'.pdf',sep=''))
  p <- plot(roc_curve, main = paste("ROC Curve for Gene Expression of ",gene_name,sep=''), col = "#FFC107", lwd = 4,print.auc = TRUE, auc.polygon = TRUE,auc.polygon.col = "#FFECB3")
  print(p)
  dev.off()
  
  #输入dev.off()直到报错
  plot(roc_curve, main = paste("",gene_name,sep=''), col = "#FFC107", lwd = 4,print.auc = TRUE, auc.polygon = TRUE,auc.polygon.col = "#FFECB3")
  eval(parse(text = paste(gene_name1,'<- recordPlot()')))
  out_put_list <- append(out_put_list,c(gene_name1,NULL))

  
  auc_normal_cold_high <- c(auc_normal_cold_high,as.numeric(p$auc))
  
}

pdf("auc_normal_cold_high.pdf",width=40,height=60)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=8))',sep='')))

dev.off()


auc_normal_hot_low <- c()
out_put_list <- c()

for(gene_name in gene_names){
  
  dat <- fread(paste('D:/R/ML_PAAD/data1/alldata/exp_',sample_name,'.csv',sep=''), header = T)
  dat$V1 <- gsub("-", "_", dat$V1)
  temp_exp <- as.data.frame(dat)
  temp_exp <- temp_exp[which(temp_exp[,c(1)]==gene_name),]
  temp_exp <- temp_exp[,-1]
  expc0 <- temp_exp[,match(hot_low,colnames(temp_exp))]
  exp_hot_low <- as.numeric(as.vector(t(expc0)))
  
  expc0 <- temp_exp[,match(cold_high,colnames(temp_exp))]
  exp_cold_high <- as.numeric(as.vector(t(expc0)))
  
  normal_expression <- exp_hot_low
  disease_expression <- exp_cold_high
  
  # 合并两组表达量数据，并创建对应的类别标签 (0表示正常组，1表示疾病组)
  all_expression <- c(normal_expression, disease_expression)
  labels <- factor(c(rep(0, length(normal_expression)), rep(1, length(disease_expression))))
  
  # 计算 ROC 曲线的真阳性率和假阳性率
  roc_curve <- roc(labels, all_expression)
  
  # 绘制 ROC 曲线
  pdf(paste('hot_low/normal_hot_low_',gene_name,'.pdf',sep=''))
  p <- plot(roc_curve, main = paste("ROC Curve for Gene Expression of ",gene_name,sep=''), col = "#FFC107", lwd = 4,print.auc = TRUE, auc.polygon = TRUE,auc.polygon.col = "#FFECB3")
  print(p)
  dev.off()
  
  #输入dev.off()直到报错
  plot(roc_curve, main = paste("",gene_name,sep=''), col = "#FFC107", lwd = 4,print.auc = TRUE, auc.polygon = TRUE,auc.polygon.col = "#FFECB3")
  eval(parse(text = paste(gene_name,'<- recordPlot()')))
  out_put_list <- append(out_put_list,c(gene_name,NULL))
  
  
  auc_normal_hot_low <- c(auc_normal_hot_low,as.numeric(p$auc))
  
}

pdf("auc_normal_hot_low.pdf",width=40,height=60)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=8))',sep='')))

dev.off()



dat <- cbind(Gene = gene_names,auc_normal_cold_high,auc_normal_hot_low)

write.csv(dat,file=paste("auc.csv",sep=''),row.names = T)
