
setwd('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/linkET_gene')

library(linkET)
#devtools::install_github ("Hy4m/linkET", force = TRUE)
library(ggplot2)
library(dplyr)
library(GSVA)
library(data.table)

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


genedat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/KM/single_coxTCGA_ICGC_up_p_HR.csv",sep=''),header = T)
geneup <- genedat$genename
genedat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/KM/single_coxTCGA_ICGC_down_p_HR.csv",sep=''),header = T)
genedown <- genedat$genename
geneupdown <- c(geneup,genedown)
riskscore_geneset <- genedown



##################exp
sample_name <- "TCGA_ICGC"
CIBER <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name,'.csv',sep=''), header = T)
CIBER <- as.data.frame(CIBER)
colnames(CIBER)[1] <- 'SampleName'
CIBER <- CIBER[,match(c('SampleName',riskscore_geneset),colnames(CIBER))]


######################riskscore

sample_name <- "TCGA_ICGC"
riskscore <- fread("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/queue/risk/TCGA_ICGC/risk.csv",header = T)
riskscore <- as.data.frame(riskscore)
dat_riskscore <- riskscore[,c(1,2)]
colnames(dat_riskscore)[1] <- 'SampleName'

###################ssGSVA
sample_name <- "TCGA_ICGC"
dat <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/ssgsva/',sample_name,'_genedown.csv',sep=''), header = T)
dat <- as.data.frame(dat)
dat <- dat[!duplicated(dat[,c(1)]),]
rownames(dat) <- dat[,c(1)]
dat <- dat[,-c(1)]
dat <- na.omit(dat)




#all
######################################################
#ssgsea <- dat[,which(colnames(dat)%in%hot_low)]
ssgsea <- dat
ssgsea <-as.data.frame(t(ssgsea))
ssgsea$SampleName <- rownames(ssgsea)

out <- merge(CIBER,ssgsea,by='SampleName')
out <- merge(out,dat_riskscore,by='SampleName')
ssgsea <- data.frame(risk = out$riskscore)
riskscore <- data.frame(risk = out$genedown)
rownames(out) <- out$SampleName
out <- subset(out,select=-c(riskscore,SampleName,genedown))


###################
data.corr <- qcorrplot(correlate(out), type = "lower", diag = FALSE)
data_cor <- data.corr$data  # 相关系数

write.csv(data_cor,file= paste("cor_all_data.csv",sep=''),quote=F)




mantel1 <- mantel_test(ssgsea, out,
                       spec_select = list(ssGSVA = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色


mantel2 <- mantel_test(riskscore, out,
                       spec_select = list(riskscore = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色

mantel <- rbind(mantel1,mantel2)

pdf(paste('linkET_gene_all.pdf',sep=''),width=length(colnames(CIBER))/4,height=length(colnames(CIBER))/4)
qcorrplot(correlate(out), type = "lower", diag = FALSE,grid_col = "grey88",grid_size = 0.1,) +#热图绘制
  geom_square() +#热图绘制
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes里面是线条格式，data对应的是mantel test 计算结果，curvature控制线条曲率
  scale_fill_gradientn(colours = c( "white","#F0FF8A","#FFCE74","#FFA78A","#C90500"),limits = c(0, 1)) +
  #scale_fill_gradientn(colours = c( "#2C9F00","#C7FF94","#FFA78A","#FF2300"),limits = c(0.2, 1)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",##guides()函数调整标签样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()







#hot_low
######################################################
ssgsea <- dat[,which(colnames(dat)%in%hot_low)]
ssgsea <-as.data.frame(t(ssgsea))
ssgsea$SampleName <- rownames(ssgsea)

ssg <- ssgsea
colnames(ssg)[1] <- "ssgsea"
write.csv(ssg,file= paste("ssgsea_hot_low_data.csv",sep=''),quote=F)


out <- merge(CIBER,ssgsea,by='SampleName')
out <- merge(out,dat_riskscore,by='SampleName')
ssgsea <- data.frame(risk = out$riskscore)
riskscore <- data.frame(risk = out$genedown)

risks <- riskscore
rownames(risks) <- out$SampleName
write.csv(risks,file= paste("riskscore_hot_low_data.csv",sep=''),quote=F)

rownames(out) <- out$SampleName
out <- subset(out,select=-c(riskscore,SampleName,genedown))

###################
data.corr <- qcorrplot(correlate(out), type = "lower", diag = FALSE)
data_cor <- data.corr$data  # 相关系数

write.csv(data_cor,file= paste("cor_hot_low_data.csv",sep=''),quote=F)


mantel1 <- mantel_test(ssgsea, out,
                       spec_select = list(ssGSVA = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色


mantel2 <- mantel_test(riskscore, out,
                       spec_select = list(riskscore = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色

mantel <- rbind(mantel1,mantel2)

pdf(paste('linkET_gene_hot_low.pdf',sep=''),width=length(colnames(CIBER))/4,height=length(colnames(CIBER))/4)
qcorrplot(correlate(out), type = "upper", diag = FALSE,grid_col = "grey88",grid_size = 0.1,) +#热图绘制
  geom_square() +#热图绘制
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes里面是线条格式，data对应的是mantel test 计算结果，curvature控制线条曲率
  scale_fill_gradientn(colours = c( "white","#F0FF8A","#FFCE74","#FFA78A","#C90500"),limits = c(0, 1)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",##guides()函数调整标签样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()

ssgsea <- dat[,which(colnames(dat)%in%cold_high)]
ssgsea <-as.data.frame(t(ssgsea))
ssgsea$SampleName <- rownames(ssgsea)

ssg <- ssgsea
colnames(ssg)[1] <- "ssgsea"
write.csv(ssg,file= paste("ssgsea_cold_high_data.csv",sep=''),quote=F)

out <- merge(CIBER,ssgsea,by='SampleName')
out <- merge(out,dat_riskscore,by='SampleName')
ssgsea <- data.frame(risk = out$riskscore)
riskscore <- data.frame(risk = out$genedown)

risks <- riskscore
rownames(risks) <- out$SampleName
write.csv(risks,file= paste("riskscore_cold_high_data.csv",sep=''),quote=F)






#cold_high
######################################################

rownames(out) <- out$SampleName
out <- subset(out,select=-c(riskscore,SampleName,genedown))


###################
data.corr <- qcorrplot(correlate(out), type = "lower", diag = FALSE)
data_cor <- data.corr$data  # 相关系数

write.csv(data_cor,file= paste("cor_cold_high_data.csv",sep=''),quote=F)

mantel1 <- mantel_test(ssgsea, out,
                       spec_select = list(ssGSVA = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色


mantel2 <- mantel_test(riskscore, out,
                       spec_select = list(riskscore = 1:1)) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
         pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色

mantel <- rbind(mantel1,mantel2)

pdf(paste('linkET_gene_cold_high.pdf',sep=''),width=length(colnames(CIBER))/4,height=length(colnames(CIBER))/4)
qcorrplot(correlate(out), type = "lower", diag = FALSE,grid_col = "grey88",grid_size = 0.1,) +#热图绘制
  geom_square() +#热图绘制
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes里面是线条格式，data对应的是mantel test 计算结果，curvature控制线条曲率
  scale_fill_gradientn(colours = c( "white","#F0FF8A","#FFCE74","#FFA78A","#C90500"),limits = c(0, 1)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",##guides()函数调整标签样式
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()
