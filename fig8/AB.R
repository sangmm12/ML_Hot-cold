setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/clusterhc_riskscore/KM4")


library(data.table)
riskscore <- fread("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/queue/risk/TCGA_ICGC/risk.csv",header = T)
riskscore <- as.data.frame(riskscore)


cluster <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/CIBERSORT/clusterresult/clusterresult_hc.csv",sep=''),header =F)
cluster <- as.data.frame(cluster)

sur <-  fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_TCGA_ICGC.csv',sep=''),header=T)
sur <- as.data.frame(sur)


#根据SampleName列进行匹配
all_name <- names(which(table(c(riskscore$V1,cluster[,1],sur$SampleName ))==3))

#处理原表格
dat_cluster <- cluster[match(all_name,cluster[,1]),]
dat_riskscore<- riskscore[match(all_name,riskscore$V1),]
dat_sur<- sur[match(all_name,sur$SampleName),]

#合并为新表格
dat<- data.frame(dat_riskscore[,-4:-5],cluster = dat_cluster[,2],Time = dat_sur[,2],Status = dat_sur[,3])



#########################################################################
high_high <- subset(dat,dat$group=="High" & dat$cluster=="hot")
high_low <- subset(dat,dat$group=="High" & dat$cluster=="cold")
low_high <- subset(dat,dat$group=="Low" & dat$cluster=="hot")
low_low <- subset(dat,dat$group=="Low" & dat$cluster=="cold")

class1 <- rep("hot-high", times=length(rownames(high_high)))
class2 <- rep("cold-high", times=length(rownames(high_low)))
class3 <- rep("hot-low", times=length(rownames(low_high)))
class4 <- rep("cold-low", times=length(rownames(low_low)))

length(rownames(high_high))
length(rownames(high_low))
length(rownames(low_high))
length(rownames(low_low))

group<- c(class1,class2,class3,class4)


KM4 <- rbind(high_high,high_low,low_high,low_low)

KM4 <- cbind(KM4,group)

write.csv(KM4,file=paste("sur_clusterhc_riskscore_KM4.csv",sep=''),row.names = F)



