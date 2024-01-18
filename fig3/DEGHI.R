
library(dplyr)
library(survival)
library(Matrix)
library(glmnet)
library(plsRcox)
library(data.table)
library(dplyr)
library(randomForestSRC)
library(survival)
library(xgboost)
library(Matrix)
library(xgboostExplainer)

setwd('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/queue')

Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

final_list <- data.frame()
seed <- 123456

#sample_name <- 'TCGA_ICGC'
for(sample_name in c('TCGA_ICGC'))
{
  name_fold = paste0('risk/',sample_name)
  if (!dir.exists(name_fold)){
    dir.create(name_fold)
  } else {
    print("Dir already exists!")
  }
  
  all_sur_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name,'.csv',sep=''),header=T)
  my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name,'_HRp.csv',sep=''), header = T)
  all_sur_data <- as.data.frame(all_sur_data)
  gene_exp <- as.data.frame(my_data)
  #gene_exp = gene_exp[,which(!colnames(gene_exp)=='ADGRE5')]
  mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
  colnames(mixed) <- sub('Time','Time',colnames(mixed))
  colnames(mixed) <- sub('Status','Status',colnames(mixed))
  temp_row <- mixed$SampleName
  mixed<- mixed %>% 
    subset(select = -c(SampleName)) %>% 
    lapply(as.numeric) %>% 
    as.data.frame()
  mixed$SampleName <- temp_row
  mixed <- mixed %>% 
    select( SampleName,Status, Time, everything())
  mixed <- na.omit(mixed)
  
  
  set.seed(123456)
  p.obj <- rfsrc(Surv(Time,Status)~.,data = mixed[,-c(1)],
                 ntree = 1000,nodesize = 5,
                 splitrule = 'logrank',
                 proximity = T,
                 forest = T,
                 seed = seed)
  cv.plsRcox.res=cv.plsRcox(list(x=mixed[,-c(1,2,3)],time=mixed$Time,status=mixed$Status),nt=10,verbStatuse = FALSE)
  fit_plsr <- plsRcox(mixed[,-c(1,2,3)],time=mixed$Time,event=mixed$Status,nt=as.numeric(cv.plsRcox.res[5]))
  save(p.obj,file='p.obj.Rda')
  save(cv.plsRcox.res,file='cv.plsRcox.res.Rda')
  save(fit_plsr,file='fit_plsr.Rda')
  
  x_p = Fun(as.numeric(predict(p.obj,mixed[,-c(1,2,3)])$predicted))
  
  y_p = Fun(abs(min(as.numeric(predict(fit_plsr,type="lp",newdata=mixed[,-c(1,2,3)])))) + as.numeric(predict(fit_plsr,type="lp",newdata=mixed[,-c(1,2,3)])))
  
  threshold <- median(Risk_scroe_temp)
  
  risk <- data.frame(riskscore = Fun(Risk_scroe_temp)*100,group = risk_group <- ifelse((Risk_scroe_temp) > threshold, "High", "Low"),Status=mixed$Status,Time=mixed$Time/365)
  rownames(risk)=mixed$SampleName
  write.csv(risk,file=paste(name_fold,'/risk.csv',sep=''))
  final_list <- rbind(final_list,data.frame(CODE=rep(sample_name,length(risk$group)),SampleName=mixed$SampleName,Risk=risk$group,Status=risk$Status))
  
  rt=risk[order(risk$riskscore),]
  riskClass=rt$group
  lowLength=length(which(riskClass=="Low"))
  highLength=length(which(riskClass=="High"))
  lowMax=max(rt[which(rt$group=="Low"),]$riskscore)
  line=rt$riskscore
  line[line>1000]=1000
  pdf(file=paste(name_fold,'/risk line.pdf',sep=''),width = 6,height = 6)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("#009E73",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","#009E73"),cex=1.2)
  dev.off()
  color=as.vector(rt$Status)
  color[color==1]="red"
  color[color==0]="#009E73"
  pdf(file=paste(name_fold,'/risk point.pdf',sep=''),width = 6,height = 6)
  plot(rt$Time, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","#009E73"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  library(timeROC)
  predict_full <- risk$riskscore
  full_time <- mixed$Time/365
  ROC_rt=timeROC(T=full_time,delta=mixed$Status,
                 marker=predict_full,cause=1,
                 weighting='aalen',
                 times=c(1,2,3),ROC=TRUE)
  pdf(file=paste(name_fold,'/roc.pdf',sep=''),width = 6,height = 6)
  plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
  plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
  dev.off()
  
  

}
#write.csv(final_list,file='riskRFTcell.csv',row.names=F)





#sample_name <- 'TCGA'
for(sample_name in c("TCGA","ICGC_CA","ICGC_AU"))
{
  name_fold = paste0('risk/',sample_name)
  if (!dir.exists(name_fold)){
    dir.create(name_fold)
  } else {
    print("Dir already exists!")
  }
  
  all_sur_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name,'.csv',sep=''),header=T)
  my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name,'_HRp.csv',sep=''), header = T)
  all_sur_data <- as.data.frame(all_sur_data)
  gene_exp <- as.data.frame(my_data)
  #gene_exp = gene_exp[,which(!colnames(gene_exp)=='ADGRE5')]
  mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
  colnames(mixed) <- sub('Time','Time',colnames(mixed))
  colnames(mixed) <- sub('Status','Status',colnames(mixed))
  temp_row <- mixed$SampleName
  mixed<- mixed %>% 
    subset(select = -c(SampleName)) %>% 
    lapply(as.numeric) %>% 
    as.data.frame()
  mixed$SampleName <- temp_row
  mixed <- mixed %>% 
    select( SampleName,Status, Time, everything())
  mixed <- na.omit(mixed)
  
  
  set.seed(123456)
  # p.obj <- rfsrc(Surv(Time,Status)~.,data = mixed[,-c(1)],
  #                ntree = 1000,nodesize = 5,
  #                splitrule = 'logrank',
  #                proximity = T,
  #                forest = T,
  #                seed = seed)
  # cv.plsRcox.res=cv.plsRcox(list(x=mixed[,-c(1,2,3)],time=mixed$Time,status=mixed$Status),nt=10,verbStatuse = FALSE)
  # fit_plsr <- plsRcox(mixed[,-c(1,2,3)],time=mixed$Time,event=mixed$Status,nt=as.numeric(cv.plsRcox.res[5]))
  # save(p.obj,file='p.obj.Rda')
  # save(cv.plsRcox.res,file='cv.plsRcox.res.Rda')
  # save(fit_plsr,file='fit_plsr.Rda')

  x_p = Fun(as.numeric(predict(p.obj,mixed[,-c(1,2,3)])$predicted))
  y_p = Fun(abs(min(as.numeric(predict(fit_plsr,type="lp",newdata=mixed[,-c(1,2,3)])))) + as.numeric(predict(fit_plsr,type="lp",newdata=mixed[,-c(1,2,3)])))
  
  rs <- data.frame(Time=mixed$Time,Status=mixed$Status,RS1=x_p,RS2=y_p)
  
  Risk_scroe_temp <- summary(coxph(Surv(Time,Status)~RS1+RS2,rs))$coef[3]*x_p+summary(coxph(Surv(Time,Status)~RS1+RS2,rs))$coef[4]*y_p
  
  threshold <- median(Risk_scroe_temp)
  
  risk <- data.frame(riskscore = Fun(Risk_scroe_temp)*100,group = risk_group <- ifelse((Risk_scroe_temp) > threshold, "High", "Low"),Status=mixed$Status,Time=mixed$Time/365)
  rownames(risk)=mixed$SampleName
  write.csv(risk,file=paste(name_fold,'/risk.csv',sep=''))
  final_list <- rbind(final_list,data.frame(CODE=rep(sample_name,length(risk$group)),SampleName=mixed$SampleName,Risk=risk$group,Status=risk$Status))
  
  rt=risk[order(risk$riskscore),]
  riskClass=rt$group
  lowLength=length(which(riskClass=="Low"))
  highLength=length(which(riskClass=="High"))
  lowMax=max(rt[which(rt$group=="Low"),]$riskscore)
  line=rt$riskscore
  line[line>1000]=1000
  pdf(file=paste(name_fold,'/risk line.pdf',sep=''),width = 6,height = 6)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("#009E73",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","#009E73"),cex=1.2)
  dev.off()
  color=as.vector(rt$Status)
  color[color==1]="red"
  color[color==0]="#009E73"
  pdf(file=paste(name_fold,'/risk point.pdf',sep=''),width = 6,height = 6)
  plot(rt$Time, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","#009E73"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  library(timeROC)
  predict_full <- risk$riskscore
  full_time <- mixed$Time/365
  ROC_rt=timeROC(T=full_time,delta=mixed$Status,
                 marker=predict_full,cause=1,
                 weighting='aalen',
                 times=c(1,2,3),ROC=TRUE)
  pdf(file=paste(name_fold,'/roc.pdf',sep=''),width = 6,height = 6)
  plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
  plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
  dev.off()
  
  
  
  library(survminer)
  rt=risk[order(risk$riskscore),]
  diff=survdiff(Surv(Time, Status) ~ group,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(Time, Status) ~ group, data = rt)
  surPlot=ggsurvplot(fit,
                     data=rt,
                     #font.title = paste(connam[i],sep=''),
                     #ggtitle = paste(connam[i],sep=''),
                     #conf.int=TRUE,
                     legend.labs=c( "L","H"),
                     legend = "top",
                     legend.title="Risk",
                     pval=paste0("p=",pValue),
                     pval.size=5,
                     xlab="Time(years)",
                     break.time.by = ceiling((max(rt$Time))/4),
                     risk.table.title="",
                     palette=c("#009E73","red"),
                     risk.table=T,
                     risk.table.height=.25,)
  pdf(file=paste(name_fold,'/survival.pdf',sep=''),onefile = FALSE,width = 6,height =8)
  print(surPlot)
  dev.off()
}
#write.csv(final_list,file='riskRFTcell.csv',row.names=F)






sample_name1 <- 'TCGA_ICGC'
sample_name2 <- 'GSE62452'
for(sample_name2 in c("GSE85916","GSE28735","GSE57495","GSE62452","GSE71729","GSE78229"))
{
  all_sur_data1 <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name1,'.csv',sep=''),header=T)
  all_sur_data2 <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name2,'.csv',sep=''),header=T)
  my_data1 <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name1,'.csv',sep=''), header = T)
  my_data2 <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name2,'.csv',sep=''), header = T)
  
  
  all_sur_data1 <- as.data.frame(all_sur_data1)
  all_sur_data2 <- as.data.frame(all_sur_data2)
  gene_exp1 <- as.data.frame(my_data1)
  gene_exp2 <- as.data.frame(my_data2)
  
  all_name <- names(which(table(c(colnames(gene_exp1),colnames(gene_exp2) ))==2))
  gene_exp1 <- gene_exp1[,match(all_name,colnames(gene_exp1))]
  gene_exp2 <- gene_exp2[,match(all_name,colnames(gene_exp2))]
  
  all_sur_data <- all_sur_data1
  gene_exp <- gene_exp1

    #gene_exp = gene_exp[,which(!colnames(gene_exp)=='ADGRE5')]
  mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
  colnames(mixed) <- sub('Time','Time',colnames(mixed))
  colnames(mixed) <- sub('Status','Status',colnames(mixed))
  temp_row <- mixed$SampleName
  mixed<- mixed %>% 
    subset(select = -c(SampleName)) %>% 
    lapply(as.numeric) %>% 
    as.data.frame()
  mixed$SampleName <- temp_row
  mixed <- mixed %>% 
    select( SampleName,Status, Time, everything())
  mixed <- na.omit(mixed)
  
  
  set.seed(123456)
  p.obj <- rfsrc(Surv(Time,Status)~.,data = mixed[,-c(1)],
                  ntree = 1000,nodesize = 5,
                  splitrule = 'logrank',
                  proximity = T,
                  forest = T,
                  seed = seed)

  cv.plsRcox.res=cv.plsRcox(list(x=mixed[,-c(1,2,3)],time=mixed$Time,status=mixed$Status),nt=10,verbStatuse = FALSE)
  fit_plsr <- plsRcox(mixed[,-c(1,2,3)],time=mixed$Time,event=mixed$Status,nt=as.numeric(cv.plsRcox.res[5]))
  save(p.obj,file='p.obj.Rda')
  save(cv.plsRcox.res,file='cv.plsRcox.res.Rda')
  save(fit_plsr,file='fit_plsr.Rda')
  
  all_sur_data <- all_sur_data2
  gene_exp <- gene_exp2
  
  sample_name <- sample_name2
  name_fold = paste0('risk/',sample_name)
  if (!dir.exists(name_fold)){
    dir.create(name_fold)
  } else {
    print("Dir already exists!")
  }
  
  #gene_exp = gene_exp[,which(!colnames(gene_exp)=='ADGRE5')]
  mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
  colnames(mixed) <- sub('Time','Time',colnames(mixed))
  colnames(mixed) <- sub('Status','Status',colnames(mixed))
  temp_row <- mixed$SampleName
  mixed<- mixed %>% 
    subset(select = -c(SampleName)) %>% 
    lapply(as.numeric) %>% 
    as.data.frame()
  mixed$SampleName <- temp_row
  mixed <- mixed %>% 
    select( SampleName,Status, Time, everything())
  mixed <- na.omit(mixed)
  
  
  set.seed(123456)
  # p.obj <- rfsrc(Surv(Time,Status)~.,data = mixed[,-c(1)],
  #                mtry=3,
  #                nodesize=1,
  #                ntree=3
  # )
  # cv.plsRcox.res=cv.plsRcox(list(x=mixed[,-c(1,2,3)],time=mixed$Time,status=mixed$Status),nt=10,verbStatuse = FALSE)
  # fit_plsr <- plsRcox(mixed[,-c(1,2,3)],time=mixed$Time,event=mixed$Status,nt=as.numeric(cv.plsRcox.res[5]))
  # save(p.obj,file='p.obj.Rda')
  # save(cv.plsRcox.res,file='cv.plsRcox.res.Rda')
  # save(fit_plsr,file='fit_plsr.Rda')
  
  x_p = Fun(as.numeric(predict(p.obj,mixed[,-c(1,2,3)])$predicted))
  y_p = Fun(abs(min(as.numeric(predict(fit_plsr,type="lp",newdata=mixed[,-c(1,2,3)])))) + as.numeric(predict(fit_plsr,type="lp",newdata=mixed[,-c(1,2,3)])))
  
  rs <- data.frame(Time=mixed$Time,Status=mixed$Status,RS1=x_p,RS2=y_p)
  
  Risk_scroe_temp <- summary(coxph(Surv(Time,Status)~RS1+RS2,rs))$coef[3]*x_p+summary(coxph(Surv(Time,Status)~RS1+RS2,rs))$coef[4]*y_p
  
  threshold <- median(Risk_scroe_temp)
  
  risk <- data.frame(riskscore = Fun(Risk_scroe_temp)*100,group = risk_group <- ifelse((Risk_scroe_temp) > threshold, "High", "Low"),Status=mixed$Status,Time=mixed$Time/365)
  rownames(risk)=mixed$SampleName
  write.csv(risk,file=paste(name_fold,'/risk.csv',sep=''))
  final_list <- rbind(final_list,data.frame(CODE=rep(sample_name,length(risk$group)),SampleName=mixed$SampleName,Risk=risk$group,Status=risk$Status))
  
  rt=risk[order(risk$riskscore),]
  riskClass=rt$group
  lowLength=length(which(riskClass=="Low"))
  highLength=length(which(riskClass=="High"))
  lowMax=max(rt[which(rt$group=="Low"),]$riskscore)
  line=rt$riskscore
  line[line>1000]=1000
  pdf(file=paste(name_fold,'/risk line.pdf',sep=''),width = 6,height = 6)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("green4",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","green4"),cex=1.2)
  dev.off()
  color=as.vector(rt$Status)
  color[color==1]="red"
  color[color==0]="green4"
  pdf(file=paste(name_fold,'/risk point.pdf',sep=''),width = 6,height = 6)
  plot(rt$Time, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green4"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  library(timeROC)
  predict_full <- risk$riskscore
  full_time <- mixed$Time/365
  ROC_rt=timeROC(T=full_time,delta=mixed$Status,
                 marker=predict_full,cause=1,
                 weighting='aalen',
                 times=c(1,2,3),ROC=TRUE)
  pdf(file=paste(name_fold,'/roc.pdf',sep=''),width = 6,height = 6)
  plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
  plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
  dev.off()
  
  
  
  library(survminer)
  rt=risk[order(risk$riskscore),]
  diff=survdiff(Surv(Time, Status) ~ group,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  pValue=signif(pValue,4)
  pValue=format(pValue, scientific = TRUE)
  fit <- survfit(Surv(Time, Status) ~ group, data = rt)
  surPlot=ggsurvplot(fit,
                     data=rt,
                     #font.title = paste(connam[i],sep=''),
                     #ggtitle = paste(connam[i],sep=''),
                     #conf.int=TRUE,
                     legend.labs=c( "L","H"),
                     legend = "top",
                     legend.title="Risk",
                     pval=paste0("p=",pValue),
                     pval.size=5,
                     xlab="Time(years)",
                     break.time.by = ceiling((max(rt$Time))/4),
                     risk.table.title="",
                     palette=c("#009E73","red"),
                     risk.table=T,
                     risk.table.height=.25,)
  pdf(file=paste(name_fold,'/survival.pdf',sep=''),onefile = FALSE,width = 6,height =8)
  print(surPlot)
  dev.off()
}
#write.csv(final_list,file='riskRFTcell.csv',row.names=F)
