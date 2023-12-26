#source('/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/queue.R')


##设置目录
#setwd('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.3/up/queue')
setwd('/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/')
##载入R包
library(dplyr)
library(tibble)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(neuralnet)
library(xgboost)
library(Matrix)
library(xgboostExplainer)
library(ggplot2)
library(ggrepel)
library(data.table)

Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}

##因为有些GSE数据库问题，cox的最大深度以及SPC的fold设置为3
coxphitermax = 5
SPCnfold = 10
#final_result <- data.frame(TCGA_ICGC=c(1),TCGA_ICGC_hot=c(1),TCGA_ICGC_cold=c(1),TCGA=c(1),ICGC_CA=c(1),ICGC_AU=c(1),GSE85916=c(1),GSE28735=c(1),GSE57495=c(1),GSE62452=c(1),GSE71729=c(1),GSE78229=c(1))
final_result <- data.frame(TCGA_ICGC=c(1),TCGA_ICGC_hot=c(1),TCGA_ICGC_cold=c(1),TCGA=c(1),ICGC_CA=c(1),ICGC_AU=c(1))
list_data = list()

#sample_name <- 'GSE85916'
#sample_name <- 'TCGA_ICGC_hot'
#sample_name <- 'TCGA_ICGC'
#for(sample_name in c("TCGA_ICGC","TCGA_ICGC_hot","TCGA_ICGC_cold","TCGA","ICGC_CA","ICGC_AU","GSE85916","GSE28735","GSE57495","GSE62452","GSE71729","GSE78229"))
for(sample_name in c("TCGA_ICGC","TCGA_ICGC_hot","TCGA_ICGC_cold","TCGA","ICGC_CA","ICGC_AU"))
{
  # all_sur_data <- read.table(paste('~/machine_learning/data/',sample_name,'/LAML_patient.txt',sep=''),header=T)
  # gene_exp <- read.csv(paste('~/machine_learning/data/',sample_name,'/survival_out.csv',sep=''),header=T)
  
  #all_sur_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name,'.csv',sep=''),header=T)
  #my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name,'_HRp.csv',sep=''), header = T)
  
  all_sur_data <- fread(paste('/home/sangmm/projects/ML_PAAD/queue/sur/sur_',sample_name,'.csv',sep=''),header=T)
  my_data <- fread(paste('/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/exp_',sample_name,'_HRp.csv',sep=''), header = T)
  
  all_sur_data <- as.data.frame(all_sur_data)
  gene_exp <- as.data.frame(my_data)
  
  print(length(colnames(gene_exp)))
  
  #gene_exp = gene_exp[,which(!colnames(gene_exp)=='ADGRE5')]
  #if(max(sapply(gene_exp[,-c(1)],max))<10){gene_exp[,-c(1)] = 2^gene_exp[,-c(1)]}
  # temp_sm = gene_exp[,c(1)]
  # temp = abs(as.numeric(min(sapply(gene_exp[,-c(1)],min))))+1
  # gene_exp = log2(gene_exp[,-c(1)]+temp)
  # gene_exp = Fun(gene_exp)
  # gene_exp$SampleName = temp_sm
  mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
  colnames(mixed) <- sub('Time','OS.time',colnames(mixed))
  colnames(mixed) <- sub('Status','OS',colnames(mixed))
  temp_row <- mixed$SampleName
  mixed<- mixed %>% 
    subset(select = -c(SampleName)) %>% 
    lapply(as.numeric) %>% 
    as.data.frame()
  mixed$SampleName <- temp_row
  mixed <- select(mixed, SampleName,OS, OS.time, everything())
  mixed <- na.omit(mixed)
  eval(parse(text=paste0('list_data$',sample_name,'=mixed')))
}
mixed = list_data$TCGA_ICGC


#70% 30%划分
# train_idx <- sample(1:nrow(mixed), 0.7 * nrow(mixed)) # 70% 的数据作为训练集
# mixed_test <- mixed[-train_idx, ]
# mixed <- mixed[train_idx, ]

##所有预测变量
#RS_COXBOOST
#RS_STEPWISE
#RS_RIDGE
#RS_LASSO
#RS_SVM
#RS_GBDT
#RS_SPC
#RS_PLS
#RS_ANN
#RS_XGBOOST

seed <- 123456
predict_list = c("RS_COXBOOST","RS_STEPWISE","RS_RIDGE","RS_LASSO","RS_SVM","RS_GBDT","RS_SPC","RS_PLS","RS_ANN","RS_XGBOOST")
predict_list_name = c("CoxBoost","Stepwise Cox","Ridge","Lasso","Survival SVM","GBDT","Supervised PCA","plsRcox","ANN","XGboost")
comb <- combn(predict_list,2)

### CoxBoost ####
set.seed(seed)
pen <- optimCoxBoostPenalty(mixed$OS.time,mixed$OS,as.matrix(mixed[,-c(1,2,3)]),
                            trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(mixed$OS.time,mixed$OS,as.matrix(mixed[,-c(1,2,3)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(mixed$OS.time,mixed$OS,as.matrix(mixed[,-c(1,2,3)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
predict_comb = function(x)
{
  as.numeric(predict(fit,newdata=x[,-c(1,2,3)],newtime=x[,3], newstatus=x[,2], type="lp"))
}
predict_x = function(x)
{
  rs <- cbind(x[,2:3],RS = as.numeric(predict(fit,newdata=x[,-c(1,2,3)],newtime=x[,3], newstatus=x[,2], type="lp")))
  cc <- summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[[1]]
}
temp_vector = c()
RS_COXBOOST = lapply(list_data,predict_comb)
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],'CoxBoost')
rm(list=c('fit','predict_comb','temp_vector'))




#### Stepwise Cox ####
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,mixed[,-c(1)],iter.max=coxphitermax),direction = direction)
  if(direction == 'forward') 
  {
    predict_comb = function(x)
    {
      as.numeric(predict(fit,type = 'risk',newdata = x[,-c(1)]))
    }
    RS_STEPWISE = lapply(list_data,predict_comb)
  }
  predict_x = function(x)
  {
    rs <- cbind(x[,c(2,3)],RS = as.numeric(predict(fit,type = 'risk',newdata = x[,-c(1)])))
    cc <- summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
  }
  temp_vector = c()
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Stepwise Cox','[',direction,']')))
  rm(list=c('fit','predict_comb','temp_vector'))
}


#### Lasso,Ridge,Enet ####
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit <- cv.glmnet(as.matrix(mixed[,-c(1,2,3)]), as.matrix(Surv(mixed$OS.time,mixed$OS)),family = "cox",alpha=alpha,nfolds = 10)
  predict_comb = function(x)
  {
    as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2,3)]),s=fit$lambda.min))
  }
  predict_x = function(x)
  {
    rs <- cbind(x[,2:3],RS = as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2,3)]),s=fit$lambda.min)))
    summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
  }
  temp_vector = c()
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Enet','[a=',alpha,']',sep='')))
  if(alpha == 0)
  {
    rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Ridge'))
    RS_RIDGE <- lapply(list_data,predict_comb)
  }
  if(alpha == 1) 
  {
    rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Lasso'))
    RS_LASSO <- lapply(list_data,predict_comb)
  }
  rm(list=c('fit','predict_comb','temp_vector'))
}


#### Survival SVM ####
fit <- survivalsvm(Surv(OS.time,OS)~., data= mixed[,-c(1)], gamma.mu = 1)
predict_comb = function(x)
{
  as.numeric(predict(fit, x[,-c(1)])$predicted)
}
predict_x = function(x)
{
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit, x[,-c(1)])$predicted))
  summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
}
RS_SVM = lapply(list_data,predict_comb)
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Survival SVM'))
rm(list=c('fit','predict_comb','temp_vector'))


#### GBDT ####
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
predict_comb = function(x)
{
  as.numeric(predict(fit,x[,-c(1)],n.trees = best,type = 'link'))
}
predict_x = function(x)
{
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit,x[,-c(1)],n.trees = best,type = 'link')))
  cc <- summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
}
RS_GBDT = lapply(list_data,predict_comb)
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('GBDT'))
rm(list=c('fit','predict_comb','temp_vector'))


#### Supervised principal components ####
data <- list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
set.seed(seed)
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                     n.fold = SPCnfold,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
predict_comb = function(x)
{
  as.numeric(superpc.predict(fit,data,list(x=t(x[,-c(1,2,3)]),y=x$OS.time,censoring.status=x$OS,featurenames=colnames(x)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred)
}
predict_x = function(x)
{
  rs <- cbind(x[,2:3],RS=as.numeric(superpc.predict(fit,data,list(x=t(x[,-c(1,2,3)]),y=x$OS.time,censoring.status=x$OS,featurenames=colnames(x)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred))
  cc <- summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
}
RS_SPC = lapply(list_data,predict_comb)
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Supervised PCA'))
rm(list=c('fit','predict_comb','temp_vector'))



#### plsRcox ####
set.seed(seed)
pdf('plsRcox.pdf')
cv.plsRcox.res=cv.plsRcox(list(x=mixed[,-c(1,2,3)],time=mixed$OS.time,status=mixed$OS),nt=10,verbose = FALSE)
fit <- plsRcox(mixed[,-c(1,2,3)],time=mixed$OS.time,event=mixed$OS,nt=as.numeric(cv.plsRcox.res[5]))
predict_comb = function(x)
{
  as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2,3)]))
}
predict_x = function(x)
{
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2,3)])))
  summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
RS_PLS = lapply(list_data,predict_comb)
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('plsRcox'))
rm(list=c('fit','predict_comb','temp_vector'))
dev.off()


#ANN
set.seed(seed)
for (max in seq(5,15,1)) {
  fit <- neuralnet(OS ~ ., data=mixed[,-c(1,3)], hidden=max,learningrate=0.1,startweights = "random",linear.output = T, lifesign = "full")
  while(length(fit)<14)
  {
    fit <- neuralnet(OS ~ ., data=mixed[,-c(1,3)], hidden=max)
  }
  predict_comb = function(x)
  {
    predict(fit,newdata = x[,-c(1,3)])
  }
  predict_x = function(x)
  {
    rs <- cbind(x[,2:3],RS=predict(fit,newdata = x[,-c(1,3)]))
    as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  }
  if(max==10) {
    RS_ANN = lapply(list_data,predict_comb)
  }
  temp_vector = c()
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('ANN [hidden=',max,']')))
  rm(list=c('fit','predict_comb','temp_vector'))
}


#XGboost
set.seed(seed)
dtrain <- xgb.DMatrix(as.matrix(mixed[,-c(1,2,3)]), label = mixed$OS.time,weight = mixed$OS)
for (max in seq(1,10,1)) {
  params <- list(
    objective = "survival:cox",
    eval_metric = "cox-nloglik",
    eta = 0.01,
    max_depth = max,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
  predict_comb = function(x)
  {
    dtest <- xgb.DMatrix(as.matrix(x[,-c(1,2,3)]), label = x$OS.time,weight = x$OS)
    predict(xgb_model,dtest)
  }
  predict_x = function(x)
  {
    dtest <- xgb.DMatrix(as.matrix(x[,-c(1,2,3)]), label = x$OS.time,weight = x$OS)
    rs <- cbind(x[,2:3],RS=predict(xgb_model,dtest))
    summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
  }
  if(max==10) {RS_XGBOOST=lapply(list_data,predict_comb)}
  temp_vector = c()
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('XGboost [max_depth=',max,']')))
  rm(list=c('predict_comb'))
}


for(i in 1:ncol(comb))
{
  temp_vector = c()
  for(data_name in names(list_data))
  {
    eval(parse(text = paste('rs <- cbind(list_data$',data_name,'[,2:3],RS1=',comb[1,i],'$',data_name,',RS2=',comb[2,i],'$',data_name,')',sep='')))
    temp_vector = append(temp_vector,summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1])
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0(predict_list_name[which(predict_list==comb[1,i])],' + ',predict_list_name[which(predict_list==comb[2,i])])))
}


#### Random survival forest ####
set.seed(seed)
fit_rf <- rfsrc(Surv(OS.time,OS)~.,data = mixed[,-c(1)],
                ntree = 1000,nodesize = 5,
                splitrule = 'logrank',
                proximity = T,
                forest = T,
                seed = seed)
predict_x = function(x)
{
  rs <- cbind(x[,2:3],RS=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted))
  summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]
}
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF')))


#ANN + Random survival forest
set.seed(seed)
for (max in seq(5,15,1)) {
  fit <- neuralnet(OS ~ ., data=mixed[,-c(1,3)], hidden=max)
  while(length(fit)<14)
  {
    fit <- neuralnet(OS ~ ., data=mixed[,-c(1,3)], hidden=max)
  }
  predict_x = function(x)
  {
    rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=predict(fit,newdata = x[,-c(1,3)]))
    summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
  }
  temp_vector = c()
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste('Survival RF + ANN [hidden=',max,']',sep='')))
  rm(list=c('fit'))
}


#GBDT + Random survival forest
set.seed(seed)
fit_GB <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
              n.trees = 10000,
              interaction.depth = 3,
              n.minobsinnode = 10,
              shrinkage = 0.001,
              cv.folds = 10,n.cores = 6)
best <- which.min(fit_GB$cv.error)
set.seed(seed)
fit_GB <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
              n.trees = best,
              interaction.depth = 3,
              n.minobsinnode = 10,
              shrinkage = 0.001,
              cv.folds = 10,n.cores = 8)
predict_x = function(x)
{
  rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=as.numeric(predict(fit_GB,x[,-c(1)],n.trees = best,type = 'link')))
  summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
}
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c('Survival RF + GBDT'))
rm(list=c('fit_GB'))





#Stepwise Cox + Random survival forest
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,mixed[,-c(1)],iter.max=coxphitermax),direction = direction)
  predict_x = function(x)
  {
    tryCatch({
      rs = cbind(x[,c(2,3)],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=as.numeric(predict(fit,type = 'risk',newdata = x[,-c(1)])))
      summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
    }, error = function(e){
      rs = cbind(x[,c(2,3)],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted))
      summary(coxph(Surv(OS.time,OS)~RS1,rs))$concordance[1]
    })
  }
  temp_vector = c()
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Stepwise Cox',' [',direction,']')))
  rm(list=c('fit'))
}



#### Supervised principal components ####
data <- list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
set.seed(seed)
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                     n.fold = SPCnfold,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
RS_SPC <- data.frame(RS = as.numeric(superpc.predict(fit,data,list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred),name = "Supervised principal components")
predict_x = function(x)
{
  rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=as.numeric(superpc.predict(fit,data,list(x=t(x[,-c(1,2,3)]),y=x$OS.time,censoring.status=x$OS,featurenames=colnames(x)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred))
  summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
}
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Supervised PCA')))
rm(list=c('fit'))




#plsRcox + Random survival forest
set.seed(seed)
pdf('pls.pdf')
cv.plsRcox.res=cv.plsRcox(list(x=mixed[,-c(1,2,3)],time=mixed$OS.time,status=mixed$OS),nt=10,verbose = FALSE)
fit <- plsRcox(mixed[,-c(1,2,3)],time=mixed$OS.time,event=mixed$OS,nt=as.numeric(cv.plsRcox.res[5]))
predict_x = function(x)
{
  rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2,3)])))
  summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
}
temp_vector = c()
for(i in names(list_data))
{
  temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
}
final_result = rbind(final_result,temp_vector)
rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + plsRcox')))
rm(list=c('fit'))
dev.off()




#XGboost + Random survival forest
set.seed(seed)
dtrain <- xgb.DMatrix(as.matrix(mixed[,-c(1,2,3)]), label = mixed$OS.time,weight = mixed$OS)
for (max in seq(1,5,1)) {
  params <- list(
    objective = "survival:cox",
    eval_metric = "cox-nloglik",
    eta = 0.01,
    max_depth = max,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
  predict_x = function(x)
  {
    dtest <- xgb.DMatrix(as.matrix(x[,-c(1,2,3)]), label = x$OS.time,weight = x$OS)
    rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2=predict(xgb_model,dtest))
    summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
  }
  temp_vector = c()
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + XGboost [max_depth=',max,']')))
}


#Lasso + Random survival forest
for (alpha in seq(0,1,0.1)) {
  set.seed(seed)
  fit <- cv.glmnet(as.matrix(mixed[,-c(1,2,3)]), as.matrix(Surv(mixed$OS.time,mixed$OS)),family = "cox",alpha=alpha,nfolds = 10)
  predict_x = function(x)
  {
    rs <- cbind(x[,2:3],RS1=as.numeric(predict(fit_rf,newdata = x[,-c(1)])$predicted),RS2 = as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2,3)]),s=fit$lambda.min)))
    summary(coxph(Surv(OS.time,OS)~RS1+RS2,rs))$concordance[1]
  }
  temp_vector = c()
  for(i in names(list_data))
  {
    temp_vector = append(temp_vector,eval(parse(text=paste0('lapply(list_data,predict_x)$',i))))
  }
  final_result = rbind(final_result,temp_vector)
  rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Enet',' [a=',alpha,']')))
  if(alpha == 0) {rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Ridge')))}
  if(alpha == 1) {rownames(final_result) = c(rownames(final_result)[1:(length(rownames(final_result))-1)],c(paste0('Survival RF + Lasso')))}
  rm(list=c('fit'))
}

save(final_result,file='result_final.Rdata')

load(file='result_final.Rdata')

final_result = final_result[-which(rownames(final_result)=='1'),]
# 
# final_result$Total <-rowSums(final_result[,c(4,5,6,7,8,9)])/6
# 
# final_result = final_result[order(final_result$Total),]

# library(ggplot2)
# library(ggrepel)
# 
# data <- data.frame(
#   x = rep(1:9, each = length(final_result$TCGA_ICGC)),
#   y = rep(1:length(final_result$TCGA_ICGC), times = 9),
#   value = c(round(final_result$TCGA_ICGC, 2),
#             round(final_result$TCGA_ICGC_hot, 2),
#             round(final_result$TCGA_ICGC_cold, 2),
#             round(final_result$TCGA, 2),
#             round(final_result$ICGC_CA, 2),
#             round(final_result$ICGC_AU, 2),
#             round(final_result$GSE57495, 2),
#             round(final_result$GSE62452, 2),
#             round(final_result$GSE78229, 2)),
#   text = rownames(final_result),
#   sample = rep(colnames(final_result)[1:9]),
#   fill_column = rep(c("no_fill", "no_fill","no_fill", "no_fill", "no_fill", "no_fill", "no_fill","no_fill", "fill"), each = length(final_result$TCGA_ICGC))
# )
# 
# # 创建颜色标尺
# #color_scale <- scale_fill_gradient(low = "#FFCDD2", high = "#B71C1C")
# 
# color_scale <- scale_fill_gradient(low = "white", high = "#f47983",limits = c(0.5, 1))
# 
# pdf('Q.pdf',width=50,height=22)
# 
# ggplot(data, aes(x = x, y = y, fill = value)) +
#   geom_tile(width = 0.8, height = 0.9) +
#   geom_line(aes(x=1,y=1))+
#   theme_void()+geom_text(aes(label = value), color = "black", size = 4)+
#   theme_void()+geom_text(aes(label = text,x=x+1.61,y=y),data = subset(data, fill_column == "fill"), color = "black", size = 4)+
#   #theme_void()+geom_text(aes(label = sample,x=x,y=y+0.01),data = subset(data,fill_row), color = "black", size = 6)+
#   scale_x_continuous(expand = c(4, 0))+
#   color_scale
# 
# dev.off()
# 
# pdf('C.pdf',width=40,height=10)
# 
# ggplot(final_result, aes(x = reorder(rownames(final_result), -Total), y = Total)) +
#   geom_bar(stat = "identity",
#            show.legend = FALSE,
#            width = 0.8) + aes(fill=Total)+theme_classic()+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,size=20),axis.text.y = element_text( hjust = 1,size=20))+
#   scale_fill_gradient(low = "#a4e2c6", high = "#057748")
# 
# dev.off()

write.csv(final_result,file=paste("final_result.csv",sep=''),row.names = T)
