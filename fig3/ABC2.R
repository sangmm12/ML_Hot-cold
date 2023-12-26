#source('/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/queue_GSE.R')

setwd(paste("/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/",sep=''))
#setwd(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/queue",sep=''))


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


GSEs <- c("GSE85916","GSE28735","GSE57495","GSE62452","GSE71729","GSE78229")


#aa=1

for (aa in 1:6) {
  
  #setwd(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/queue",sep=''))
  setwd(paste("/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/",sep=''))
  
  if(dir.exists(GSEs[aa])==TRUE){
    print(GSEs[aa])
  }else {
    dir.create(GSEs[aa])
  }
  
  #setwd(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/queue/",GSEs[aa],sep=''))
  setwd(paste("/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/",GSEs[aa],sep=''))
  #因为有些GSE数据库问题，cox的最大深度以及SPC的fold设置为3
  coxphitermax = 5
  SPCnfold = 10
  #final_result <- data.frame(TCGA_ICGC=c(1),TCGA_ICGC_hot=c(1),TCGA_ICGC_cold=c(1),TCGA=c(1),ICGC_CA=c(1),ICGC_AU=c(1),GSE85916=c(1),GSE28735=c(1),GSE57495=c(1),GSE62452=c(1),GSE71729=c(1),GSE78229=c(1))
  final_result <- data.frame(TCGA_ICGC=c(1),c(1))
  colnames(final_result)[2] <- GSEs[aa]
  list_data = list()
  
  sample_name1 <- "TCGA_ICGC"
  sample_name2 <- GSEs[aa]
  
  
  # all_sur_data1 <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name1,'.csv',sep=''),header=T)
  # all_sur_data2 <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name2,'.csv',sep=''),header=T)
  # my_data1 <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name1,'.csv',sep=''), header = T)
  # my_data2 <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/datadown/exp_',sample_name2,'.csv',sep=''), header = T)

  all_sur_data1 <- fread(paste('/home/sangmm/projects/ML_PAAD/queue/sur/sur_',sample_name1,'.csv',sep=''),header=T)
  all_sur_data2 <- fread(paste('/home/sangmm/projects/ML_PAAD/queue/sur/sur_',sample_name2,'.csv',sep=''),header=T)
  my_data1 <- fread(paste('/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/exp_',sample_name1,'.csv',sep=''), header = T)
  my_data2 <- fread(paste('/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/exp_',sample_name2,'.csv',sep=''), header = T)

  
  all_sur_data1 <- as.data.frame(all_sur_data1)
  all_sur_data2 <- as.data.frame(all_sur_data2)
  gene_exp1 <- as.data.frame(my_data1)
  gene_exp2 <- as.data.frame(my_data2)
  
  
  all_name <- names(which(table(c(colnames(gene_exp1),colnames(gene_exp2) ))==2))
  gene_exp1 <- gene_exp1[,match(all_name,colnames(gene_exp1))]
  gene_exp2 <- gene_exp2[,match(all_name,colnames(gene_exp2))]
  
  all_sur_data <- all_sur_data1
  gene_exp <- gene_exp1
  
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
  eval(parse(text=paste0('list_data$',sample_name1,'=mixed')))
  
  
  
  
  
  
  all_sur_data <- all_sur_data2
  gene_exp <- gene_exp2
  
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
  eval(parse(text=paste0('list_data$',sample_name2,'=mixed')))
  
  
  
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
  
  save(final_result,file=paste("result_final.Rdata",sep=''))
  
  #load(file='result_final.Rdata')
  
  final_result = final_result[-which(rownames(final_result)=='1'),]
  write.csv(final_result,file=paste("final_result.csv",sep=''),row.names = T)
  write.csv(final_result,file=paste("/home/sangmm/projects/ML_PAAD/0.6_0.2/queue/down/final_result_",GSEs[aa],".csv",sep=''),row.names = T)
  
}
