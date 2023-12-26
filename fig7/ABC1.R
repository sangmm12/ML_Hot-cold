  #source("/home/sangmm/projects/ML_PAAD/algorithm/algorithm_rdata.R")
  
  
  for(sample_name in c("TCGA_ICGC","TCGA_ICGC_hot","TCGA_ICGC_cold")){
    for(genelist in c("down")){
      
      library(data.table)
      ##########ann
      library(neuralnet)
      library(dplyr)
      
      #sample_name <- 'TCGA_ICGC'
      #genelist <- 'down'
      setwd(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/rstudio-export/',sample_name,'/',genelist,sep=''))
      all_sur_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name,'.csv',sep=''),header=T)
      #all_sur_data <- fread(paste('/home/sangmm/projects/ML_PAAD/queue/sur/sur_',sample_name,'.csv',sep=''),header=T)
      all_sur_data <- as.data.frame(all_sur_data)
      my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/',genelist,'/data',genelist,'/exp_',sample_name,'_HRp.csv',sep=''), header = T)
      #my_data <- fread(paste('/home/sangmm/projects/ML_PAAD/queue/',genelist,'/exp_',sample_name,'.csv',sep=''), header = T)
      gene_exp <- as.data.frame(my_data)
      mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
      rownames(mixed) <- mixed$SampleName
      mixed<- subset(mixed, select = -c(SampleName,Time))
      colnames(mixed) <- sub('Status','target',colnames(mixed))
      data <- na.omit(as.data.frame(lapply(mixed,as.numeric)))
      
      
      
      # 划分训练集和测试集
      set.seed(123) # 设置随机种子，确保结果可复现
      train_idx <- sample(1:nrow(data), 0.8 * nrow(data))
      train_data <- data[train_idx, ]
      test_data <- data[-train_idx, ]
      
      set.seed(123)
      #100%
      nn <- neuralnet(target ~ ., data=data, hidden=9)
      #%80
      nn <- neuralnet(target ~ ., data=train_data, hidden=9)
      summary(nn)
      
      train_pred <- neuralnet::compute(nn, train_data[, -ncol(train_data)])$net.result
      test_pred <- neuralnet::compute(nn, test_data[, -ncol(test_data)])$net.result
      
      train_accuracy <- sum(train_data$target == round(train_pred)) / nrow(train_data)
      test_accuracy <- sum(test_data$target == round(test_pred)) / nrow(test_data)
      
      cat("Train Accuracy:", train_accuracy, "\n")
      cat("Test Accuracy:", test_accuracy, "\n")
      weight <- data.frame(sum = rowSums(abs(nn$weights[[1]][[1]]))[-1],row.names = colnames(gene_exp)[-1])
      
      #排名
      # ANN_order <- data.frame(ANN=rep(length(weight[,1]):1),gene=rownames(weight)[order(weight$sum)])
      #权重
      Fun <- function(x){
        return ((x - min(x))/(max(x)-min(x)))
      }
      
      ANN_order <- data.frame(ANN=Fun(weight[order(weight$sum),]),gene=rownames(weight)[order(weight$sum)])
      save(ANN_order,file='ANN_order.Rdata')
      #save(ANN_order,file=paste('/home/sangmm/projects/ML_PAAD/algorithm/',sample_name,'/',genelist,'/ANN_order.Rdata',sep=''))
      
      library(ggplot2)
      pdf("ann_rank.pdf",width=15)
      p <- ggplot(weight, aes(x = reorder(rownames(weight), -sum), y = sum)) +
        geom_bar(stat = "identity",   
                 show.legend = FALSE,   
                 width = .8,fill = 'darkblue') +
        labs(x = "FeatureName", y = "weight") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      p
      print(p)
      dev.off()
      
      
      
      
      ########Boruta
      library(Boruta)
      library(survival)
      
      set.seed(123)
      library(dplyr)
      #sample_name <- 'TCGA_ICGC'
      #genelist <- 'up'
      all_sur_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name,'.csv',sep=''),header=T)
      all_sur_data <- as.data.frame(all_sur_data)
      my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/',genelist,'/data',genelist,'/exp_',sample_name,'_HRp.csv',sep=''), header = T)
      gene_exp <- as.data.frame(my_data)
      
      
      mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
      
      rownames(mixed) <- mixed$SampleName
      
      time <- mixed$Time
      status <- mixed$Status
      
      mixed <- select(mixed,-c(SampleName,Time,Status))
      
      features <- mixed
      
      surv_object <- Surv(time = time, event = status)
      
      boruta_result <- Boruta(features, surv_object)
      
      selected_features <- getSelectedAttributes(boruta_result,withTentative = TRUE)
      
      Fun <- function(x){
        return ((x - min(x))/(max(x)-min(x)))
      }
      
      Boruta_imp <- attStats(boruta_result)[order(attStats(boruta_result)$meanImp),]
      BORUTA_order <- data.frame(BORUTA=Fun(Boruta_imp$meanImp),gene=rownames(Boruta_imp))
      
      save(BORUTA_order,file='BORUTA_order.Rdata')
      
      pdf('boruta_importance.pdf',width = 20,height=10)
      par(oma=c(3,3,3,3)) 
      plot(boruta_result,las=2,xlab='')
      legend(x = 'topleft', 
             legend = c(paste('P-value:',boruta_result$pValue),sep=''),
             lty = 0,
             bty = 'n')
      dev.off()
      
      pdf('boruta_history.pdf',width = 24,height=10)
      par(oma=c(3,3,3,3)) 
      plot(plotImpHistory(boruta_result),las=2)
      legend(x = 'topleft', 
             legend = c(paste('P-value:',boruta_result$pValue),sep=''),
             lty = 0,
             bty = 'n')
      dev.off()
      
      
      
      #############randomforest
      #sample_name <- 'TCGA_ICGC'
      
      name_fold = sample_name
      
      all_sur_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name,'.csv',sep=''),header=T)
      all_sur_data <- as.data.frame(all_sur_data)
      my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/',genelist,'/data',genelist,'/exp_',sample_name,'_HRp.csv',sep=''), header = T)
      gene_exp <- as.data.frame(my_data)
      
      mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
      aa <- colnames(mixed)[-1]
      rownames(mixed) <- mixed$SampleName
      
      mixed<- subset(mixed, select = -c(SampleName))
      mixed <- as.data.frame(lapply(mixed,as.numeric))
      colnames(mixed) <- aa
      #########################################################################
      
      mixed_train <- mixed#[1:86,]
      mixed_test <- mixed#[87:nrow(mixed),]
      
      
      library(randomForestSRC)
      library(survival)
      
      set.seed(123456)
      
      
      
      v.obj <- rfsrc(Surv(Time,Status)~.,data = mixed_train,
                     mtry=3,
                     nodesize=5,
                     ntree=2000,
                     tree.err = TRUE,
                     importance = TRUE
      )
      
      out.rf <- var.select(object=v.obj,conservative = "high",)
      
      test <- predict(v.obj, mixed_test,importance = TRUE)
      full <- predict(v.obj, mixed,importance = TRUE)
      pdf('randomforest_out_train.pdf',height = 21,width=15)
      plot(v.obj)
      dev.off()
      pdf('randomforest_out_test.pdf',height = 21,width=15)
      plot(test)
      dev.off()
      pdf('randomforest_full_matrix.pdf',height=21,width=15)
      plot(full)
      dev.off()
      
      Fun <- function(x){
        return ((x - min(x))/(max(x)-min(x)))
      }
      
      RF_order <- data.frame(RF=Fun(out.rf$varselect$vimp),gene=rownames(out.rf$varselect))
      
      save(RF_order,file='RF_order.Rdata')
      colnames(mixed)
      var_out<- data.frame(vimp=out.rf$varselect$vimp,group=rownames(out.rf$varselect))
      
      pdf('randomforest_vimp.pdf',width=13,height=5)
      
      library(ggplot2)
      library(hrbrthemes)
      library(showtext)
      showtext_auto()
      p <- ggplot(var_out, aes(x = reorder(group, -vimp), y = vimp)) + 
        geom_bar(stat = "identity",   
                 show.legend = FALSE,   
                 width = .7,fill='darkred') + aes(fill=vimp)+
        xlab("Gene") + 
        ylab("Vimp")+  theme_classic()+
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      print(p)
      dev.off()
      
      
      var_out<- data.frame(depth=out.rf$varselect$depth,group=rownames(out.rf$varselect))
      
      pdf('randomforest_depth.pdf',width=24,height=10)
      
      library(ggplot2)
      library(hrbrthemes)
      library(showtext)
      showtext_auto()
      p <- ggplot(var_out, aes(x = reorder(group, -depth), y = depth)) + 
        geom_bar(stat = "identity",   
                 show.legend = FALSE,   
                 width = .7,fill='darkblue') + aes(fill=depth)+
        xlab("Gene") + 
        ylab("Depth")+  theme_classic()+
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      print(p)
      dev.off()
      
      data <- data.frame(
        PointName = rownames(out.rf$varselect),
        XValue = out.rf$varselect$depth,
        YValue = out.rf$varselect$vimp
      )
      
      library(ggrepel)
      pdf('randomforest_point.pdf',height=10,width=10)
      p <- ggplot(data, aes(x = XValue, y = YValue, label = PointName)) +
        geom_point() +
        geom_text_repel(vjust = -0.5) +
        labs(x = "Depth", y = "Vimp")
      print(p)
      dev.off()
      
      
      
      
      
      ############SVM
      library(dplyr)
      
      all_sur_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name,'.csv',sep=''),header=T)
      all_sur_data <- as.data.frame(all_sur_data)
      my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/',genelist,'/data',genelist,'/exp_',sample_name,'_HRp.csv',sep=''), header = T)
      gene_exp <- as.data.frame(my_data)
      
      mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
      rownames(mixed) <- mixed$SampleName
      mixed <- na.omit(mixed)
      gene_list <- colnames(mixed)
      
      mixed <- as.data.frame(apply(mixed,2,function(x) as.numeric(as.character(x))))
      colnames(mixed) <- gene_list
      mixed <- cbind(data.frame(status = mixed$Status),mixed)
      mixed <- select(mixed,-c(SampleName,Time,Status))
      
      library(e1071)
      
      source("D:/R/ML_PAAD/data0/AUCAtcga0/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.3_2/up/rstudio-export/SVM/msvmRFE.R")
      
      nfold = 10 #10倍交叉验证
      nrows = nrow(mixed)
      folds = rep(1:nfold, len=nrows)[sample(nrows)]
      folds = lapply(1:nfold, function(x) which(folds == x))
      
      
      results = lapply(folds, svmRFE.wrap, mixed, k=10, halve.above=100)
      top.features = WriteFeatures(results, mixed, save=F)
      featsweep = lapply(1:5, FeatSweep.wrap, results, mixed)
      
      Fun <- function(x){
        return ((x - min(x))/(max(x)-min(x)))
      }
      
      SVM_order <- data.frame(SVM=Fun(length(top.features$AvgRank)-top.features$AvgRank),gene=top.features$FeatureName)
      save(SVM_order,file='SVM_order.Rdata')
      
      no.info = min(prop.table(table(mixed[,1])))
      errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
      
      pdf("svm_rfe3.pdf", height = 8, width = 10)
      plot(errors,type=c("b"),col = 4, lty = 2)
      text(x = errors,y = sprintf("%.5f", errors),label  = sprintf("%.5f", errors), pos = 3, offset = 0.5, col = "red")
      dev.off()
      
      
      library(ggplot2)
      pdf("svm_rank.pdf",width=13)
      p <- ggplot(top.features, aes(x = reorder(FeatureName, AvgRank), y = AvgRank)) +
        geom_bar(stat = "identity",   
                 show.legend = FALSE,   
                 width = .8,fill = 'lightblue') +
        labs(x = "FeatureName", y = "rank") + theme_classic()+
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      p
      print(p)
      dev.off()
      
      
      
      
      
      #########XGboost
      library(survival)
      library(xgboost)
      library(Matrix)
      library(xgboostExplainer)
      library(timeROC)
      library(dplyr)
      library(survminer)
      
      set.seed(2)
      all_sur_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/sur/sur_',sample_name,'.csv',sep=''),header=T)
      all_sur_data <- as.data.frame(all_sur_data)
      my_data <- fread(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/',genelist,'/data',genelist,'/exp_',sample_name,'_HRp.csv',sep=''), header = T)
      gene_exp <- as.data.frame(my_data)
      
      mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
      rownames(mixed) <- mixed$SampleName
      mixed <- na.omit(mixed)
      gene_list <- colnames(mixed)
      
      mixed <- as.data.frame(apply(mixed,2,function(x) as.numeric(as.character(x))))
      colnames(mixed) <- gene_list
      mixed <- cbind(data.frame(status = mixed$Status),mixed)
      mixed <- select(mixed,-c(SampleName))
      
      mixed_train <- mixed#[ 1:86, ]
      mixed_test <- mixed#[ 87:nrow(mixed), ]
      mixed_full <- mixed
      
      train_status <- mixed_train$Status
      test_status <- mixed_test$Status
      full_status <- mixed_full$Status
      
      train_time <- mixed_train$Time
      test_time <- mixed_test$Time
      full_time <- mixed_full$Time
      
      mixed_train <- select(mixed_train,-c(Status,Time))
      mixed_test <- select(mixed_test,-c(Status,Time))
      mixed_full <- select(mixed_full,-c(Status,Time))
      #xgb_data <- xgb.DMatrix(data = as.matrix(features), label = as.matrix(status))
      
      dtrain <- xgb.DMatrix(as.matrix(mixed_train), label = train_time,weight = train_status)
      dtest <- xgb.DMatrix(as.matrix(mixed_test), label = test_time, weight = test_status)
      dfull <- xgb.DMatrix(as.matrix(mixed_full), label = full_time, weight = full_status)
      
      
      params <- list(
        objective = "survival:cox",
        eval_metric = "cox-nloglik",
        eta = 0.01,
        max_depth = 3,
        subsample = 0.8,
        colsample_bytree = 0.8
      )
      
      # 运行XGBoost模型
      set.seed(123)
      xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
      
      # pre_xgb = round(predict(xgb_model,newdata = dtest))
      # table(test_status,pre_xgb,dnn=c("true","pre"))
      # xgboost_roc <- roc(test_status,as.numeric(pre_xgb))
      # plot(xgboost_roc, print.auc=TRUE, auc.polygon=TRUE, 
      #      grid=c(0.1, 0.2),grid.col=c("green", "red"), 
      #      max.auc.polygon=TRUE,auc.polygon.col="skyblue", 
      #      print.thres=TRUE,main='ROC curve')
      
      feature_importance <- xgb.importance(model = xgb_model)
      
      predict_test <- predict(xgb_model,dtest)
      test_time <- test_time/365
      ROC_rt=timeROC(T=test_time,delta=test_status,
                     marker=predict_test,cause=1,
                     weighting='aalen',
                     times=c(1,2,3),ROC=TRUE)
      pdf(file='XGboost_roc_test.pdf',width=5,height=5)
      plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
      plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
      plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
      legend('bottomright',
             c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
               paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
               paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
             col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
      dev.off()
      
      predict_full <- predict(xgb_model,dfull)
      full_time <- full_time/365
      ROC_rt=timeROC(T=full_time,delta=full_status,
                     marker=predict_full,cause=1,
                     weighting='aalen',
                     times=c(1,2,3),ROC=TRUE)
      pdf(file='XGboost_roc_full.pdf',width=5,height=5)
      plot(ROC_rt,time=1,col='green3',title=FALSE,lwd=2)
      plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
      plot(ROC_rt,time=3,col='darkred',add=TRUE,title=FALSE,lwd=2)
      legend('bottomright',
             c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
               paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
               paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
             col=c("green3",'blue4','darkred'),lwd=2,bty = 'n')
      dev.off()
      
      Fun <- function(x){
        return ((x - min(x))/(max(x)-min(x)))
      }
      
      feature_importance <- feature_importance[order(feature_importance$Gain),]
      XGBOOST_order <- data.frame(XGBOOST=Fun(feature_importance$Gain),gene=feature_importance$Feature)
      save(XGBOOST_order,file='XGBOOST_order.Rdata')
      print(colnames(mixed_full)[which(!colnames(mixed_full)%in%XGBOOST_order$gene)])
      pdf("XGboost_rank.pdf",width=13)
      p <- ggplot(feature_importance, aes(x = reorder(feature_importance$Feature, -feature_importance$Gain), y = feature_importance$Gain)) +
        geom_bar(stat = "identity",   
                 show.legend = FALSE,   
                 width = .8,fill = 'orange2') +
        labs(x = "FeatureName", y = "Gain") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      p
      print(p)
      dev.off()
      
      
      
    }
    
  }
  
