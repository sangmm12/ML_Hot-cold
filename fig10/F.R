
setwd(paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/zong/genescore_drug_importance',sep=''))
library(data.table)
library(GSVA)
library(randomForestSRC)
library(survival)
library(neuralnet)
library(Boruta)
library(xgboost)
library(Matrix)
library(xgboostExplainer)
library(ggplot2)
library(dplyr)
library(e1071)


sample_name <- 'TCGA_ICGC'
genelist <- 'down'

my_data <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/out_put_",sample_name,"_DEGs_down_zong.csv",sep=''))
drug_exp <- as.data.frame(my_data)
datar <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/drug/drugCTRP2_2/DEGs_downRMA/zong/cor_gene/data.r2.csv",sep=''))
drug_exp <- drug_exp[,match(colnames(datar),colnames(drug_exp))]
colnames(drug_exp)[1]<- "SampleName"


dat <- fread(paste('D:/R/ML_PAAD/data1/alldata/exp_TCGA_ICGC.txt',sep=''))
dat <- as.data.frame(dat)
dat <- dat[!duplicated(dat[,c(1)]),]
rownames(dat) <- dat[,c(1)]
dat <- dat[,-c(1)]
dat <- na.omit(dat)

genedat <- fread(paste("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/KM/single_coxTCGA_ICGC_down_p_HR.csv",sep=''),header = T)
genedown <- genedat$genename
expr1 <- as.matrix(dat[which(rownames(dat)%in%genedown),])
ssample <- list('risk'=genedown)

ssgsea <- gsva(as.matrix(dat),ssample, method='ssgsea', kcdf='Gaussian',abs.ranking=TRUE)#基因集需要是list为对象。#默认情况下，kcdf="Gaussian"，适用于输入表达式值连续的情况，如对数尺度的微阵列荧光单元、RNA-seq log-CPMs、log-RPKMs或log-TPMs。#当输入表达式值是整数计数时，比如那些从RNA-seq实验中得到的值，那么这个参数应该设置为kcdf="Poisson"
write.csv(ssgsea,file = 'GSVA.csv')
temp <- data.frame(risk = ssgsea[c(1),],SampleName = colnames(ssgsea))


Fun <- function(x){
  return ((x - min(x))/(max(x)-min(x)))
}



train_geneset_data <- merge(drug_exp,temp,by='SampleName')[,-c(1)]




p.obj <- rfsrc(risk ~.,data = train_geneset_data,
               mtry=3,
               nodesize=5,
               ntree=2000,
               tree.err = TRUE,
               importance = TRUE
)
out.rf <- var.select(object=p.obj,conservative = "high",)
RF_order <- data.frame(RF=Fun(out.rf$varselect$vimp),gene=rownames(out.rf$varselect))
full <- predict(p.obj,train_geneset_data ,importance = TRUE)
pdf('train_var_important.pdf',height=length(colnames(train_geneset_data))/3)
plot(full)
dev.off()


train_geneset_data0 <- train_geneset_data

train_geneset_data <- train_geneset_data0
colnames(train_geneset_data) <- gsub("[^A-Za-z0-9_]", "_", colnames(train_geneset_data))

as.numeric(train_geneset_data)
train_geneset_data[] <- lapply(train_geneset_data, as.numeric)
nn <- neuralnet(risk ~ ., data=train_geneset_data, hidden=c(100,100),stepmax=1e6,learningrate=0.1,startweights = "random",linear.output = T, lifesign = "full")
pdf(paste0("ANN_struct.pdf"),width=40,height=40)
p <- plot(nn,col.text = "black",col.entry='black',col.intercept='#50026E',col.out.synapse='#008209',col.hidden.synapse='#0F4FA8',col.entry.synapse='#A40004',col.out='#008209',fontsize=30, col.hidden = "#0F4FA8",radius=0.1, show.weights = FALSE,information=F)
print(p)
dev.off()
weight <- data.frame(sum = rowSums(abs(nn$weights[[1]][[1]]))[-1],row.names = colnames(train_geneset_data)[-ncol(train_geneset_data)])
rownames(weight) <- colnames(train_geneset_data0)[-length(train_geneset_data0)]
ANN_order <- data.frame(ANN=Fun(weight[order(weight$sum),]),gene=rownames(weight)[order(weight$sum)])





boruta_result <- Boruta(train_geneset_data[,c(-ncol(train_geneset_data))], train_geneset_data[,ncol(train_geneset_data)])
Boruta_imp <- attStats(boruta_result)[order(attStats(boruta_result)$meanImp),]
BORUTA_order <- data.frame(BORUTA=Fun(Boruta_imp$meanImp),gene=rownames(Boruta_imp))
pdf('boruta_importance.pdf',width = 20,height=10)
par(oma=c(3,3,3,3)) 
plot(boruta_result,las=2,xlab='')
legend(x = 'topleft', 
       legend = c(paste('P-value:',boruta_result$pValue),sep=''),
       lty = 0,
       bty = 'n')
dev.off()

pdf('boruta_history.pdf',width = 20,height=10)
par(oma=c(3,3,3,3)) 
plot(plotImpHistory(boruta_result),las=2)
legend(x = 'topleft', 
       legend = c(paste('P-value:',boruta_result$pValue),sep=''),
       lty = 0,
       bty = 'n')
dev.off()





dtrain <- xgb.DMatrix(as.matrix(train_geneset_data[,c(-ncol(train_geneset_data))]), label = train_geneset_data[,ncol(train_geneset_data)])

params <- list(
  eta = 0.01,
  max_depth = 3,
  subsample = 0.8,
  colsample_bytree = 0.8
)
xgb_model <- xgb.train(params = params,nrounds = 100,dtrain)
feature_importance <- xgb.importance(model = xgb_model)
XGBOOST_order <- data.frame(XGBOOST=Fun(feature_importance$Gain),gene=feature_importance$Feature)

colnames(train_geneset_data)[-length(colnames(train_geneset_data))]

out_gene <- colnames(train_geneset_data)[-length(colnames(train_geneset_data))][which(!(colnames(train_geneset_data)[-length(colnames(train_geneset_data))]%in%XGBOOST_order$gene))]
XGBOOST_order <- rbind(XGBOOST_order,data.frame(gene=out_gene,XGBOOST=rep(0,length(out_gene))))





data = cbind(risk=train_geneset_data$risk,train_geneset_data[,-ncol(train_geneset_data)])
data$risk = ifelse(data$risk >  median(data$risk), 1, 0)
source("D:/R/ML_PAAD/data0/AUCAtcga0/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.3_2/up/rstudio-export/SVM/msvmRFE.R")
nfold = 10
nrows = nrow(data)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))

results = lapply(folds, svmRFE.wrap, data, k=10, halve.above=100)
top.features = WriteFeatures(results, data, save=F)
SVM_order <- data.frame(SVM=Fun(length(top.features$AvgRank)-top.features$AvgRank),gene=top.features$FeatureName)








out <- ANN_order %>% left_join(BORUTA_order,by='gene') %>% left_join(RF_order,by='gene') %>% left_join(XGBOOST_order,by='gene') %>% left_join(SVM_order,by='gene') 
rownames(out) <- out$gene
out <- subset(out,select=-c(gene))
out <- cbind(Mean = rowSums(out)/length(colnames(out)),out)
out <- out[order(-out$Mean),]*3
value <- c()

for(i in 1:length(colnames(out)))
{
  eval(parse(text=paste('value<-append(value,',out[,i],')',sep='')))
}

# color <- c()
# COLORS <- c("#70f3ff","#44cef6","#3eede7","#1685a9","#177cb0","#065279","#003472","#4b5cc4","#a1afc9","#2e4e7e","#3b2e7e","#4a4266","#426666","#425066","#574266","#8d4bbb","#815463","#815476","#4c221b","#003371","#56004f","#801dae","#4c8dae","#b0a4e3","#cca4e3","#edd1d8","#e4c6d0","#ff461f","#ff2d51","#f36838","#ed5736","#ff4777","#f00056","#ffb3a7","#f47983","#db5a6b","#c93756","#f9906f","#f05654","#ff2121","#f20c00","#8c4356","#c83c23","#9d2933","#ff4c00","#ff4e20","#f35336","#dc3023","#ff3300","#cb3a56","#a98175","#b36d61","#ef7a82","#ff0097","#c32136","#be002f","#c91f37","#bf242a","#c3272b","#9d2933","#60281e","#622a1d","#bce672","#c9dd22","#bddd22","#afdd22","#a3d900","#9ed900","#9ed048","#96ce54","#00bc12","#0eb83a","#0eb83a","#0aa344","#16a951","#21a675","#057748","#0c8918","#00e500","#40de5a","#00e079","#00e09e","#3de1ad","#2add9c","#2edfa3","#7fecad","#a4e2c6","#7bcfa6","#1bd1a5","#48c0a3","#549688","#789262","#758a99","#50616d","#424c50","#41555d","#eaff56","#fff143","#faff72","#ffa631","#ffa400","#fa8c35","#ff8c31","#ff8936","#ff7500","#ffb61e","#ffc773","#ffc64b","#f2be45","#f0c239","#e9bb1d","#d9b611","#eacd76","#eedeb0","#d3b17d","#e29c45","#a78e44","#c89b40","#ae7000","#ca6924","#b25d25","#b35c44","#9b4400","#9c5333","#a88462","#896c39","#827100","#6e511e","#7c4b00","#955539","#845a33","#ffffff","#e9e7ef")
# 
# 
# for(j in 1:length(colnames(out)))
# {
#   eval(parse(text=paste('color <- append(color,c(',colnames(out)[j],'="',sample(COLORS,size=1),'"))',sep='')))
# }


widelength <- length(rownames(XGBOOST_order))
color <- c("#4b5cc4","#dc3023","#057748","#fff143","#758a99","#177cb0")

my_data <- data.frame(
  Category = rep(c(colnames(out)[1], colnames(out)[2], colnames(out)[3], colnames(out)[4],colnames(out)[5],colnames(out)[6]), each = length(out[,1])),
  Sample = rep(1:(length(colnames(train_geneset_data))-1), times = length(colnames(out))),
  Value = value,
  Color = color,
  name = rownames(out)
)

my_data$Category = factor(my_data$Category,levels = colnames(out))

color[1] = "#4c221b"

pdf('algorithm.pdf',width=widelength/3+2,height=10)
ggplot(my_data, aes(x = as.factor(Sample), y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_grid(Category ~ ., scales = "free_y", space = "free") +
  scale_fill_manual(values = color) +
  labs(title = "", x = "Drug", y = "Important") +
  scale_x_discrete(labels = my_data$name) +  
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,size=14,color='black', hjust = 1),axis.text.y = element_text(size=14,color='black'))
dev.off()

#write.csv(out,file=paste('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/rstudio-export/',sample_name,'_',genelist,'_out.csv',sep=''),row.names = T)


