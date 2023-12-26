#source('/home/sangmm/projects/ML_PAAD/WGCNA/3data1020/3data_WGCNA.R')

setwd("/home/sangmm/projects/ML_PAAD/WGCNA/3data1020/")

###########

library(data.table)
data <- fread("/home/sangmm/projects/ML_PAAD/WGCNA/3data1020/3data_TPM.txt")
#data <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/3data_TPM.txt")
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- data[,-1]


cluster <- read.csv("/home/sangmm/projects/ML_PAAD/WGCNA/3data1020/clusterresult_hc.csv",header = F)
# cluster <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/clusterresult_hc.csv",header = F)


c1 <- cluster[which(cluster$V2=="hot"),]
c2 <- cluster[which(cluster$V2=="cold"),]



data_c1 <- as.matrix(data)[,na.omit(match(c1[,1],colnames(data)))]
data_c2 <- as.matrix(data)[,na.omit(match(c2[,1],colnames(data)))]


data_merge <- cbind(data_c1,data_c2)
data_merge=apply(data_merge,2,as.numeric)


rownames(data_merge) <- rownames(data)

##########

#BiocManager::install("impute")



library("WGCNA")

MAIT_data  <- data_merge

enableWGCNAThreads()

dataExpr<-MAIT_data
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad >0),]

dataExpr <- as.data.frame(t(dataExprVar))

gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)


sampleTree = hclust(dist(dataExpr,method='manhattan'), method = "average");

#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)



pdf("Plots_sampleClustering.pdf", width = 10, height = 8);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#abline(h = 15, col = "red");
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#别忘了dev.off()
dev.off()



powers = c(c(1:10), seq(from = 12, to=100, by=2))
powers

sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5,RsquaredCut = 0.85)
sft$powerEstimate







pdf("Scale independence.pdf", width = 8, height =6);

par(mfrow = c(1,2))

plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab="Soft Threshold (power)",
  ylab="Scale Free Topology Model Fit,signed R^2",type="n",
  main = paste("Scale independence")  
)

text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels=powers,
  cex=0.85,
  col="red"
)

abline(h=0.85, col="red")
cex1 = 0.85

plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  xlab="Soft Threshold (power)",
  ylab="Mean Connectivity",
  type="n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels=powers,
  cex=cex1,
  col="red"
)

dev.off()






net = blockwiseModules(dataExpr,power = sft$powerEstimate, maxBlockSize = 60000,TOMType = "unsigned", 
                       minModuleSize = 80,reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE,
                       pamRespectsDendro = FALSE,verbose = 3)

mergedColors = labels2colors(net$colors)
table(mergedColors)


pdf("plotDendroAndColors.pdf", width = 10, height =8)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()



MEList = moduleEigengenes(dataExpr, colors = mergedColors) 
###计算得到的模块特征值
MEs = MEList$eigengenes
MET = orderMEs(MEs)
rownames(MET) <- rownames(dataExpr)
head(MET)[,1:8]

pdf("plotEigengeneNetworks.pdf", width =12, height =16)
plotEigengeneNetworks(MET, '', marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), 
                      cex.lab =0.8, xLabelsAngle= 90,printAdjacency = T,cex.adjacency = 0.4)
dev.off()







nSelect = 1000
set.seed(10);
dissTOM = 1-TOMsimilarityFromExpr(dataExpr, power = sft$powerEstimate);


nGenes = ncol(dataExpr)
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = 'average')
selectColors = mergedColors[select];
plotDiss = selectTOM^7;
diag(plotDiss) = NA;

pdf("TOMplot.pdf", width = 8, height =8)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()





###构建样品矩阵
nSamples <- nrow(dataExpr)
sampleMatrix <- as.data.frame(diag(x=1,nrow=nSamples))

###设置样品矩阵行列名称
rownames(sampleMatrix) <- rownames(dataExpr)
colnames(sampleMatrix) <- rownames(dataExpr)

###查看样品矩阵和模块特征值矩阵
sampleMatrix[1:9,1:9]

###通过cor 函数计算模块特征值(MET)和样品矩阵直接的相关性
moduleSampleCor <- cor(MET,sampleMatrix,use="p")
moduleSampleCor[1:9,1:9]




###通过corPvalueStudent计算pval 值
moduleSamplePvalue <- corPvalueStudent(moduleSampleCor, nSamples)

###将corr 和 pval 整合到一起
textMatrix <- paste(signif(moduleSampleCor,2), "\n(", signif(moduleSamplePvalue,1), ")", sep="")

write.csv(moduleSampleCor,file="moduleSampleCor.csv")
write.csv(moduleSamplePvalue,file="moduleSamplePvalue.csv")




pdf("labeledHeatmap.pdf", width =12, height =12)
par(mar = c(4,7, 4, 4))
labeledHeatmap(Matrix=moduleSampleCor, xLabels=names(sampleMatrix), yLabels=names(MET),
               ySymbols=names(MET), xLabelsAngle = 90, cex.lab = 0.6,
               colorLabels=FALSE, colors=blueWhiteRed(100),
               textMatrix=textMatrix, setStdMargins=FALSE,
               cex.text=0.5, zlim=c(-1,1),
               yColorWidth=0.03, xColorWidth=0.05, main= paste("Module-Sample relationship"))

dev.off()



save(net,file="net.RData")
save(MET,file="MET.RData")
save(nSamples,file="nSamples.RData")
save(dataExpr,file="dataExpr.RData")


#########################################################
moduleColors <- labels2colors(net$colors)
modules <- moduleColors

###black 模块基因表达模式
#par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 3))
#plotMat(t(scale(dataExpr[,modules=="pink"])), nrgcols=30, rlabels=F, clabels=F,rcols="pink", main="pink", cex.main=1.5)
#par(mar=c(5, 2.5, 0, 0.7))
#ME <- MEs[,"MEred"]
#names(ME) <- rownames(dataExpr)

#pdf("ME_pink.pdf", width =12, height =12)
#barplot(ME, col="red", main="", cex.main=0.6,ylab="eigengene expression", xlab="array sample")
#dev.off()






moduleColors <- labels2colors(net$colors)
TOM <- TOMsimilarityFromExpr(dataExpr, power=sft$powerEstimate)


###Select modules 选择所有模块导出
modules <- moduleColors
allModules <- unique(moduleColors)
allModules


###Select module genes
Genes <- colnames(dataExpr)
inModule <- is.finite(match(moduleColors, modules))
modGenes <- Genes[inModule]

###Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modGenes, modGenes)

cyt <- exportNetworkToCytoscape(modTOM, edgeFile="CytoscapeInput.edge.all.txt",nodeFile="CytoscapeInput.node.all.txt",weighted=TRUE, threshold=0.2, nodeNames=modGenes, altNodeNames=modGenes, nodeAttr=moduleColors[inModule])






#############################################

for(i in 1:length(allModules)){
  ###Select modules
  modules <- allModules[i]
  
  Genes <- colnames(dataExpr)
  inModule <- is.finite(match(moduleColors, modules))
  modGenes <- Genes[inModule]
  
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modGenes, modGenes)
  
  name1<-paste("CytoscapeInput.edge.",modules,".txt",sep="")
  name2<-paste("CytoscapeInput.node.",modules,".txt",sep="")
  cyt <- exportNetworkToCytoscape(modTOM, edgeFile=name1,nodeFile=name2,weighted=TRUE, threshold=0.2, nodeNames=modGenes, altNodeNames=modGenes, nodeAttr=moduleColors[inModule])
  
  #name3<-paste("MAIT/CytoscapeInput.edge.",modules,".all.txt",sep="")
  #name4<-paste("MAIT/CytoscapeInput.node.",modules,".all.txt",sep="")
  
  #cyt <- exportNetworkToCytoscape(modTOM, edgeFile=name3,nodeFile=name4,weighted=TRUE, threshold = 0, nodeNames=modGenes, altNodeNames=modGenes, nodeAttr=moduleColors[inModule])
  
  gene <- rownames(modTOM)
  name5<-paste("gene.",modules,".txt",sep="")
  write.csv(gene ,file=name5,row.names = F,quote = F)
  
}



cell <- fread("3data.csv")
cell <- as.matrix(cell)
rownames(cell) <- cell[,1]
cell <- cell[,-1]

cell <- t(cell)



cell_type <- cell[match(colnames(data_merge),rownames(cell)),]

dataTrait<- cell_type


moduleTraitCor = cor(MET, dataTrait, use = "p")
moduleTraitCor



moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")

write.csv(moduleTraitCor,file="1_moduleTraitCor.csv")
write.csv(moduleTraitPvalue,file="1_moduleTraitPvalue.csv")



pdf("1_labeledHeatmap.pdf", width =12, height =9)
par(mar = c(6,9, 4, 4))
labeledHeatmap(Matrix = moduleTraitCor,xLabels = colnames(dataTrait),yLabels = names(MET),ySymbols = names(MET),
               colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.4,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()





MEs_col <- net$MEs

cancer_type <- cluster[match(colnames(data_merge),cluster$V1),]
rownames(cancer_type) <- cancer_type[,1]
cancer_type  <- as.matrix(cancer_type )
cancer_type <- cancer_type[,2]
cancer_type  <- as.matrix(cancer_type )
colnames(cancer_type) <- c("Cancer Type")


design <- model.matrix(~0 + cancer_type)
dimnames(design) <- list(colnames(data_merge), sort(unique(cancer_type)))
design <- design[rownames(MEs_col),]

# 计算 pearson 相关性和显著性
modTraitCor <- cor(MEs_col, design, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, dim(cancer_type)[1])


textMatrix1 = paste(signif(modTraitCor, 2), "\n(",signif(modTraitP, 1), ")", sep = "")

pdf("1_labeledHeatmap_cancer_type.pdf", width =4, height =6)
par(mar = c(6,9, 4, 4))
labeledHeatmap(Matrix = modTraitCor,xLabels = c("Cold","Hot"),yLabels = names(MET),ySymbols = names(MET),
               colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix1,setStdMargins = FALSE,cex.text = 0.4,zlim = c(-1,1),main = paste("Module-trait relationships"))
dev.off()




### 计算模块与基因的相关性矩阵
geneModuleMembership = as.data.frame(cor(dataExpr, MET, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(
  as.matrix(geneModuleMembership), nSamples))

write.csv(geneModuleMembership,file=paste("geneModuleMembership.csv",sep=''),row.names = T)
# 计算性状与基因的相关性矩阵
## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。
geneTraitCor = as.data.frame(cor(dataExpr, design, use = "p"))
geneTraitP = as.data.frame(corPvalueStudent(
  as.matrix(geneTraitCor), nSamples))


#save.image(file = "rdata.RData")
#load(file = "data/projectimage.RData")

module = "brown"
pheno = "hot"
modNames = substring(colnames(MET), 3)

# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
# 获取模块内的基因
moduleGenes = moduleColors == module

# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf("1_verboseScatterplot_brown_hot.pdf", width =6, height =6)
#sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(geneModuleMembership[moduleGenes, module_column],
                   geneTraitCor[moduleGenes, pheno_column],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 20)
dev.off()

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_brown_hot_0.5_0.2.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.25)]
length(all_name)
write.csv(all_name,file=paste("all_name_brown_hot_0.5_0.25.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.3)]
length(all_name)
write.csv(all_name,file=paste("all_name_brown_hot_0.5_0.3.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.6 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_brown_hot_0.6_0.2.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.6 & geneTraitCor[moduleGenes, pheno_column] > 0.3)]
length(all_name)
write.csv(all_name,file=paste("all_name_brown_hot_0.6_0.3.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.7 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_brown_hot_0.7_0.2.csv",sep=''),row.names = F)




module = "pink"
pheno = "cold"
modNames = substring(colnames(MET), 3)

# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
# 获取模块内的基因
moduleGenes = moduleColors == module

# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf("1_verboseScatterplot_pink_cold.pdf", width =6, height =6)
#sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(geneModuleMembership[moduleGenes, module_column],
                   geneTraitCor[moduleGenes, pheno_column],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 20)
dev.off()

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_pink_cold_0.5_0.2.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.25)]
length(all_name)
write.csv(all_name,file=paste("all_name_pink_cold_0.5_0.25.csv",sep=''),row.names = F)


gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.3)]
length(all_name)
write.csv(all_name,file=paste("all_name_pink_cold_0.5_0.3.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.6 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_pink_cold_0.6_0.2.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.6 & geneTraitCor[moduleGenes, pheno_column] > 0.3)]
length(all_name)
write.csv(all_name,file=paste("all_name_pink_cold_0.6_0.3.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.7 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_pink_cold_0.7_0.2.csv",sep=''),row.names = F)





module = "turquoise"
pheno = "cold"
modNames = substring(colnames(MET), 3)

# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
# 获取模块内的基因
moduleGenes = moduleColors == module

# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf("1_verboseScatterplot_turquoise_cold.pdf", width =6, height =6)
#sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(geneModuleMembership[moduleGenes, module_column],
                   geneTraitCor[moduleGenes, pheno_column],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 20)
dev.off()

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_turquoise_cold_0.5_0.2.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.25)]
length(all_name)
write.csv(all_name,file=paste("all_name_turquoise_cold_0.5_0.25.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.3)]
length(all_name)
write.csv(all_name,file=paste("all_name_turquoise_cold_0.5_0.3.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.6 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_turquoise_cold_0.6_0.2.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.6 & geneTraitCor[moduleGenes, pheno_column] > 0.3)]
length(all_name)
write.csv(all_name,file=paste("all_name_turquoise_cold_0.6_0.3.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.7 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_turquoise_cold_0.7_0.2.csv",sep=''),row.names = F)







module = "black"
pheno = "hot"
modNames = substring(colnames(MET), 3)

# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
# 获取模块内的基因
moduleGenes = moduleColors == module

# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf("1_verboseScatterplot_black_hot.pdf", width =6, height =6)
#sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(geneModuleMembership[moduleGenes, module_column],
                   geneTraitCor[moduleGenes, pheno_column],
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 20)
dev.off()

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_black_hot_0.5_0.2.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.25)]
length(all_name)
write.csv(all_name,file=paste("all_name_black_hot_0.5_0.25.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.5 & geneTraitCor[moduleGenes, pheno_column] > 0.3)]
length(all_name)
write.csv(all_name,file=paste("all_name_black_hot_0.5_0.3.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.6 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_black_hot_0.6_0.2.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.6 & geneTraitCor[moduleGenes, pheno_column] > 0.3)]
length(all_name)
write.csv(all_name,file=paste("all_name_black_hot_0.6_0.3.csv",sep=''),row.names = F)

gene_brown <- colnames(dataExpr)[moduleGenes]
all_name <- gene_brown[which(geneModuleMembership[moduleGenes, module_column] > 0.7 & geneTraitCor[moduleGenes, pheno_column] > 0.2)]
length(all_name)
write.csv(all_name,file=paste("all_name_black_hot_0.7_0.2.csv",sep=''),row.names = F)


