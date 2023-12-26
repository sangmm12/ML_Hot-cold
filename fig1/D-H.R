setwd("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1")

library(data.table)
data <- fread("D:/R/ML_PAAD/data1/AUCAtcga/3data_count/3data_count.txt")
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- data[,-1]


cluster <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/CIBERSORT/clusterresult/clusterresult_hc.csv",header = F)


c1 <- cluster[which(cluster$V2=="cold"),]
c2 <- cluster[which(cluster$V2=="hot"),]



data_c1 <- as.matrix(data)[,na.omit(match(c1[,1],colnames(data)))]
data_c2 <- as.matrix(data)[,na.omit(match(c2[,1],colnames(data)))]


data_merge <- cbind(data_c1,data_c2)
data_merge=apply(data_merge,2,as.numeric)


rownames(data_merge) <- rownames(data)

data_merge <- data.frame(data_merge)


###########################

library(limma)
library(edgeR)


d0 <- DGEList(data_merge)
# 注意： calcNormFactors并不会标准化数据，只是计算标准化因子
d0 <- calcNormFactors(d0, method="TMM")
#logCPM<-cpm(d0,log=T,prior.count = 3)

# 过滤低表达
cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

group_list <- factor(c(rep("c1",dim(data_c1)[2]),
                       rep("c2",dim(data_c2)[2])), levels = c("c1","c2"))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(date)


v <- voom(d, design, plot=F)

fit <- lmFit(v, design)
fit <- eBayes(fit, trend=F)
output <- topTable(fit, coef=2,n=Inf)

res<-subset(output, adj.P.Val<0.05 & abs(logFC) >0.5 )


write.csv(res,file= "DEG_voom-limma.csv")







library(ggpubr)
library(ggthemes)
library(ggrepel)

deg.data <- output

#rownames(deg.data) <- gene.name

deg.data <-as.data.frame(deg.data)

deg.data$log10P <-  (-log10(deg.data$adj.P.Val))
deg.data$Symbol <- rownames(deg.data)

deg.data$Group <- "not-significant"
deg.data$Group[which((deg.data$adj.P.Val<0.01) & (deg.data$logFC > 1) )]="up-regulated"
deg.data$Group[which((deg.data$adj.P.Val<0.01) & (deg.data$logFC < -1) )]="down-regulated"
table(deg.data$Group)

deg.data$Label=""
deg.data<-deg.data[order(deg.data$adj.P.Val),]



up.genes<-head(deg.data$Symbol[which(deg.data$Group=="up-regulated")],10)
down.genes<-head(deg.data$Symbol[which(deg.data$Group=="down-regulated")],10)

deg.top10.genes<-c(as.character(up.genes),as.character(down.genes))

deg.data$Label[match(deg.top10.genes,deg.data$Symbol)]<-deg.top10.genes


pdf("volcano.pdf",width = 10)
p <- ggscatter(deg.data, x="logFC", y="log10P",
               color="Group",
               palette=c("#2f5688","#BBBBBB","#CC0000"),
               alpha=0.8, size=2.5,
               label=deg.data$Label,
               font.label=8,repel=T,
               xlab="logFC",ylab="-log10(padj)",)+
  theme_base()+
  geom_hline(yintercept=-log10(0.01),linetype="dashed")+
  geom_vline(xintercept=c(-1,1),linetype="dashed") +
  xlim(c(-3, 3)) + 
  labs(x="logFC",y="-log10(padj)")
p
dev.off()





###############################################################

library(clusterProfiler)
library(org.Hs.eg.db)



gene_name<-deg.data$Symbol[which(deg.data$Group=="up-regulated")]

write.csv(gene_name,file="DEGs_up.csv")
gene.df <- bitr(gene_name, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)


ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "all",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=0.5,
                   readable= TRUE)#"none"


write.csv(ego_plot,file="GO_up_plot.csv")


pdf("GO_up-dotplot.pdf",height = 12,width=8)
dotplot(ego_plot, split="ONTOLOGY",showCategory=5)+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()

#BP#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "BP",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=0.5,
                   readable= TRUE)#"none"

write.csv(ego_plot,file="GO_up_BP_plot.csv")


pdf("GO_up_BP-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()

#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "CC",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=0.5,
                   readable= TRUE)#"none"

write.csv(ego_plot,file="GO_up_CC_dotplot.csv")


pdf("GO_up_CC-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()

#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "MF",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=0.5,
                   readable= TRUE)#"none"

write.csv(ego_plot,file="GO_up_MF_dotplot.csv")


pdf("GO_up_MF-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()





ekk_plot<-enrichKEGG(gene=gene.df$ENTREZID,organism="hsa",
                     keyType = "kegg",pAdjustMethod="none",pvalueCutoff =0.05,qvalueCutoff=0.5)#"none"


ekk2 <- setReadable(ekk_plot, 'org.Hs.eg.db', 'ENTREZID')
write.csv(ekk2,file="KEGG_up.csv")


pdf('KEGG_up_dotplot.pdf')
dotplot(ekk2)
dev.off()







gene_name<-deg.data$Symbol[which(deg.data$Group=="down-regulated")]
write.csv(gene_name,file="DEGs_down.csv")

gene.df <- bitr(gene_name, fromType = "SYMBOL",
                toType = c("ENSEMBL","ENTREZID"),
                OrgDb = org.Hs.eg.db)


ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "all",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=0.5,
                   readable= TRUE)#"none"


write.csv(ego_plot,file="GO_down_plot.csv")


pdf("GO_down-dotplot.pdf",height = 12,width=7.7)
dotplot(ego_plot, split="ONTOLOGY",showCategory=5)+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()

#BP#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "BP",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=0.5,
                   readable= TRUE)#"none"


write.csv(ego_plot,file="GO_down_BP_plot.csv")


pdf("GO_down_BP-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()

#CC#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "CC",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"


write.csv(ego_plot,file="GO_down_CC_plot.csv")


pdf("GO_down_CC-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()

#MF
ego_plot<-enrichGO(gene=gene.df$ENTREZID,OrgDb=org.Hs.eg.db,ont = "MF",
                   pAdjustMethod ="BH", pvalueCutoff =0.05,qvalueCutoff=1,
                   readable= TRUE)#"none"


write.csv(ego_plot,file="GO_down_MF_plot.csv")


pdf("GO_down_MF-dotplot.pdf",height = 7,width=8)
dotplot(ego_plot)
dev.off()



ekk_plot<-enrichKEGG(gene=gene.df$ENTREZID,organism="hsa",
                     keyType = "kegg",pAdjustMethod="none",pvalueCutoff =0.05,qvalueCutoff=1)#"none"


ekk2 <- setReadable(ekk_plot, 'org.Hs.eg.db', 'ENTREZID')
write.csv(ekk2,file="KEGG_down.csv")


pdf('KEGG_down_dotplot.pdf')
dotplot(ekk2)
dev.off()
