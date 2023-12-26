
setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/GSVA_hallmark")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
#BiocManager::install("GSVA",type="binary")
#BiocManager::install("msigdbr")
library(GSVA)
library(GSEABase)
library(limma)
library(Seurat)
library(msigdbr)

#all_gene_sets = msigdbr(species = "Mus musculus")
human <- msigdbr(species = "Homo sapiens")
#hallmark
human_GO_bp = msigdbr(species = "human",
                      category = "H") %>% 
  dplyr::select(gs_name,gene_symbol)


human_GO_bp_Set = human_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)


s.sets <- human_GO_bp_Set


###########

library(data.table)
data <- fread("D:/R/ML_PAAD/data1/AUCAtcga/3data_count/3data_count.txt")
data <- as.data.frame(data)
rownames(data) <- data[,1]
data <- data[,-1]

rf_risk <- read.csv('D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec/0.6_0.2/down/clusterhc_riskscore/KM4/sur_clusterhc_riskscore_KM4.csv')
rf_risk <- rf_risk[,c(1,7)]
cluster <- subset(rf_risk, group %in% c("hot-low", "cold-high"))
cluster <- as.data.frame(cluster)
rownames(cluster) <- cluster[,1]
colnames(cluster) <- c("V1","V2")



c1 <- cluster[which(cluster$V2=="hot-low"),]
c2 <- cluster[which(cluster$V2=="cold-high"),]



data_c1 <- as.matrix(data)[,na.omit(match(c1[,1],colnames(data)))]
data_c2 <- as.matrix(data)[,na.omit(match(c2[,1],colnames(data)))]


data_merge <- cbind(data_c1,data_c2)
data_merge=apply(data_merge,2,as.numeric)


rownames(data_merge) <- rownames(data)

##########

expr1 <- as.matrix(data_merge)


expr <- expr1



es.matrix <- gsva(
  expr,
  s.sets,
  min.sz = 10,
  max.sz = Inf,
  tau = 1,
  method = "gsva",
  abs.ranking = FALSE,
  verbose = TRUE,
  parallel.sz = 1
)

saveRDS(es.matrix,file="hallmark.rds")





n1 <- 1:dim(data_c1)[2]
#grep("male",seurat_obj@meta.data$sex)

n2 <- dim(data_c1)[2]+1:dim(data_c2)[2]

#grep("female",seurat_obj@meta.data$sex)

es.matrix.1 <-
  as.data.frame(es.matrix[, n1],
                row.names = row.names(es.matrix))
es.matrix.2 <-
  as.data.frame(es.matrix[, n2],
                row.names = row.names(es.matrix))


es.matrix.f <- cbind(es.matrix.1, es.matrix.2)

grouP <-
  c(rep("case", dim(es.matrix.1)[2]),
    rep("control", dim(es.matrix.2)[2]))

grouP <- as.factor( grouP)
design <- model.matrix(~ grouP + 0)



row.names(design)<-c(colnames(es.matrix.1), colnames(es.matrix.2))

comparE <-
  makeContrasts(grouPcase - grouPcontrol, levels = design)

fit <- lmFit(es.matrix, design)
fit2 <- contrasts.fit(fit, comparE)
fit3 <- eBayes(fit2)


diff <- topTable(fit3, coef = 1, number = dim(es.matrix)[1])

t_results <-
  as.data.frame(diff$t, row.names = rownames(es.matrix))
head(t_results)
colnames(t_results) <- c("t_value")




saveRDS(t_results,file="t_results_hallmark.rds")
write.csv(t_results,file="t_results_hallmark.csv")

#t_results <- readRDS("t_results_c5.go.bp_Elongating.rds")

library(ggplot2)

library(pheatmap)
rownames(t_results) <- gsub("HALLMARK_", "", rownames(t_results))
rownames(t_results) <- gsub("_", " ", rownames(t_results))
focus.cluster <- "t_value"
sub_t_results <- as.data.frame(t_results[, focus.cluster],
                               row.names = rownames(t_results))
sub_t_results$hallmark <- rownames(sub_t_results)
colnames(sub_t_results) <- c("t", "hallmark")

sub_t_results$hallmark = with(sub_t_results, reorder(hallmark, t))
sub_t_results$fill <- ""
sub_t_results[sub_t_results$t >= 2.58,]$fill <-
  "up"
# sub_t_results[sub_t_results$t <= -2.58,]$fill <-
#   "down"
# sub_t_results[abs(sub_t_results$t) < 2.58,]$fill <-
#   "no"
 sub_t_results$color <- ""
# sub_t_results[abs(sub_t_results$t) < 2.58,]$color <-
#   "n"
# sub_t_results[abs(sub_t_results$t) >= 2.58,]$color <-
#   "y"


sub_t_results <- sub_t_results[c(1:50),]

p <-
  ggplot(sub_t_results, aes(x = hallmark, y = t, fill = fill)) +
  geom_bar(stat = "identity", width = .8) +
  scale_fill_manual(
    values = c(
      #"down" = "#36648b",
      "up" = "#e94644"
      #"no" = "#cccccc"
    ),
    guide = F
  ) + ylim(-5,10)+
  geom_hline(
    yintercept = c(-2.58, 2.58),
    color = "white",
    linetype = "dotted",
    size = 0.5
  ) +
  coord_flip() +
  xlab("") +
  geom_text(
    data = subset(sub_t_results, t < 0),
    aes(
      x = hallmark,
      y = 0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 1.8,
    hjust = "inward"
  ) +geom_text(
    data = subset(sub_t_results, t > 0),
    aes(
      x = hallmark,
      y = -0.1,
      label = paste0(" ", hallmark),
      color = color
    ),
    size = 1.8,
    hjust = "outward"
  ) +
  scale_colour_manual(values = c("y" = "black"),
                      guide = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 0.5
    ),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.ticks.x = element_line(colour = "black", size = 0.5),
  )
#p
ggsave(
  filename = "hallmark.pdf",
  plot = p,
  height = 6.5,
  width =5
)

