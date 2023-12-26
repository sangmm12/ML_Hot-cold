setwd("D:/R/ML_PAAD/data1/AUCAtcga/DEGs&WGCNA/TPM_0.3_0.01_1_sec")


# modlist <- c("gene.pink","gene.yellow","gene.grey")
# WGCNA <- c()
# 
# #modu <- c("gene.royalblue")
# for(modu in modlist){
#   
#   library(data.table)
#   data <- read.table( paste("D:/R/ML_PAAD/data/AUCAtcga/WGCNA/TPM/",modu,".txt",sep=''),header = T)
#   
#   WGCNA <- c(WGCNA,data$x)
#   
# }
# 
# WGCNA <- unique(WGCNA)


library(data.table)
dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_black_hot_0.5_0.2.csv",header = T)
WGCNA_grown_hot <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_pink_cold_0.5_0.2.csv",header = T)
WGCNA_pink_cold <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_turquoise_cold_0.5_0.2.csv",header = T)
WGCNA_turquoise_cold <- unlist(dat[,1])

WGCNA_up <- WGCNA_grown_hot
WGCNA_down <- c(WGCNA_pink_cold,WGCNA_turquoise_cold)

DEGsup <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_up.csv",header = T)

DEGsup <- unique(DEGsup$x)

DEGsdown <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_down.csv",header = T)

DEGsdown <- unique(DEGsdown$x)


library(VennDiagram)

venn.diagram(list(WGCNA_black_hot=WGCNA_grown_hot,WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             alpha=c(0.5,0.5,0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180,180,180),
             cat.dist = c(0.03,0.01,0.05,0.05,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.2/WGCNA_DEGs.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)





venn.diagram(list(WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_down=DEGsdown), 
             fill=c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.2/WGCNA_DEGsadj0.01FC1_down0.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.2/WGCNA_DEGsadj0.01FC1_down.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5), 
             col = c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.2/WGCNA_DEGsadj0.01FC1_onlydown.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.2/WGCNA_DEGsadj0.01FC1_up.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup), 
             fill=c("#ED665D" ,"#FF9E4A" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.2/WGCNA_DEGsadj0.01FC1_onlyup.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C"  ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,0,0),
             cat.dist = c(0.1,0.1,0.1,0.1),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.2/WGCNA_0.5_0.2_DEGsadj0.01FC1_4.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)







library(data.table)
dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_black_hot_0.6_0.2.csv",header = T)
WGCNA_grown_hot <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_pink_cold_0.6_0.2.csv",header = T)
WGCNA_pink_cold <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_turquoise_cold_0.6_0.2.csv",header = T)
WGCNA_turquoise_cold <- unlist(dat[,1])

WGCNA_up <- WGCNA_grown_hot
WGCNA_down <- c(WGCNA_pink_cold,WGCNA_turquoise_cold)

DEGsup <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_up.csv",header = T)

DEGsup <- unique(DEGsup$x)

DEGsdown <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_down.csv",header = T)

DEGsdown <- unique(DEGsdown$x)


library(VennDiagram)

venn.diagram(list(WGCNA_black_hot=WGCNA_grown_hot,WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             alpha=c(0.5,0.5,0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180,180,180),
             cat.dist = c(0.03,0.01,0.05,0.05,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.2/WGCNA0.6_0.2_DEGs.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)





venn.diagram(list(WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_down=DEGsdown), 
             fill=c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.2/WGCNA0.6_0.2_DEGsadj0.01FC1_down0.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.2/WGCNA0.6_0.2_DEGsadj0.01FC1_down.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)

venn.diagram(list(CRGs=WGCNA_down,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,120),
             cat.dist = c(0.01,0.01,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.2/WGCNA0.6_0.2_DEGsadj0.01FC1_down1.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5), 
             col = c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.2/WGCNA0.6_0.2_DEGsadj0.01FC1_onlydown.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.2/WGCNA0.6_0.2_DEGsadj0.01FC1_up.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(HRGs=WGCNA_up,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.2/WGCNA0.6_0.2_DEGsadj0.01FC1_up1.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup), 
             fill=c("#ED665D" ,"#FF9E4A" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.2/WGCNA0.6_0.2_DEGsadj0.01FC1_onlyup.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)



venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C"  ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,0,0),
             cat.dist = c(0.1,0.1,0.1,0.1),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.2/WGCNA_0.6_0.2_DEGsadj0.01FC1_4.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)







library(data.table)
dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_black_hot_0.6_0.3.csv",header = T)
WGCNA_grown_hot <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_pink_cold_0.6_0.3.csv",header = T)
WGCNA_pink_cold <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_turquoise_cold_0.6_0.3.csv",header = T)
WGCNA_turquoise_cold <- unlist(dat[,1])

WGCNA_up <- WGCNA_grown_hot
WGCNA_down <- c(WGCNA_pink_cold,WGCNA_turquoise_cold)

DEGsup <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_up.csv",header = T)

DEGsup <- unique(DEGsup$x)

DEGsdown <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_down.csv",header = T)

DEGsdown <- unique(DEGsdown$x)


library(VennDiagram)

venn.diagram(list(WGCNA_black_hot=WGCNA_grown_hot,WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             alpha=c(0.5,0.5,0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180,180,180),
             cat.dist = c(0.03,0.01,0.05,0.05,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.3/WGCNA0.6_0.3_DEGs.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)





venn.diagram(list(WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_down=DEGsdown), 
             fill=c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.3/WGCNA0.6_0.3_DEGsadj0.01FC1_down0.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,90),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.3/WGCNA0.6_0.3_DEGsadj0.01FC1_down.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5), 
             col = c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.3/WGCNA0.6_0.3_DEGsadj0.01FC1_onlydown.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.3/WGCNA0.6_0.3_DEGsadj0.01FC1_up.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup), 
             fill=c("#ED665D" ,"#FF9E4A" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.3/WGCNA0.6_0.3_DEGsadj0.01FC1_onlyup.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)



venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C"  ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,0,0),
             cat.dist = c(0.1,0.1,0.1,0.1),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.6_0.3/WGCNA_0.6_0.3_DEGsadj0.01FC1_4.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)










library(data.table)
dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_black_hot_0.5_0.3.csv",header = T)
WGCNA_grown_hot <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_pink_cold_0.5_0.3.csv",header = T)
WGCNA_pink_cold <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_turquoise_cold_0.5_0.3.csv",header = T)
WGCNA_turquoise_cold <- unlist(dat[,1])

WGCNA_up <- WGCNA_grown_hot
WGCNA_down <- c(WGCNA_pink_cold,WGCNA_turquoise_cold)

DEGsup <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_up.csv",header = T)

DEGsup <- unique(DEGsup$x)

DEGsdown <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_down.csv",header = T)

DEGsdown <- unique(DEGsdown$x)


library(VennDiagram)

venn.diagram(list(WGCNA_black_hot=WGCNA_grown_hot,WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             alpha=c(0.5,0.5,0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180,180,180),
             cat.dist = c(0.03,0.01,0.05,0.05,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.3/WGCNA0.5_0.3_DEGs.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)





venn.diagram(list(WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_down=DEGsdown), 
             fill=c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.3/WGCNA0.5_0.3_DEGsadj0.01FC1_down0.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.3/WGCNA0.5_0.3_DEGsadj0.01FC1_down.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5), 
             col = c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.3/WGCNA0.5_0.3_DEGsadj0.01FC1_onlydown.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.3/WGCNA0.5_0.3_DEGsadj0.01FC1_up.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup), 
             fill=c("#ED665D" ,"#FF9E4A" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.3/WGCNA0.5_0.3_DEGsadj0.01FC1_onlyup.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)




venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C"  ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,0,0),
             cat.dist = c(0.1,0.1,0.1,0.1),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.3/WGCNA_0.5_0.3_DEGsadj0.01FC1_4.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)









library(data.table)
dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_black_hot_0.5_0.25.csv",header = T)
WGCNA_grown_hot <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_pink_cold_0.5_0.25.csv",header = T)
WGCNA_pink_cold <- unlist(dat[,1])

dat <- fread("D:/R/ML_PAAD/data1/AUCAtcga/WGCNA/TPM/all_name_turquoise_cold_0.5_0.25.csv",header = T)
WGCNA_turquoise_cold <- unlist(dat[,1])

WGCNA_up <- WGCNA_grown_hot
WGCNA_down <- c(WGCNA_pink_cold,WGCNA_turquoise_cold)

DEGsup <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_up.csv",header = T)

DEGsup <- unique(DEGsup$x)

DEGsdown <- read.csv("D:/R/ML_PAAD/data1/AUCAtcga/limma/adj0.01FC1/DEGs_down.csv",header = T)

DEGsdown <- unique(DEGsdown$x)


library(VennDiagram)

venn.diagram(list(WGCNA_black_hot=WGCNA_grown_hot,WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             alpha=c(0.5,0.5,0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ,"#ED665D","#AD8BC9"),
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180,180,180),
             cat.dist = c(0.03,0.01,0.05,0.05,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.25/WGCNA0.5_0.25_DEGs.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)





venn.diagram(list(WGCNA_pink_cold=WGCNA_pink_cold,WGCNA_turquoise_cold=WGCNA_turquoise_cold,DEGs_down=DEGsdown), 
             fill=c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#67BF5C" ,"#ED665D","#AD8BC9"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.25/WGCNA0.5_0.25_DEGsadj0.01FC1_down0.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5,0.5), 
             col = c("#729ECE","#FF9E4A","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.25/WGCNA0.5_0.25_DEGsadj0.01FC1_down.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             alpha=c(0.5,0.5), 
             col = c("#729ECE","#67BF5C" ),#,"#67BF5C" ,"#ED665D","#AD8BC9"
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.25/WGCNA0.5_0.25_DEGsadj0.01FC1_onlydown.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,180),
             cat.dist = c(0.03,0.01,0.05),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.25/WGCNA0.5_0.25_DEGsadj0.01FC1_up.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)


venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup), 
             fill=c("#ED665D" ,"#FF9E4A" ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A"),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180),
             cat.dist = c(0.03,0.01),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.25/WGCNA0.5_0.25_DEGsadj0.01FC1_onlyup.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)



venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup,WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
             fill=c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C"  ),#"#729ECE","#FF9E4A",
             alpha=c(0.5,0.5,0.5,0.5), 
             col = c("#ED665D" ,"#FF9E4A","#729ECE","#67BF5C" ),#"#729ECE","#FF9E4A",
             cex=1,
             # cat.pos = 10,
             cat.pos = c(180,180,0,0),
             cat.dist = c(0.1,0.1,0.1,0.1),
             
             cat.fontface=4, 
             # fontfamily=1,
             filename =paste("0.5_0.25/WGCNA_0.5_0.25_DEGsadj0.01FC1_4.tiff",sep=''),
             height = 1600, 
             width = 1600, resolution = 500)













# 
# 
# 
# 
# 
# all_nameupz <- names(which(table(c(WGCNA_up,DEGsup))==2))
# all_namedownz <- names(which(table(c(WGCNA_down,DEGsdown))==2))
# all_namez <- c(all_nameupz,all_namedownz)
# 
# if(dir.exists("up")==TRUE){
#   print("up")
# }else {
#   dir.create("up")
# }
# 
# # if(dir.exists("down")==TRUE){
# #   print(cancer)
# # }else {
# #   dir.create("down")
# # }
# 
# 
# 
# if(dir.exists("updown")==TRUE){
#   print("updown")
# }else {
#   dir.create("updown")
# }
# 
# # if(dir.exists("down")==TRUE){
# #   print(cancer)
# # }else {
# #   dir.create("down")
# # }
# 
# 
# 
# #################################################################################
# 
# library(dplyr)
# library(survival)
# library(data.table)
# 
# 
# tcga_list = c()
# 
# store_data = list()
# 
# #sample_name <- "tcga_PAAD"
# for(sample_name in c("TCGA_ICGC"))
# {
#   ####
#   #dat <- read.csv(paste('~/machine_learning/data/',sample_name,'/survival_out.csv',sep=''),header=T,row.name=1,check.names = F)
#   ####
#   
#   ###
#   #dat <- fread(paste('/home/sangmm/projects/ML_PAAD/single_cox/exp_',sample_name,'.csv',sep=''), header = T)
#   dat <- fread(paste('D:/R/ML_PAAD/data/alldata/exp_',sample_name,'.csv',sep=''), header = T)
#   dat <- as.data.frame(dat)
#   dat <- dat[!duplicated(dat[,c(1)]),]
#   rownames(dat) <- dat[,c(1)]
#   dat <- dat[,-c(1)]
#   dat <- na.omit(dat)
#   dat <- as.data.frame(t(dat))
#   aa <- colnames(dat)
#   ###
#   
#   #all_sur_data <- fread(paste('/home/sangmm/projects/ML_PAAD/single_cox/sur_',sample_name,'.csv',sep=''),header=T)
#   all_sur_data <- fread(paste('D:/R/ML_PAAD/data/alldata/sur_',sample_name,'.csv',sep=''),header=T)
#   all_sur_data <- as.data.frame(all_sur_data)
#   colnames(all_sur_data) <- c("SampleName","Time","Status")
#   output <- data.frame()
#   #gene_name <- "AC005304.1"
#   for(gene_name in all_namez)
#   {
#     dat1 <- dat[which(colnames(dat)==gene_name)]
#     
#     dat1$SampleName = rownames(dat1)
#     
#     erged_df <- merge(dat1, all_sur_data, by = "SampleName")
#     
#     erged_df <- erged_df[,-c(1)]
#     
#     cox_model <- coxph(Surv(Time,Status) ~ as.numeric(erged_df[,1]),data=erged_df)
#     
#     lower = round(summary(cox_model)$conf.int[,3],2)
#     upper = round(summary(cox_model)$conf.int[,4],2)
#     hr = round(summary(cox_model)$conf.int[,1],2)
#     p = summary(cox_model)$coefficients[,5]
#     # if(is.na(p))
#     # {
#     #   next
#     # }
#     # if(p>0.05||hr<1)
#     # {
#     #   next
#     # }
#     ci = paste(hr,'(',lower,'|',upper,')',sep='')
#     temp_data <- data.frame(genename = c(gene_name),pvalue = c(p),HR=c(hr),Lower=c(lower),Upper=c(upper),paste0(hr," (",lower,"-",lower,")"),CI=c(ci))
#     output <- rbind(output,temp_data)
#     eval(parse(text = paste('store_data$',sample_name,"<-append(store_data$",sample_name,",'",gene_name,"')",sep='')))
#   }
#   colnames(output) <- c("genename","pvalue","HR","Lower","Upper","Hazard Ratio(95%CI)","CI")
#   write.csv(output,file=paste("updown/single_cox",sample_name,".csv",sep=''),row.names = T)
#   
#   
#   # library(ggplot2)
#   # 
#   # if(sample_name=='tcga_PAAD')
#   # {
#   #   data <- data.frame(
#   #     group = output$genename,
#   #     HR = output$HR,
#   #     p_value = output$pvalue,
#   #     lower=output$Lower,
#   #     upper = output$Upper
#   #   )
#   #   tcga_list <- output$HR
#   #   data <- data[order(data$HR),]
#   # } else {
#   #   data <- data.frame(
#   #     group = output$genename,
#   #     HR = output$HR,
#   #     tcga = tcga_list,
#   #     p_value = output$pvalue,
#   #     lower=output$Lower,
#   #     upper = output$Upper
#   #   )
#   #   data <- data[order(data$tcga),]
#   # }
#   # data$lower[which(data$HR>=1)] = 'Risky'
#   # data$lower[which(data$HR<1)] = 'Protect'
#   # data$lower[which(data$p_value>=0.05)] = 'Not'
#   # eval(parse(text = paste('store_data$',sample_name,"<-data$lower",sep='')))
#   # 
#   
# }
# 
# save(store_data,file='updown/store.Rdata')
# 
# 
# aaaaup <- output[match(all_nameupz,output$genename),]
# write.csv(aaaaup,file =paste("updown/single_cox",sample_name,"_up.csv",sep=''),quote=T)
# 
# aaaadown <- output[match(all_namedownz,output$genename),]
# write.csv(aaaadown,file =paste("updown/single_cox",sample_name,"_down.csv",sep=''),quote=T)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# modlist <- c("gene.grey")
# WGCNA_up <- c()
# 
# #modu <- c("gene.royalblue")
# for(modu in modlist){
#   
#   library(data.table)
#   data <- read.table( paste("D:/R/ML_PAAD/data/AUCAtcga/WGCNA/TPM/",modu,".txt",sep=''),header = T)
#   
#   WGCNA_up <- c(WGCNA_up,data$x)
#   
# }
# 
# WGCNA <- unique(WGCNA)
# 
# DEGsup <- read.csv("D:/R/ML_PAAD/data/AUCAtcga/limma/adj0.01FC1/DEGs_up.csv",header = T)
# 
# DEGsup <- unique(DEGsup$x)
# 
# 
# library(VennDiagram)
# 
# venn.diagram(list(WGCNA_up=WGCNA_up,DEGs_up=DEGsup), 
#              fill=c("#729ECE","#FF9E4A"), #,"#67BF5C","#ED665D","#AD8BC9"
#              alpha=c(0.5,0.5), 
#              col = c("#729ECE","#FF9E4A"),
#              cex=1,
#              # cat.pos = 10,
#              cat.pos = c(180,180),
#              cat.dist = c(0.03,0.01),
#              
#              cat.fontface=4, 
#              # fontfamily=1,
#              filename =paste("WGCNA_DEGs_up.tiff",sep=''),
#              height = 1600, 
#              width = 1600, resolution = 500)
# 
# 
# 
# 
# 
# modlist <- c("gene.pink","gene.yellow")
# WGCNA_down <- c()
# 
# #modu <- c("gene.royalblue")
# for(modu in modlist){
#   
#   library(data.table)
#   data <- read.table( paste("D:/R/ML_PAAD/data/AUCAtcga/WGCNA/TPM/",modu,".txt",sep=''),header = T)
#   
#   WGCNA_down <- c(WGCNA_down,data$x)
#   
# }
# 
# 
# DEGsdown <- read.csv("D:/R/ML_PAAD/data/AUCAtcga/limma/adj0.01FC1/DEGs_down.csv",header = T)
# 
# DEGsdown <- unique(DEGsdown$x)
# 
# 
# library(VennDiagram)
# 
# venn.diagram(list(WGCNA_down=WGCNA_down,DEGs_down=DEGsdown), 
#              fill=c("#729ECE","#67BF5C"), #,"#67BF5C","#ED665D","#AD8BC9"
#              alpha=c(0.5,0.5), 
#              col = c("#729ECE","#67BF5C"),
#              cex=1,
#              # cat.pos = 10,
#              cat.pos = c(180,180),
#              cat.dist = c(0.03,0.05),
#              
#              cat.fontface=4, 
#              # fontfamily=1,
#              filename =paste("WGCNA_DEGs_down.tiff",sep=''),
#              height = 1600, 
#              width = 1600, resolution = 500)
# 
# 
# 
# all_nameup <- names(which(table(c(WGCNA_up,DEGsup))==2))
# 
# 
# library(dplyr)
# library(survival)
# library(data.table)
# 
# 
# tcga_list = c()
# 
# store_data = list()
# 
# #sample_name <- "tcga_PAAD"
# for(sample_name in c("TCGA_ICGC"))
# {
#   ####
#   #dat <- read.csv(paste('~/machine_learning/data/',sample_name,'/survival_out.csv',sep=''),header=T,row.name=1,check.names = F)
#   ####
#   
#   ###
#   #dat <- fread(paste('/home/sangmm/projects/ML_PAAD/single_cox/exp_',sample_name,'.csv',sep=''), header = T)
#   dat <- fread(paste('D:/R/ML_PAAD/data/alldata/exp_',sample_name,'.csv',sep=''), header = T)
#   dat <- as.data.frame(dat)
#   dat <- dat[!duplicated(dat[,c(1)]),]
#   rownames(dat) <- dat[,c(1)]
#   dat <- dat[,-c(1)]
#   dat <- na.omit(dat)
#   dat <- as.data.frame(t(dat))
#   
#   ###
#   
#   #all_sur_data <- fread(paste('/home/sangmm/projects/ML_PAAD/single_cox/sur_',sample_name,'.csv',sep=''),header=T)
#   all_sur_data <- fread(paste('D:/R/ML_PAAD/data/alldata/sur_',sample_name,'.csv',sep=''),header=T)
#   all_sur_data <- as.data.frame(all_sur_data)
#   colnames(all_sur_data) <- c("SampleName","Time","Status")
#   output <- data.frame()
#   #gene_name <- "AC005304.1"
#   for(gene_name in all_nameup)
#   {
#     dat1 <- dat[which(colnames(dat)==gene_name)]
#     
#     dat1$SampleName = rownames(dat1)
#     
#     erged_df <- merge(dat1, all_sur_data, by = "SampleName")
#     
#     erged_df <- erged_df[,-c(1)]
#     
#     cox_model <- coxph(Surv(Time,Status) ~ as.numeric(erged_df[,1]),data=erged_df)
#     
#     lower = round(summary(cox_model)$conf.int[,3],2)
#     upper = round(summary(cox_model)$conf.int[,4],2)
#     hr = round(summary(cox_model)$conf.int[,1],2)
#     p = summary(cox_model)$coefficients[,5]
#     # if(is.na(p))
#     # {
#     #   next
#     # }
#     # if(p>0.05||hr<1)
#     # {
#     #   next
#     # }
#     ci = paste(hr,'(',lower,'|',upper,')',sep='')
#     temp_data <- data.frame(genename = c(gene_name),pvalue = c(p),HR=c(hr),Lower=c(lower),Upper=c(upper),paste0(hr," (",lower,"-",lower,")"),CI=c(ci))
#     output <- rbind(output,temp_data)
#     eval(parse(text = paste('store_data$',sample_name,"<-append(store_data$",sample_name,",'",gene_name,"')",sep='')))
#   }
#   colnames(output) <- c("genename","pvalue","HR","Lower","Upper","Hazard Ratio(95%CI)","CI")
#   write.csv(output,file=paste("up/single_cox",sample_name,"_up.csv",sep=''),row.names = T)
# }