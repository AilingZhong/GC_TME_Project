## Bulk RNAseq Analyses

### Part1. showing the heatmap of TM,TMP,and TMPC  specific genes 

```R
library(Seurat)
library(ggplot2)
library(pheatmap)
require(RColorBrewer)
source("MyBestFunction_scRNA.R")

spcific_genes_in_9_subcotous <- read.csv("./Subcutaneous_RNAseq.csv")
spcific_genes_in_9_subcotous <- na.omit(spcific_genes_in_9_subcotous)
TM_group <- subset(spcific_genes_in_9_subcotous,TM_VS_others_pvalue <=0.05 & TM_VS_others_log2FoldChange >0 )
TMP_group <- subset(spcific_genes_in_9_subcotous,TMP_VS_others_pvalue <=0.05 &  TMP_VS_others_log2FoldChange >0)
TMPCd_group <- subset(spcific_genes_in_9_subcotous,TMPC_VS_others_pvalue <=0.05 &  TMPC_VS_others_log2FoldChange >0)

genes_sel <- c("Abcc3","Ace","Adrb2","Agt","Areg","Muc13","Muc5ac","Ntn1","Rb1","Sftpd","Sulf1","Tff2")
rownames(TM_group) <- TM_group$Symbol
TM_group <- TM_group[order(-TM_group$TM_VS_others_log2FoldChange),]
TM_group <- TM_group[,c("DESeq2_Sub_TM_1","DESeq2_Sub_TM_2","DESeq2_Sub_TM_3","DESeq2_Sub_TMP_1","DESeq2_Sub_TMP_2","DESeq2_Sub_TMP_3","DESeq2_Sub_TMPC_1","DESeq2_Sub_TMPC_2","DESeq2_Sub_TMPC_3")]
tmp_data <- log(TM_group+1,2)
chonglai_zscore_1 <- t(apply(tmp_data, 1, function(x) (x-mean(x))/sd(x)))
range(chonglai_zscore_1)
chonglai_zscore_1[chonglai_zscore_1 > 1] <- 1
chonglai_zscore_1[chonglai_zscore_1 < -1] <- -1
SeuratObject <- CreateSeuratObject(counts = chonglai_zscore_1, project = "TM")
gene <- rownames(SeuratObject)
SeuratObject@meta.data$group <- rownames(SeuratObject@meta.data)
pdf("Fig3_TM_gene_module.pdf",width=10,height=3)
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident",gene = gene,all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),
  min_and_max_cut=2,show_row_names=FALSE,mark_gene=genes_sel,label_size=0,scale = FALSE)
dev.off()          

genes_sel <- c("Gdf5","Adamts8","Amhr2","Col11a2","Rell2","Mgp","Col10a1","Lgr6","Lox","Pou3f4","Ahsg","Gpr173","Igf2")
rownames(TMP_group) <- TMP_group$Symbol
TMP_group <- TMP_group[order(-TMP_group$TMP_VS_others_log2FoldChange),]
TMP_group <- TMP_group[,c("DESeq2_Sub_TM_1","DESeq2_Sub_TM_2","DESeq2_Sub_TM_3","DESeq2_Sub_TMP_1","DESeq2_Sub_TMP_2","DESeq2_Sub_TMP_3","DESeq2_Sub_TMPC_1","DESeq2_Sub_TMPC_2","DESeq2_Sub_TMPC_3")]
tmp_data <- log(TMP_group+1,2)
chonglai_zscore_1 <- t(apply(tmp_data, 1, function(x) (x-mean(x))/sd(x)))
range(chonglai_zscore_1)
chonglai_zscore_1[chonglai_zscore_1 > 1] <- 1
chonglai_zscore_1[chonglai_zscore_1 < -1] <- -1
SeuratObject <- CreateSeuratObject(counts = chonglai_zscore_1, project = "TMP")
gene <- rownames(SeuratObject)
SeuratObject@meta.data$group <- rownames(SeuratObject@meta.data)
pdf("Fig3_TMP_gene_module.pdf",width=10,height=3)
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident",gene = gene,all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),
  min_and_max_cut=2,show_row_names=FALSE,mark_gene=genes_sel,label_size=0,scale = FALSE)
dev.off()          

genes_sel <- c("Twist1","Top2a","Mcm2","Mcm10","Lef1","Isl1","Ezh2","Cdc45","Brca1","Actn3")
rownames(TMPCd_group) <- TMPCd_group$Symbol
TMPCd_group <- TMPCd_group[order(-TMPCd_group$TMPC_VS_others_log2FoldChange),]
TMPCd_group <- TMPCd_group[,c("DESeq2_Sub_TM_1","DESeq2_Sub_TM_2","DESeq2_Sub_TM_3","DESeq2_Sub_TMP_1","DESeq2_Sub_TMP_2","DESeq2_Sub_TMP_3","DESeq2_Sub_TMPC_1","DESeq2_Sub_TMPC_2","DESeq2_Sub_TMPC_3")]
tmp_data <- log(TMPCd_group+1,2)
chonglai_zscore_1 <- t(apply(tmp_data, 1, function(x) (x-mean(x))/sd(x)))
range(chonglai_zscore_1)
chonglai_zscore_1[chonglai_zscore_1 > 1] <- 1
chonglai_zscore_1[chonglai_zscore_1 < -1] <- -1
range(chonglai_zscore_1)
chonglai_zscore_1[chonglai_zscore_1 > 1] <- 1
chonglai_zscore_1[chonglai_zscore_1 < -1] <- -1
SeuratObject <- CreateSeuratObject(counts = chonglai_zscore_1, project = "TMPC")
gene <- rownames(SeuratObject)
SeuratObject@meta.data$group <- rownames(SeuratObject@meta.data)
pdf("Fig3_TMPC_gene_module.pdf",width=10,height=3)
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident",gene = gene,all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=colorRampPalette(brewer.pal(10, "RdBu"))(101),
  min_and_max_cut=2,show_row_names=FALSE,mark_gene=genes_sel,label_size=0,scale = FALSE)
dev.off()       
```

![](Bulk_RNAseq.assets/Fig3_TM_gene_module.png)

![](Bulk_RNAseq.assets/Fig3_TMP_gene_module-01.png)

![Fig3_TMPC_gene_module](Bulk_RNAseq.assets/Fig3_TMPC_gene_module.png)





### Part2. showing the GO enrichment results of TM,TMP,and TMPC specific genes 

```R
spcific_genes_in_9_subcotous <- read.csv("./Subcutaneous_RNAseq.csv")
spcific_genes_in_9_subcotous <- na.omit(spcific_genes_in_9_subcotous)
sig_high_TM <- subset(spcific_genes_in_9_subcotous,TM_pvalue <=0.05 & TM_VS_others_log2FoldChange >1 )
sig_high_TMP <- subset(spcific_genes_in_9_subcotous,TMP_pvalue <=0.05 &  TMP_VS_others_log2FoldChange >1)
sig_high_TMPCd <- subset(spcific_genes_in_9_subcotous,TMPCd_pvalue <=0.05 &  TMPCd_VS_others_log2FoldChange >1)
GO_p005_TM <- enrichGO(gene = as.character(sig_high_TM$ENTREZID), 
             OrgDb = org.Mm.eg.db,
        ont = "all", 
                 pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
GO_p005_TMP <- enrichGO(gene = as.character(sig_high_TMP$ENTREZID), 
             OrgDb = org.Mm.eg.db,
        ont = "all", 
                 pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
GO_p005_TMPCd <- enrichGO(gene = as.character(sig_high_TMPCd$ENTREZID), 
             OrgDb = org.Mm.eg.db,
        ont = "all", 
                 pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)

res_GO_TMPCd <- read.csv("./results_GO_TMPCd.csv")
res_GO_TMP <- read.csv("./results_GO_TMP.csv")
res_GO_TM <- read.csv("./results_GO_TM.csv")
library(ggpubr)
res_GO_TM <- as.data.frame(res_GO_TM)
res_GO_TM$log10_p.adjust <- -log(res_GO_TM$p.adjust,10)
res_GO_TM_10 <- head(res_GO_TM,10)
p1 <- ggbarplot(res_GO_TM_10, 
  x = "Description", 
  y = "log10_p.adjust",
  color = "#5B8FCF",            # Set bar border colors to white
  fill ="#5B8FCF",
  sort.val = "asc",          # Sort the value in dscending order
  x.text.angle = 90,           # Rotate vertically x axis texts
  rotate = TRUE,
  title="atac methylation overlap")
ggsave(p1,file="Fig3_GO_TM.png",width =9, height = 2.5,dpi=1080)

res_GO_TMP <- as.data.frame(res_GO_TMP)
res_GO_TMP$log10_p.adjust <- -log(res_GO_TMP$p.adjust,10)
res_GO_TMP_10 <- head(res_GO_TMP,10)
p1 <- ggbarplot(res_GO_TMP_10, 
  x = "Description", 
  y = "log10_p.adjust",
  color = "#5B8FCF",            # Set bar border colors to white
  fill ="#5B8FCF",
  sort.val = "asc",          # Sort the value in dscending order
  x.text.angle = 90,           # Rotate vertically x axis texts
  rotate = TRUE,
  title="atac methylation overlap")
ggsave(p1,file="Fig3_GO_TMP.png",width =9, height = 2.5,dpi=1080)

res_GO_TMPCd <- as.data.frame(res_GO_TMPCd)
res_GO_TMPCd$log10_p.adjust <- -log(res_GO_TMPCd$p.adjust,10)
res_GO_TMPCd_10 <- head(res_GO_TMPCd,10)
p1 <- ggbarplot(res_GO_TMPCd_10, 
  x = "Description", 
  y = "log10_p.adjust",
  color = "#5B8FCF",            # Set bar border colors to white
  fill ="#5B8FCF",
  sort.val = "asc",          # Sort the value in dscending order
  x.text.angle = 90,           # Rotate vertically x axis texts
  rotate = TRUE,
  title="atac methylation overlap")
ggsave(p1,file="Fig3_GO_TMPCd.png",width =9, height = 2.5,dpi=1080)
```

![Fig3_GO_TM](Bulk_RNAseq.assets/Fig3_GO_TM.png)

![Fig3_GO_TMP](Bulk_RNAseq.assets/Fig3_GO_TMP.png)

![Fig3_GO_TMPCd](Bulk_RNAseq.assets/Fig3_GO_TMPCd.png)



### Part3. The cell cycle and poorly differentiated scores visualization

```R
library(ggplot2)
library(ggpubr)
library(car)
library(plyr)

spcific_genes_in_9_subcotous <- read.csv("./Subcutaneous_RNAseq.csv")
spcific_genes_in_9_subcotous <- na.omit(spcific_genes_in_9_subcotous)
CellCycle_Undifferentiated_signatures <- read.csv("./CellCycle_Undifferentiated_signatures.csv")
CellCycle_sig <- subset(CellCycle_Undifferentiated_signatures,Signature=="CellCycle")
CellCycle_sig_Normalized_data <- merge(CellCycle_sig,spcific_genes_in_9_subcotous,by="Symbol")
CellCycle_sig_Normalized_data <- CellCycle_sig_Normalized_data[!duplicated(CellCycle_sig_Normalized_data$Symbol),]
rownames(CellCycle_sig_Normalized_data) <- CellCycle_sig_Normalized_data$Symbol
CellCycle_sig_Normalized_data <- CellCycle_sig_Normalized_data[,c("DESeq2_Sub_TM_1","DESeq2_Sub_TM_2","DESeq2_Sub_TM_3","DESeq2_Sub_TMP_1","DESeq2_Sub_TMP_2","DESeq2_Sub_TMP_3","DESeq2_Sub_TMPC_1","DESeq2_Sub_TMPC_2","DESeq2_Sub_TMPC_3")]
log2 <- log(CellCycle_sig_Normalized_data+1,2)
col_mean_log2 <- apply(log2,2,mean)
col_mean <- data.frame(col_mean_log2)
col_mean$group <- c("TM","TM","TM","TMP","TMP","TMP","TMPCd","TMPCd","TMPCd")

my_comparisons <- list(c("TM", "TMP"),c("TMP","TMPCd"),c("TM","TMPCd"))
ff <- ggboxplot(col_mean, x = "group", y = "col_mean_log2",
               color = "group",
                add = "jitter") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="t.test")+ylim(9,11)
ggsave(ff,file="Fig3_Cellcycle_boxplot.png",width=5,height=7,dpi=1080)
```

![Fig3_Cellcycle_boxplot](Bulk_RNAseq.assets/Fig3_Cellcycle_boxplot.png)

```R
Poorly_Differentiated_sig <- subset(CellCycle_Undifferentiated_signatures,Signature=="Poorly_Differentiated")
Poorly_Differentiated_sig_Normalized_data <- merge(Poorly_Differentiated_sig,spcific_genes_in_9_subcotous,by="Symbol")
Poorly_Differentiated_sig_Normalized_data <- Poorly_Differentiated_sig_Normalized_data[!duplicated(Poorly_Differentiated_sig_Normalized_data$Symbol),]
rownames(Poorly_Differentiated_sig_Normalized_data) <- Poorly_Differentiated_sig_Normalized_data$Symbol
Poorly_Differentiated_sig_Normalized_data <- Poorly_Differentiated_sig_Normalized_data[,c("DESeq2_Sub_TM_1","DESeq2_Sub_TM_2","DESeq2_Sub_TM_3","DESeq2_Sub_TMP_1","DESeq2_Sub_TMP_2","DESeq2_Sub_TMP_3","DESeq2_Sub_TMPC_1","DESeq2_Sub_TMPC_2","DESeq2_Sub_TMPC_3")]
log2 <- log(Poorly_Differentiated_sig_Normalized_data+1,2)
col_mean_log2 <- apply(log2,2,mean)
col_mean <- data.frame(col_mean_log2)
col_mean$group <- c("TM","TM","TM","TMP","TMP","TMP","TMPCd","TMPCd","TMPCd")
library(ggplot2)
library(ggpubr)
library(car)
library(plyr)
my_comparisons <- list(c("TM", "TMP"),c("TMP","TMPCd"),c("TM","TMPCd"))
ff <- ggboxplot(col_mean, x = "group", y = "col_mean_log2",
               color = "group",
                add = "jitter") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="t.test",paired=FALSE)+ylim(9,11.8)
ggsave(ff,file="Fig3_Poorly_Differentiated_boxplot.png",width=5,height=7,dpi=1080)
```

![Fig3_Poorly_Differentiated_boxplot](Bulk_RNAseq.assets/Fig3_Poorly_Differentiated_boxplot.png)



### Part4. The transcriptome differences between subcutaneous and orthotopic TMPC tumors

```R
insitu_VS_sub_TMPCd_RNAseq <- read.csv("./insitu_VS_sub_TMPCd_DEseq2normalized_allsummry.csv")
insitu_VS_sub_TMPCd_RNAseq <- na.omit(insitu_VS_sub_TMPCd_RNAseq)
insitu_VS_sub_TMPCd_RNAseq <- subset(insitu_VS_sub_TMPCd_RNAseq,pvalue < 0.05 & abs(log2FoldChange) >=1)
rownames(insitu_VS_sub_TMPCd_RNAseq) <- insitu_VS_sub_TMPCd_RNAseq$X
insitu_VS_sub_TMPCd_RNAseq <- insitu_VS_sub_TMPCd_RNAseq[order(-insitu_VS_sub_TMPCd_RNAseq$log2FoldChange),]
insitu_VS_sub_TMPCd_RNAseq <- insitu_VS_sub_TMPCd_RNAseq[,c("DESeq2_insitu_TMPCd_1", "DESeq2_insitu_TMPCd_2","DESeq2_insitu_TMPCd_3", "DESeq2_sub_TMPCd_1", "DESeq2_sub_TMPCd_2", "DESeq2_sub_TMPCd_3")]

library(Seurat)
genes_sel <- c("Ccr3","Sfrp1","Dcn","Eln","Cxcl12","Fgf2","Cd74","Mmp2","Wnt5a","Hoxa7","Krt19","Pik3cd","Csf1r","Pik3cg","Ccr2","Itgb4","Krt7","Tgfbr2","Egfr","Cxcl10","Itga6")
tmp_data <- log(insitu_VS_sub_TMPCd_RNAseq+1,2)
chonglai_zscore_1 <- t(apply(tmp_data, 1, function(x) (x-mean(x))/sd(x)))
range(chonglai_zscore_1)
chonglai_zscore_1[chonglai_zscore_1 > 1] <- 1
chonglai_zscore_1[chonglai_zscore_1 < -1] <- -1
SeuratObject <- CreateSeuratObject(counts = chonglai_zscore_1, project = "")
gene <- rownames(SeuratObject)

SeuratObject@meta.data$group <- rownames(SeuratObject@meta.data)
ff <- XY_heatmap(seurat_obj=SeuratObject,group="orig.ident",genes=gene,all_num=FALSE,new_names=NULL,labels_rot=90,assay_sel="RNA",color=c("#74a9cf","#a6bddb","#d0d1e6","#ece7f2","#fcbba1","#fb6a4a","#ef3b2c","#cb181d","#a50f15"),min_and_max_cut=2,show_row_names=FALSE,mark_gene=genes_sel,label_size=0,scale = FALSE)
ggsave(ff,file="insituTMPC_VS_subTM_gene_module.pdf",width=6,height=6,dpi=1080)
```

![](Bulk_RNAseq.assets/Fig3_insituTMPC_VS_subTM_gene_module.png)

```R
insitu_VS_Sub_TMPCd <- read.csv("./insitu_TMPCd_VS_Sub_TMPCd_RNAseq.csv")
geneList <- insitu_VS_Sub_TMPCd[,c("HGNC.symbol","log2FoldChange")]
geneList <- geneList[!duplicated(geneList$HGNC.symbol), ]
geneList <- geneList[order(-geneList$log2FoldChange),]
geneList_expr <- geneList[,2]
names(geneList_expr) <- as.character(geneList[,1])
gmtfile_c5 <- read.gmt("./c5.all.v6.1.symbols.gmt")
gmtfile_h <- read.gmt("./h.all.v6.1.symbols.gmt")
gmtfile_c2 <- read.gmt("./c2.all.v6.1.symbols.gmt")
GSEA_c5 <- GSEA(geneList_expr,  pvalueCutoff = 1, seed=10,TERM2GENE=gmtfile_c5, verbose=FALSE, minGSSize = 5,maxGSSize = 500, nPerm = 13500)
GSEA_h <- GSEA(geneList_expr, pvalueCutoff = 1, seed=10,TERM2GENE=gmtfile_h, verbose=FALSE,minGSSize = 5,maxGSSize = 500,nPerm = 13500)
GSEA_c2 <- GSEA(geneList_expr, pvalueCutoff = 1, seed=10,TERM2GENE=gmtfile_c2, verbose=FALSE,minGSSize = 5,maxGSSize = 500,nPerm = 13500)
c5 <- data.frame(GSEA_c5)
h <- data.frame(GSEA_h)
c2 <- data.frame(GSEA_c2)
mouse_GSEA_hc2c5 <- rbind(h,c2,c5)

rownames(insitu_VS_Sub_GSEA_hc2c5) <- mouse_GSEA_hc2c5$X
insitu_VS_Sub_GSEA_hc2c5 <- insitu_VS_Sub_GSEA_hc2c5[,-1]
select_insitu_VS_Sub_GSEA_hc2c5 <- subset(insitu_VS_Sub_GSEA_hc2c5,ID=="GO_CYTOKINE_RECEPTOR_ACTIVITY" |ID=="GO_CELL_CHEMOTAXIS" |ID=="GO_CILIUM_MOVEMENT" | ID=="KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"|ID=="GO_EPITHELIAL_CILIUM_MOVEMENT" | ID=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" | ID=="WU_CELL_MIGRATION"|ID=="TAVAZOIE_METASTASIS"|ID=="REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS")
GSEA_h@result <- data.frame(select_insitu_VS_Sub_GSEA_hc2c5)
library(ggplot2)
pdf("renew_selected_insitu_VS_Sub_GSEA.pdf",width=10,height=6)
ridgeplot(GSEA_h,core_enrichment = TRUE,showCategory =10,fill = "p.adjust")(50))+xlim(-0.5,8)+scale_color_gradient2(low="#fcbba1",mid="#fb6a4a",high="#67000d")
dev.off()
```

![](assets/Fig3_renew_selected_insitu_VS_Sub_GSEA-16566010123145.png)

### Part 5.  The transcriptome similarities of TM, TMP&TMPC mouse tumor tissues, human normal, and tumor tissues. 

```R
DEseq_Normal <- readRDS("./DEseq_normal.rds")
library(DESeq2)
res_Normal <- results(DEseq_Normal, contrast=c("group_2","Normal","non_Normal"))
names(res_Normal) <- c("Nor_VS_Tumor_baseMean","Nor_VS_Tumor_log2FoldChange","Nor_VS_Tumor_lfcSE","Nor_VS_Tumor_stat","Nor_VS_Tumor_pvalue","Nor_VS_Tumor_padj")

res_Normal <- data.frame(res_Normal)
ENSEMBL <- rownames(res_Normal)
res_Normal$symbol <- mapIds(x = org.Hs.eg.db,
  keys = as.character(ENSEMBL),
  keytype ="ENSEMBL",
  column ="SYMBOL",
  multiVals="first")

Normal_sig <- subset(res_Normal,Nor_VS_Tumor_pvalue <=0.05 & Nor_VS_Tumor_log2FoldChange >0 )
Normal_sig <- Normal_sig[order(Normal_sig$Nor_VS_Tumor_log2FoldChange,decreasing=TRUE),]
Normal_sig <- na.omit(Normal_sig)
dim(Normal_sig)
Normal_sig <- data.frame(Normal_sig[,7])

Tumor_sig <- subset(res_Normal,Nor_VS_Tumor_pvalue <=0.05 & Nor_VS_Tumor_log2FoldChange < 0 )
Tumor_sig <- Tumor_sig[order(Tumor_sig$Nor_VS_Tumor_log2FoldChange,decreasing=FALSE),]
Tumor_sig <- na.omit(Tumor_sig)
dim(Tumor_sig)
Tumor_sig <- data.frame(Tumor_sig[,7])

names(Normal_sig) <- "symbol"
names(Tumor_sig) <- "symbol"

spcific_genes_in_9_subcotous <- read.csv("./all_mouse_spcific_genes_in_9_subcotous_to_human.csv")
TM_sig <- subset(spcific_genes_in_9_subcotous,TM_pvalue <=0.05 & TM_log2FoldChange > 0)

TMPandCd_sig <- subset(spcific_genes_in_9_subcotous,TM_pvalue <=0.05 & TM_log2FoldChange < 0)

dim(TM_sig)
dim(TMPandCd_sig)

TM_sig <- TM_sig[order(TM_sig$TM_log2FoldChange,decreasing=TRUE),]
TM_sig <- data.frame(TM_sig[,32])

TMPandCd_sig <- TMPandCd_sig[order(TMPandCd_sig$TM_log2FoldChange,decreasing=TRUE),]
TMPandCd_sig <- data.frame(TMPandCd_sig[,32])

names(TM_sig) <- "symbol"
names(TMPandCd_sig) <- "symbol"

length(unique(intersect(Normal_sig$symbol,TM_sig$symbol)))
length(unique(intersect(Normal_sig$symbol,TMPandCd_sig$symbol)))
length(unique(intersect(Tumor_sig$symbol,TM_sig$symbol)))
length(unique(intersect(Tumor_sig$symbol,TMPandCd_sig$symbol)))

Normal_TM_complement_1 <- length(Normal_sig$symbol)-length(unique(intersect(Normal_sig$symbol,TM_sig$symbol)))
TM_Normal_complement_2 <- length(TM_sig$symbol)-length(unique(intersect(Normal_sig$symbol,TM_sig$symbol)))
Tumor_TM_complement_1 <- length(Tumor_sig$symbol)-length(unique(intersect(Tumor_sig$symbol,TM_sig$symbol)))
TM_Tumor_complement_2 <- length(TM_sig$symbol)-length(unique(intersect(Tumor_sig$symbol,TM_sig$symbol)))

Normal_TMPandCd_complement_1 <- length(Normal_sig$symbol)-length(unique(intersect(Normal_sig$symbol,TMPandCd_sig$symbol)))
TMPandCd_Normal_complement_2 <- length(TMPandCd_sig$symbol)-length(unique(intersect(Normal_sig$symbol,TMPandCd_sig$symbol)))
Tumor_TMPandCd_complement_1 <- length(Tumor_sig$symbol)-length(unique(intersect(Tumor_sig$symbol,TMPandCd_sig$symbol)))
TMPandCd_Tumor_complement_2 <- length(TMPandCd_sig$symbol)-length(unique(intersect(Tumor_sig$symbol,TMPandCd_sig$symbol)))

library(ggplot2)
library(eulerr)
fit1 <- euler(c("Normal_sig" = 5560, "TM_sig" = 1815,
                "Normal_sig&TM_sig" = 1351))
ff <- plot(fit1,quantities = TRUE )
ggsave("Normal_TM_ol.png", plot=ff,width = 5,height = 5,dpi=1080)

fit1 <- euler(c("Normal_sig" = 6182, "TMPandCd" = 2583,
                "Normal_sig&TMPandCd" = 729))
ff <- plot(fit1,quantities = TRUE )
ggsave("Normal_TMPandCd_ol.png", plot=ff,width = 5,height = 5,dpi=1080)

fit1 <- euler(c("Tumor_sig" = 7124, "TM_sig" = 2499,
                "Tumor_sig&TM_sig" = 667))
ff <- plot(fit1,quantities = TRUE )
ggsave("Tumor_TM_ol.png", plot=ff,width = 5,height = 5,dpi=1080)

fit1 <- euler(c("Tumor_sig" = 6293, "TMPandCd" = 1814,
                "Tumor_sig&TMPandCd" = 1498))
ff <- plot(fit1,quantities = TRUE )
ggsave("Tumor_TMPandCd_ol.png", plot=ff,width = 5,height = 5,dpi=1080)
```

![](Bulk_RNAseq.assets/FigS3_Venn_Diagram.png)



```R
require(gmp)
enrich_pvalue <- function(N, A, B, k)
{
    m <- A + k
    n <- B + k
    i <- k:min(m,n)

    as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}

Normal_TM_pvalue <- enrich_pvalue(24421,Normal_TM_complement_1,TM_Normal_complement_2,length(unique(intersect(Normal_sig$symbol,TM_sig$symbol))))
Tumor_TM_pvalue <- enrich_pvalue(24421,Tumor_TM_complement_1,TM_Tumor_complement_2,length(unique(intersect(Tumor_sig$symbol,TM_sig$symbol))))
Normal_TMPandCd_pvalue <- enrich_pvalue(24421,Normal_TMPandCd_complement_1,TMPandCd_Normal_complement_2,length(unique(intersect(Normal_sig$symbol,TMPandCd_sig$symbol))))
Tumor_TMPandCd_pvalue <- enrich_pvalue(24421,Tumor_TMPandCd_complement_1,TMPandCd_Tumor_complement_2,length(unique(intersect(Tumor_sig$symbol,TMPandCd_sig$symbol))))

Normal_TM_ratio <- length(unique(intersect(Normal_sig$symbol,TM_sig$symbol)))/length(TM_sig$symbol)*100
Tumor_TM_ratio <- length(unique(intersect(Tumor_sig$symbol,TM_sig$symbol)))/length(TM_sig$symbol)*100
Normal_TMPandCd_ratio <- length(unique(intersect(Normal_sig$symbol,TMPandCd_sig$symbol)))/length(TMPandCd_sig$symbol)*100
Tumor_TMPandCd_ratio <- length(unique(intersect(Tumor_sig$symbol,TMPandCd_sig$symbol)))/length(TMPandCd_sig$symbol)*100

enrichment_score <- data.frame(pvalue=c(Normal_TM_pvalue,Tumor_TM_pvalue,Normal_TMPandCd_pvalue,Tumor_TMPandCd_pvalue),ratio=c(Normal_TM_ratio,Tumor_TM_ratio,Normal_TMPandCd_ratio,Tumor_TMPandCd_ratio))
enrichment_score$human_mouse_paired <- c("Normal_TM","Normal_TMPandCd","Tumor_TM","Tumor_TMPandCd")
enrichment_score$stage <- c("Normal","Normal","Tumor","Tumor")
enrichment_score$mouse_sig <- c("TM","TMPandCd","TM","TMPandCd")

enrichment_score$stage <- factor(enrichment_score$stage,levels=c("Normal","Tumor"))
rownames(enrichment_score) <- enrichment_score$human_mouse_paired
ff <- ggplot(enrichment_score, aes(mouse_sig,stage)) + geom_point() +
geom_point(aes(size=-log10(pvalue),color=ratio))+
scale_colour_gradient(low = "#9ecae1",high = "#cb181d")+
labs(color="ratio",size="-log10(pvalue)",y="stage") + theme_classic()
ggsave(ff,file="enrichment_score.png",width=3,height=3.5,dpi=1080)
```



![](Bulk_RNAseq.assets/enrichment_score.png)

### Part 6.  The Kaplan-Meier survival curves of patients with high expression levels of TMP and TMPC signature genes in the TCGA-STAD cohort. 

```R
spcific_genes_in_9_subcotous <- read.csv("./all_mouse_spcific_genes_in_9_subcotous_to_human.csv")
TMP_sig <- subset(spcific_genes_in_9_subcotous,TMP_pvalue <=0.05 & TMP_log2FoldChange>2.725)
TMPCd_sig <- subset(spcific_genes_in_9_subcotous,TMPCd_pvalue <=0.05 & TMPCd_log2FoldChange>2.725)
TCGA_STAD <- fread("./TCGA-STAD_RNA_expr.csv",sep=",")
TCGA_STAD <- as.data.frame(TCGA_STAD)
TCGA_STAD <- na.omit(TCGA_STAD)
rownames(TCGA_STAD) <- TCGA_STAD$V1
TCGA_STAD <- TCGA_STAD[,-1]
library(org.Hs.eg.db)
library(AnnotationDbi)
TCGA_STAD$symbol <- mapIds(x = org.Hs.eg.db,
            keys = rownames(TCGA_STAD),
            keytype ="ENSEMBL",
            column ="SYMBOL",
            multiVals="first")
TCGA_STAD <- na.omit(TCGA_STAD)
TCGA_STAD <- TCGA_STAD[!duplicated(TCGA_STAD$symbol),]
rownames(TCGA_STAD) <- TCGA_STAD$symbol
TCGA_STAD <- TCGA_STAD[,-ncol(TCGA_STAD)]
colnames(TCGA_STAD) <- substring(colnames(TCGA_STAD),1,15)
TCGA_STAD <- as.data.frame(t(TCGA_STAD))
TCGA_STAD1 <- TCGA_STAD
TCGA_STAD1 <- log(TCGA_STAD+1,2)
TCGA_STAD1 <- as.data.frame(t(TCGA_STAD1))
TCGA_STAD1 <- TCGA_STAD1[apply(TCGA_STAD1,1,sd)!=0,]
TCGA_STAD_clinical_sample <- read.csv("./data_bcr_clinical_data_sample.txt",sep="\t")
TCGA_STAD_clinical_sample <- TCGA_STAD_clinical_sample[-c(1:4),]
head(TCGA_STAD_clinical_sample)
rownames(TCGA_STAD_clinical_sample) <- TCGA_STAD_clinical_sample$X.Sample.Identifier
TCGA_STAD_clinical_patient <- read.csv("./data_bcr_clinical_data_patient.txt",sep="\t")
TCGA_STAD_clinical_patient <- TCGA_STAD_clinical_patient[-c(1:4),]
head(TCGA_STAD_clinical_patient)
rownames(TCGA_STAD_clinical_patient) <- TCGA_STAD_clinical_patient$X.Patient.Identifier
length(rownames(TCGA_STAD_clinical_sample))
length(rownames(TCGA_STAD_clinical_patient))
length(intersect(rownames(TCGA_STAD_clinical_sample),rownames(TCGA_STAD_clinical_patient)))
TCGA_STAD_clinical1 <- as.data.frame(cbind(TCGA_STAD_clinical_sample[intersect(rownames(TCGA_STAD_clinical_sample),rownames(TCGA_STAD_clinical_patient)),],TCGA_STAD_clinical_patient[intersect(rownames(TCGA_STAD_clinical_sample),rownames(TCGA_STAD_clinical_patient)),]))
TCGA_STAD_clinical1 <-TCGA_STAD_clinical1[,c("X.Patient.Identifier","Sample.Identifier","Sample.Type","Cancer.Type","Neoplasm.Histologic.Grade","Sex","Lymph.Node.s..Examined.Number","American.Joint.Committee.on.Cancer.Tumor.Stage.Code","Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code","American.Joint.Committee.on.Cancer.Metastasis.Stage.Code","Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code","Overall.Survival..Months.","Overall.Survival.Status")]
colnames(TCGA_STAD_clinical1) <- c("Patient","Sample","Sample.Type","Cancer.Type","Neoplasm.Histologic.Grade","Sex","Lymph.Node.s..Examined.Number","Cancer.Tumor.Stage.Code","Neoplasm.Disease.Lymph.Node.Stage.Cancer.Code","Cancer.Metastasis.Stage.Code","Neoplasm.Disease.Stage.Cancer.Code","Overall.Survival..Months.","Overall.Survival.Status")
rownames(TCGA_STAD_clinical1) <- TCGA_STAD_clinical1$Sample

TCGA_STAD_clinical <- read.csv("./TCGA-STAD_clinical.csv")
rownames(TCGA_STAD_clinical) <- TCGA_STAD_clinical$bcr_patient_barcode
rownames(TCGA_STAD_clinical1) <- TCGA_STAD_clinical1$Patient
both_id <- intersect(rownames(TCGA_STAD_clinical1),rownames(TCGA_STAD_clinical))
TCGA_STAD_clinical <- as.data.frame(cbind(TCGA_STAD_clinical1[both_id,],TCGA_STAD_clinical[both_id,]))
TCGA_STAD_clinical <- TCGA_STAD_clinical[,c("Patient","Sample","Sample.Type","Cancer.Type","Neoplasm.Histologic.Grade","Sex","Lymph.Node.s..Examined.Number","Cancer.Tumor.Stage.Code","Neoplasm.Disease.Lymph.Node.Stage.Cancer.Code","Cancer.Metastasis.Stage.Code","Neoplasm.Disease.Stage.Cancer.Code","Overall.Survival..Months.","Overall.Survival.Status","days_to_death","vital_status","ajcc_pathologic_stage","days_to_last_follow_up")]
rownames(TCGA_STAD_clinical) <- TCGA_STAD_clinical$Sample
both_id <- intersect(colnames(TCGA_STAD1),rownames(TCGA_STAD_clinical))
TCGA_STAD_clinical <- TCGA_STAD_clinical[both_id,]
TCGA_STAD1 <- TCGA_STAD1[,both_id]
TCGA_STAD_t <- as.data.frame(t(TCGA_STAD1))

TMP_sel_genes <- as.character(TMP_sig$HGNC.symbol)
TMPCd_sel_genes <- as.character(TMPCd_sig$HGNC.symbol)
TCGA_STAD_tmp <- data.frame(TMP_sig=as.character(apply(TCGA_STAD_t[,intersect(colnames(TCGA_STAD_t),TMP_sel_genes)],1,mean)),TMPCd_sig=as.character(apply(TCGA_STAD_t[,intersect(colnames(TCGA_STAD_t),TMPCd_sel_genes)],1,mean)),row.names=rownames(TCGA_STAD_t))
TCGA_STAD_tmp$TMP_sig <- as.numeric(as.character(scale((as.numeric(as.character(TCGA_STAD_tmp$TMP_sig))))))
TCGA_STAD_tmp$TMPCd_sig <- as.numeric(as.character(scale((as.numeric(as.character(TCGA_STAD_tmp$TMPCd))))))
TCGA_STAD_tmp$group <- apply(TCGA_STAD_tmp, 1, function(t) colnames(TCGA_STAD_tmp)[which.max(t)])
TCGA_STAD_Sel <- TCGA_STAD_tmp
TCGA_STAD_clinical_sel <- TCGA_STAD_clinical[rownames(TCGA_STAD_Sel),]
All_merge_clinical <- cbind(TCGA_STAD_clinical_sel,TCGA_STAD_Sel)
meta <- All_merge_clinical
meta[is.na(meta$days_to_last_follow_up),]$days_to_last_follow_up <- "HHH"
tmp <- subset(meta,days_to_last_follow_up=="HHH")
tmp$days_to_last_follow_up <- tmp$days_to_death
no_na <- meta[setdiff(rownames(meta),rownames(tmp)),]
all_merge <- rbind(tmp,no_na)
all_merge <- subset(all_merge,days_to_last_follow_up != "HHH")
all_merge$vital_status <- as.character(all_merge$vital_status)
all_merge$status <- ifelse(all_merge$vital_status=="Alive",0,1)
all_merge$days_to_last_follow_up <- as.numeric(all_merge$days_to_last_follow_up)

library("survival")
library("survminer")
meta <- all_merge
meta <- meta[!duplicated(meta$Patient),]
meta[is.na(meta$days_to_last_follow_up),]$days_to_last_follow_up <- "HHH"
tmp <- subset(meta,days_to_last_follow_up=="HHH")
tmp$days_to_last_follow_up <- tmp$days_to_death
no_na <- meta[setdiff(rownames(meta),rownames(tmp)),]
all_merge <- rbind(tmp,no_na)
all_merge <- subset(all_merge,days_to_last_follow_up != "HHH")
all_merge$vital_status <- as.character(all_merge$vital_status)
all_merge$status <- ifelse(all_merge$vital_status=="Alive",0,1)
all_merge$days_to_last_follow_up <- as.numeric(all_merge$days_to_last_follow_up)

fit <- survfit(Surv(days_to_last_follow_up, status) ~ group, data = all_merge)
ff <- ggsurvplot(fit, data = all_merge,surv.median.line = "hv",xlim= c(0,1000), break.x.by=200,
pval = TRUE,
ggtheme = theme_bw(),
risk.table=TRUE)
ggsave(ff$plot,file="FigS3_TMP_TMPC_sig_surv.png",dpi=1080)
```

![](assets/FigS3_TMP_TMPC_sig_surv-16566009814952.png)
