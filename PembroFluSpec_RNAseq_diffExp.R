#+ fig.width=8, fig.height=8
#+ warnings = FALSE
#' ##--------- call libraries --

library("genefilter")
library("ggplot2")
library("grid");  library("reshape2")
library("ggrepel")
library("RColorBrewer")
library("DESeq2")
library("BiocParallel")
library("WGCNA")
library("Rtsne")
library("pheatmap")
library("tools")
library("viridis")
library("scales")
library("gridExtra")
library("GSVA")
library("GSEABase")
library("stringr")
library("dplyr")
library("rstatix")
library("tidyr")
sessionInfo()
source('D:/Pembro-Fluvac/Analysis/Ranalysis/PembroFluSpec_PlottingFunctions.R')

loadData <- function()    #  this function's actions are being done in the Ranalysis.R file, line 92
{
  return(readRDS(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergedData.RDS"))
}
if (! exists("mergedData"))   {   mergedData <- loadData()  }      # run if mergedData isn't already loaded into the global environment
equivTable <- mergedData[!duplicated(mergedData$Subject), c("Subject","dbCode")]
equivTable$Subject2 <- equivTable$Subject
equivTable$Subject2[grep(pattern="19616", x = equivTable$Subject2)] <- str_c("X",equivTable$Subject2[grep(pattern="19616",x=equivTable$Subject2)])
equivTable$Subject2 <- str_replace(equivTable$Subject2, pattern="-", replacement=".")

## ***     read in data and make metaData  ************** 
dataMatrix <- read.table("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/PORTresults/NORMALIZED_DATA/SPREADSHEETS/FINAL_master_list_of_gene_counts_MIN.Fastqs.txt", sep="", 
                         header=T, stringsAsFactors = F)
rownames(dataMatrix) <- make.names(dataMatrix$geneSymbol, unique=TRUE)
dataMatrix$id <- dataMatrix$geneCoordinate <- dataMatrix$geneSymbol <- NULL 

makeMetaData <- function( colnVector)
{ metaData <- data.frame(row.names= colnVector);  metaData$condition <- "empty"
  metaData$Cohort <- metaData$Subject <- metaData$TimeCategory <- metaData$Subset <- "a"
  metaData$Cohort[grep("196",rownames(metaData),value=F)] <- "anti-PD-1" ;  metaData$Cohort[grep("FS1819",rownames(metaData),value=F)] <- "Healthy"
  metaData$TimeCategory[grep("oW",rownames(metaData), value=F)] <- "oneWeek"; metaData$TimeCategory[grep("bL",rownames(metaData), value=F)] <- "baseline"; 
  metaData$Subset[grep("HiHi",rownames(metaData),value=F)] <- "HiHi_cTfh"; metaData$Subset[grep("ABC",rownames(metaData),value=F)] <- "ABC";
  metaData$Subset[grep("_PB_",rownames(metaData),value=F)] <- "PB"; metaData$Subset[grep("naiveB",rownames(metaData),value=F)] <- "navB";
  metaData$Subject <- substr(rownames(metaData), start = 1, stop = 10)
  metaData$condition <- paste0( metaData$Cohort,"_",metaData$Subset,"_",metaData$TimeCategory)
  return(metaData)      }

metaData <- makeMetaData(colnVector = colnames(dataMatrix))

metaData.dbCode <- merge(x = metaData, y=equivTable, all.x=T, by.x = "Subject", by.y="Subject2")
rownames(metaData.dbCode) <- rownames(metaData)

#' ## ***     read in log-transformed data   **************

logDataMatrix <- read.csv( file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/logDataMatrix.csv", stringsAsFactors = F, header = 1, row.names = 1)  # created in QC script


#' ##  global heatmap of relevant/interesting genes  --------
probeGenes <- logDataMatrix[ c("CXCR5", "PRDM1",  "IFNG","MAF","CCR6","CXCR3","GATA3",
               "BTLA","TNFRSF4", "CD38","TIGIT","POU2AF1","MKI67","SH2D1A", "TOX2", "BIRC5", "TBK1", "CD19","MS4A1","FUT8","ST6GAL1","B4GALT1",
               "JCHAIN", "CMTM6", "AK9", "ZNF770","CD3E","CD4", "PAX5", "SDC1") , ]  # some genes taken from genomicscape.com

metaData <- makeMetaData(colnVector = colnames(dataMatrix))
annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "anti-PD-1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression", show_colnames = F, cutree_cols = 5, cutree_rows = 4
         #, filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/SelectedGenes_allSubsets_allTimePoints_allSubjects_heatmap.pdf"
         )

## ------------- Differential expression by subgroup and age category ----------

library("BiocParallel")
register(SnowParam(10))
metaData <- makeMetaData(colnVector = colnames(dataMatrix))
fullDataset <- DESeqDataSetFromMatrix(countData=dataMatrix, colData=metaData, design= ~ condition)
dds <- estimateSizeFactors(fullDataset);   idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) > 20  # hard filter for genes with at least 20 counts in 25% of samples
fullDataset <- fullDataset[idx,]
# DESdata_fullDataset <- DESeq(fullDataset, parallel=TRUE)


# ----------- PB vs naiveB at baseline   

meta <- metaData[ grep(paste0(c("PB_bL", "naiveB_bL"), collapse="|"), rownames(metaData)),1:4]
meta$condition <- factor(meta$condition);  data <- dataMatrix[ , grep(paste0(c("PB_bL", "naiveB_bL"), collapse="|"), colnames(dataMatrix) ) ]
PB_v_nav_bL <- DESeqDataSetFromMatrix(countData=data, colData=meta, design = ~Subset + Subject)
PB_v_nav_bL <- PB_v_nav_bL[idx,]   # keep hard filter from above for consistency
DESdata_PBvNav_AH_bL <- DESeq(PB_v_nav_bL, parallel=TRUE)

PBvNav_AH_bL <- as.data.frame(results(DESdata_PBvNav_AH_bL, contrast = c("Subset", "PB", "navB") ))  # pos stats = first elem in comparison
#write.csv(PBvNav_AH_bL, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/PBvNav_AH_bL.csv")
volcanoPlot(PBvNav_AH_bL, repelThresh = 1e-50, title = "PB vs navB", leftLabel = "NaiveB", rightLabel = "PB" )

# ----------- ABC vs naiveB at baseline   

meta <- metaData[ grep(paste0(c("ABC_bL", "naiveB_bL"), collapse="|"), rownames(metaData)),1:4]
meta$condition <- factor(meta$condition);  data <- dataMatrix[ , grep(paste0(c("ABC_bL", "naiveB_bL"), collapse="|"), colnames(dataMatrix) ) ]
ABC_v_nav_bL <- DESeqDataSetFromMatrix(countData=data, colData=meta, design = ~Subset + Subject)
ABC_v_nav_bL <- ABC_v_nav_bL[idx,]   # keep hard filter from above for consistency
DESdata_ABCvNav_AH_bL <- DESeq(ABC_v_nav_bL, parallel=TRUE)

ABCvNav_AH_bL <- as.data.frame(results(DESdata_ABCvNav_AH_bL, contrast = c("Subset", "ABC", "navB") ))  # pos stats = first elem in comparison
#write.csv(ABCvNav_AH_bL, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/ABCvNav_AH_bL.csv")
volcanoPlot(ABCvNav_AH_bL, repelThresh = 1e-30, title = "ABC vs navB", leftLabel = "NaiveB", rightLabel = "ABC" )



# ----------- HiHi in aPD1 vs HC   baseline  

HiHi_aPD1_v_HC_bLmeta <- metaData[ grep("HiHi_bL", rownames(metaData)),1:2]
HiHi_aPD1_v_HC_bLmeta$condition <- factor(HiHi_aPD1_v_HC_bLmeta$condition);  HiHi_aPD1_v_HC_bLdata <- dataMatrix[ , grep("HiHi_bL", colnames(dataMatrix) ) ]
HiHi_aPD1_v_HC_bL <- DESeqDataSetFromMatrix(countData=HiHi_aPD1_v_HC_bLdata, colData=HiHi_aPD1_v_HC_bLmeta, design = ~condition)
HiHi_aPD1_v_HC_bL <- HiHi_aPD1_v_HC_bL[idx,]   # keep hard filter from above for consistency
DESdata_HiHi_AvH_bL <- DESeq(HiHi_aPD1_v_HC_bL, parallel=TRUE)

HiHi_AvH_bL <- as.data.frame(results(DESdata_HiHi_AvH_bL, contrast = c("condition", "anti-PD-1_HiHi_cTfh_baseline", "Healthy_HiHi_cTfh_baseline") ))  # pos stats = first elem in comparison
# write.csv(HiHi_AvH_bL, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/HiHi_AvH_bL.csv")
volcanoPlot(HiHi_AvH_bL, repelThresh = 0.10, title = "HiHi at baseline in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "anti-PD-1" )
HiHi_AvH_bL <- HiHi_AvH_bL[ order(HiHi_AvH_bL$stat, decreasing = F), ]; HiHi_AvH_bL <- HiHi_AvH_bL[ which(!is.na(HiHi_AvH_bL$stat)), ]

diffExpGenes <- row.names( rbind(head(HiHi_AvH_bL, 30), tail(HiHi_AvH_bL, 30)))
diffExpGenes <- diffExpGenes[-grep("DEFA3",diffExpGenes, value=F)]            # exclude because padj = NA
probeGenes <- logDataMatrix[diffExpGenes, grep("HiHi_bL",colnames(logDataMatrix),value=F)]
labels <- rownames(probeGenes); # labels[!labels %in% displayGenes] <- "" 
newnames <- lapply(labels, function(x) bquote(italic(.(x))))
annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]; 
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "anti-PD-1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  ) # 
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors, show_colnames = F, 
         labels_row = as.expression(newnames), 
         fontsize_row = 7, color=inferno(100), main = "ICOS+CD38+ cTfh at baseline", border_color = NA
        # , filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HvA_bL_heatmap_AllLabels.pdf", device="pdf"
)
dev.off()
temp <- probeGenes
names(temp) <- metaData.dbCode[names(temp),"dbCode"]
# write.csv(temp, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e5i-1.csv")


# ----------- HiHi in aPD1 vs HC   oneWeek 

HiHi_aPD1_v_HC_oWmeta <- metaData[ grep("HiHi_oW", rownames(metaData)),1:2]
HiHi_aPD1_v_HC_oWmeta$condition <- factor(HiHi_aPD1_v_HC_oWmeta$condition);  HiHi_aPD1_v_HC_oWdata <- dataMatrix[ , grep("HiHi_oW", colnames(dataMatrix) ) ]
HiHi_aPD1_v_HC_oW <- DESeqDataSetFromMatrix(countData=HiHi_aPD1_v_HC_oWdata, colData=HiHi_aPD1_v_HC_oWmeta, design = ~condition)
HiHi_aPD1_v_HC_oW <- HiHi_aPD1_v_HC_oW[idx,]   # keep hard filter from above for consistency
DESdata_HiHi_AvH_oW <- DESeq(HiHi_aPD1_v_HC_oW, parallel=TRUE)

HiHi_AvH_oW <- as.data.frame(results(DESdata_HiHi_AvH_oW, contrast = c("condition", "anti-PD-1_HiHi_cTfh_oneWeek", "Healthy_HiHi_cTfh_oneWeek") ))  # pos stats = first elem in comparison
volcanoPlot(HiHi_AvH_oW, repelThresh = 0.15, title = "HiHi at oneWeek in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/HiHi_H-v-aPD1_d7_volcano.pdf", device="pdf")
HiHi_AvH_oW <- HiHi_AvH_oW[ order(HiHi_AvH_oW$stat, decreasing = F), ];  HiHi_AvH_oW <- HiHi_AvH_oW[ which(!is.na(HiHi_AvH_oW$stat)), ]
# write.csv(HiHi_AvH_oW, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/HiHi_AvH_oW.csv")

diffExpGenes <- row.names( rbind(head(HiHi_AvH_oW, 50), tail(HiHi_AvH_oW, 50)))
probeGenes <- logDataMatrix[diffExpGenes, grep("HiHi_oW",colnames(logDataMatrix),value=F)]

displayGenes <- c("MKI67","HELLS","ESPL1","ARID1A","TNFRSF1B","DDX11", "CCDC84", "ARID1A", "DDX10", "ST6GAL1", "IFI44", "OASL", "IKZF1", 
                  "PTGES3","IFI44L", "DDX60","FOXO1"); 
labels <- rownames(probeGenes); labels[!labels %in% displayGenes] <- "" 
newnames <- lapply(labels, function(x) bquote(italic(.(x))))
annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]; 
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "anti-PD-1" = "#FFB18C"), Time = c("baseline"="grey90","oneWeek"="grey40")  ) 
annotation$Subset <- NULL;  names(annotation)[1] <- "Time"
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors, show_colnames = F, 
         labels_row = as.expression(newnames),
         fontsize_row = 13, color=inferno(100), main = "ICOS+CD38+ cTfh at one week", border_color = NA, 
        # , filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HvA_oW_heatmap_SelectLabels.pdf", device="pdf"         # ******  further modified in Illustrator
)
dev.off()
temp <- probeGenes
names(temp) <- metaData.dbCode[names(temp),"dbCode"]
# write.csv(temp, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3a.csv")


# ----------- PB in aPD1 vs HC   baseline  

PB_aPD1_v_HC_bLmeta <- metaData[ grep("PB_bL", rownames(metaData)),1:2]
PB_aPD1_v_HC_bLmeta$condition <- factor(PB_aPD1_v_HC_bLmeta$condition);  PB_aPD1_v_HC_bLdata <- dataMatrix[ , grep("PB_bL", colnames(dataMatrix) ) ]
PB_aPD1_v_HC_bL <- DESeqDataSetFromMatrix(countData=PB_aPD1_v_HC_bLdata, colData=PB_aPD1_v_HC_bLmeta, design = ~condition)
PB_aPD1_v_HC_bL <- PB_aPD1_v_HC_bL[idx,]   # keep hard filter from above for consistency
DESdata_PB_AvH_bL <- DESeq(PB_aPD1_v_HC_bL, parallel=TRUE)

PB_AvH_bL <- as.data.frame(results(DESdata_PB_AvH_bL, contrast = c("condition", "anti-PD-1_PB_baseline", "Healthy_PB_baseline") ))  # pos stats = first elem in comparison

volcanoPlot(PB_AvH_bL, repelThresh = 0.15, title = "PB at baseline in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )
PB_AvH_bL <- PB_AvH_bL[ order(PB_AvH_bL$stat, decreasing = F), ]; PB_AvH_bL <- PB_AvH_bL[ which(!is.na(PB_AvH_bL$stat)), ]
# write.csv(PB_AvH_bL, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/PB_AvH_bL.csv")

diffExpGenes <- row.names( rbind(head(PB_AvH_bL, 30), tail(PB_AvH_bL, 30)))
diffExpGenes <- diffExpGenes[-grep("IGHD",diffExpGenes, value=F)]     # exclude because padj = NA
probeGenes <- logDataMatrix[diffExpGenes, grep("PB_bL",colnames(logDataMatrix),value=F)]
labels <- rownames(probeGenes); # labels[!labels %in% displayGenes] <- "" 
newnames <- lapply(labels, function(x) bquote(italic(.(x))))
annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "anti-PD-1" = "#FFB18C"),Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  ) #  
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors, show_colnames = F, 
         labels_row = as.expression(newnames), 
         fontsize_row = 7, color=inferno(100), main = "ASC at baseline", border_color = NA
        # , filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/PB_HvA_bL_heatmap_AllLabels.pdf", device="pdf"
)
dev.off()
temp <- probeGenes
names(temp) <- metaData.dbCode[names(temp),"dbCode"]
# write.csv(temp, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e5i-2.csv")

# ----------- PB in aPD1 vs HC   oneWeek 

PB_aPD1_v_HC_oWmeta <- metaData[ grep("PB_oW", rownames(metaData)),1:2]
PB_aPD1_v_HC_oWmeta$condition <- factor(PB_aPD1_v_HC_oWmeta$condition);  PB_aPD1_v_HC_oWdata <- dataMatrix[ , grep("PB_oW", colnames(dataMatrix) ) ]
PB_aPD1_v_HC_oW <- DESeqDataSetFromMatrix(countData=PB_aPD1_v_HC_oWdata, colData=PB_aPD1_v_HC_oWmeta, design = ~condition)
PB_aPD1_v_HC_oW <- PB_aPD1_v_HC_oW[idx,]   # keep hard filter from above for consistency
DESdata_PB_AvH_oW <- DESeq(PB_aPD1_v_HC_oW, parallel=TRUE)

PB_AvH_oW <- as.data.frame(results(DESdata_PB_AvH_oW, contrast = c("condition", "anti-PD-1_PB_oneWeek", "Healthy_PB_oneWeek") ))  # pos stats = first elem in comparison
volcanoPlot(PB_AvH_oW, repelThresh = 0.15, title = "PB at oneWeek in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )
PB_AvH_oW <- PB_AvH_oW[ order(PB_AvH_oW$stat, decreasing = F), ]; PB_AvH_oW <- PB_AvH_oW[ which(!is.na(PB_AvH_oW$stat)), ]
# write.csv(PB_AvH_oW, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/PB_AvH_oW.csv")

diffExpGenes <- row.names( rbind(head(PB_AvH_oW, 50), tail(PB_AvH_oW, 50)))
probeGenes <- logDataMatrix[diffExpGenes, grep("PB_oW",colnames(logDataMatrix),value=F)]

displayGenes <- c("BUB1","ESPL1","TNIP2","TRAIP","AURKB", "KCNA3", "DUSP1", "MAFF", "IFNGR", "E2F2", "SMARCD1", "IL18RAP", "CCL2", "ADAM28","IGHV3.15","IGKV1.5"); 
labels <- rownames(probeGenes); labels[!labels %in% displayGenes] <- "" 
newnames <- lapply(labels, function(x) bquote(italic(.(x))))
annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]; annotation$Subset <- NULL;  names(annotation)[1] <- "Time"
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "anti-PD-1" = "#FFB18C"),Time = c("baseline"="grey90","oneWeek"="grey40")  ) #  Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), 
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors, show_colnames = F, 
         labels_row = as.expression(newnames), 
         fontsize_row = 13, color=inferno(100), main = "ASC at one week", border_color = NA
         # , filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/PB_HvA_oW_heatmap_SelectLabels.pdf", device="pdf"      # ******  further modified in Illustrator
)
dev.off()
temp <- probeGenes
names(temp) <- metaData.dbCode[names(temp),"dbCode"]
# write.csv(temp, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3c.csv")


# ----------- ABC in aPD1 vs HC   oneWeek 

ABC_aPD1_v_HC_oWmeta <- metaData[ grep("ABC_oW", rownames(metaData)),1:2]
ABC_aPD1_v_HC_oWmeta$condition <- factor(ABC_aPD1_v_HC_oWmeta$condition);  ABC_aPD1_v_HC_oWdata <- dataMatrix[ , grep("ABC_oW", colnames(dataMatrix) ) ]
ABC_aPD1_v_HC_oW <- DESeqDataSetFromMatrix(countData=ABC_aPD1_v_HC_oWdata, colData=ABC_aPD1_v_HC_oWmeta, design = ~condition)
ABC_aPD1_v_HC_oW <- ABC_aPD1_v_HC_oW[idx,]   # keep hard filter from above for consistency
DESdata_ABC_AvH_oW <- DESeq(ABC_aPD1_v_HC_oW, parallel=TRUE)

ABC_AvH_oW <- as.data.frame(results(DESdata_ABC_AvH_oW, contrast = c("condition", "anti-PD-1_ABC_oneWeek", "Healthy_ABC_oneWeek") ))  # pos stats = first elem in comparison
volcanoPlot(ABC_AvH_oW, repelThresh = 0.05, title = "ABC at oneWeek in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )

ABC_AvH_oW <- ABC_AvH_oW[ order(ABC_AvH_oW$stat, decreasing = F), ]; ABC_AvH_oW <- ABC_AvH_oW[ which(!is.na(ABC_AvH_oW$stat)), ]
# write.csv(ABC_AvH_oW, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/ABC_AvH_oW.csv")

diffExpGenes <- row.names( rbind(head(ABC_AvH_oW, 50), tail(ABC_AvH_oW, 50)))
probeGenes <- logDataMatrix[diffExpGenes, grep("ABC_oW",colnames(logDataMatrix),value=F)]
annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "anti-PD-1" = "#FFB18C"),Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  ) #  

pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors, show_colnames = F, # labels_row = labels, 
         fontsize_row = 4, color=inferno(100), main = "ABC at one week", border_color = NA
#         , filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/ABC_HvA_oW_heatmap_AllLabels.pdf", device="pdf"
)
dev.off()



# ************ Ellebedy subsets validation  - ASC-vs-naive ********************
pathwayPos <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/PBvNav_AH_bL_externalGenesets.GseaPreranked.1610193172180/gsea_report_for_na_pos_1610193172180.tsv", sep="\t")
pathwayNeg <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/PBvNav_AH_bL_externalGenesets.GseaPreranked.1610193172180/gsea_report_for_na_neg_1610193172180.tsv", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/PBvNav_AH_bL_externalGenesets.GseaPreranked.1610193172180/GSE68245_ASC-VS-NAVB.tsv", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "ASC-VS-NAVB", title = "GSE68245: ASC", leftLabel = "ASC", rightLabel = "Naive B", cohortCompare =F, 
                    legendUpperRight = T)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "ASC-VS-NAVB", title = "GSE68245: ASC", leftLabel = " ", rightLabel = " ", cohortCompare =F, 
                    legendUpperRight = T) + scale_x_continuous(limits=c(0,10300), breaks=seq(0,10000,2500))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/ASCvNav_AH_EllebedyReference_bL.pdf", device="pdf", width=7, height=6)
# write.csv(singlePathway, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e5d.csv")

# ************ Ellebedy subsets validation  - ABC-vs-naive ********************
pathwayPos <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/ABCvNav_AH_bL_externalGenesets.GseaPreranked.1610193286316/gsea_report_for_na_pos_1610193286316.tsv", sep="\t")
pathwayNeg <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/ABCvNav_AH_bL_externalGenesets.GseaPreranked.1610193286316/gsea_report_for_na_neg_1610193286316.tsv", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/ABCvNav_AH_bL_externalGenesets.GseaPreranked.1610193286316/GSE68245_ABC-VS-NAVB.tsv", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "ABC-VS-NAVB", title = "GSE68245: ABC", leftLabel = "ABC", rightLabel = "Naive B", cohortCompare =F, 
                    legendUpperRight = T)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "ABC-VS-NAVB", title = "GSE68245: ABC", leftLabel = " ", rightLabel = " ", cohortCompare =F, 
                    legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/ABCvNav_AH_EllebedyReference_bL.pdf", device="pdf", width=7, height=6)  
# write.csv(singlePathway, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e5e.csv")



# ------------- Gene Ontology ----------------

##' cTfh responses
GO.aPD1 <- read.csv(file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GeneOntology/HiHi_AvH_oW_enrichApd1/HiHi_AvH_oW_enrichApd1_GeneOntologies.csv")
GO.Healthy <- read.csv(file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GeneOntology/HiHi_AvH_oW_enrichHealthy/HiHi_AvH_oW_enrichHealthy_GeneOntologies.csv")
GO.aPD1$Symbols <- GO.Healthy$Symbols <- GO.aPD1$Genes <- GO.Healthy$Genes <- NULL
GO.aPD1$GroupID <- "aPD1"; GO.Healthy$GroupID <- "Healthy"
GO.aPD1$LogP <- -GO.aPD1$LogP
# GO.aPD1[nrow(GO.aPD1)+1,] <- NA             # create an empty row as a visual aid
mergedOntology <- rbind(GO.aPD1[ c(1:4,16),], GO.Healthy[ c(2,3,10,13), ] ) 

mergeResults <- mergedOntology;  title = "ICOS+CD38+ cTfh \nat one week"; leftLabel = "Healthy"; rightLabel = "aPD1" 
mergeResults$Description[grep("Signaling by Interleukins",mergeResults$Description)] <- "Signaling by\nInterleukins"
mergeResults$Description[grep("positive regulation of cytokine production",mergeResults$Description)] <- "positive regulation\nof cytokine production"
mergeResults$Description[grep("negative regulation of cell cycle",mergeResults$Description)] <- "negative regulation\nof cell cycle"
mergeResults$Description[grep("microtubule cytoskeleton organization",mergeResults$Description)] <- "microtubule\ncytoskeleton\norganization"
mergeResults$Description <- factor(mergeResults$Description, levels = mergeResults$Description[order(mergeResults$LogP, decreasing = F)])
left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.03,hjust=0, gp=gpar(col="black", fontsize=16)))
right_grob <- grobTree(textGrob(rightLabel, x=0.8,y=0.03,hjust=0, gp=gpar(col="black", fontsize=16)))
ggplot(data=mergeResults) + geom_point(aes(x=Description, y=LogP), size=6) + 
  geom_bar( data = subset(mergeResults, `GroupID` == "aPD1"), aes(x=Description, y=LogP) , stat="Identity", width=0.01, color="#FFB18C", size=1) +
  geom_bar( data = subset(mergeResults, `GroupID` == "Healthy"), aes(x=Description, y=LogP) , stat="Identity", width=0.01, color="#7FAEDB", size=1) +
  coord_flip() + theme_bw() + ggtitle(title) + ylab("-Log (P)") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=18), axis.text = element_text(size=18), title = element_text(size=18)) + 
  annotation_custom(left_grob) + annotation_custom(right_grob) + scale_y_continuous(labels = abs)
ggplot(data=mergeResults) + geom_point(aes(x=Description, y=LogP), size=7) + 
  geom_bar( data = subset(mergeResults, `GroupID` == "aPD1"), aes(x=Description, y=LogP) , stat="Identity", width=0.03, color="#FFB18C", size=1) +
  geom_bar( data = subset(mergeResults, `GroupID` == "Healthy"), aes(x=Description, y=LogP) , stat="Identity", width=0.03, color="#7FAEDB", size=1) +
  coord_flip() + theme_bw() + ggtitle(title) + ylab("-Log (P)") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=24), axis.text = element_text(size=24, color="black"), title = element_text(size=24), axis.text.x = element_text(size=24)) + 
  #annotation_custom(left_grob) + annotation_custom(right_grob) + 
  scale_y_continuous(labels = abs)  + theme(panel.grid.major=element_blank()) + geom_hline(yintercept=0, size=0.5)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GeneOntology_cTfh_AvH_oW.pdf", width=6, height=11)
# write.csv(mergeResults, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3b.csv")


##'  Plasmablasts
GO.aPD1 <- read.csv(file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GeneOntology/PB_AvH_oW_enrichApd1/PB_AvH_oW_enrichApd1_GeneOntologies.csv")
GO.Healthy <- read.csv(file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GeneOntology/PB_AvH_oW_enrichHealthy/PB_AvH_oW_enrichHealthy_GeneOntologies.csv")
GO.aPD1$Symbols <- GO.Healthy$Symbols <- GO.aPD1$Genes <- GO.Healthy$Genes <- NULL
GO.aPD1$GroupID <- "aPD1"; GO.Healthy$GroupID <- "Healthy"
GO.aPD1$LogP <- -GO.aPD1$LogP
# GO.aPD1[nrow(GO.aPD1)+1,] <- NA             # create an empty row as a visual aid
mergedOntology <- rbind(GO.aPD1[ c(1:3),], GO.Healthy[ c(1:3), ] ) 
mergedOntology$Description[4] <- substr(mergedOntology$Description[4],1,38)
mergedOntology$Description[4] <- substr(mergedOntology$Description[4],1,38)

mergeResults <- mergedOntology;  title = "ASC at one week"; leftLabel = "Healthy"; rightLabel = "aPD1" 

mergeResults$Description[grep("PID FOXM1 PATHWAY",mergeResults$Description)] <- "PID FOXM1\nPATHWAY"
mergeResults$Description[grep("TNF-alpha/Nf-kappa B signaling complex",mergeResults$Description)] <- "TNFa/Nf-kB\nsignaling complex"
mergeResults$Description[grep("regulation of chromosome segregation",mergeResults$Description)] <- "regulation of\nchromosome\nsegregation"
mergeResults$Description[grep("mitochondrion organization",mergeResults$Description)] <- "mitochondrion\norganization"

mergeResults$Description <- factor(mergeResults$Description, levels = mergeResults$Description[order(mergeResults$LogP, decreasing = F)])
left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.03,hjust=0, gp=gpar(col="black", fontsize=16)))
right_grob <- grobTree(textGrob(rightLabel, x=0.8,y=0.03,hjust=0, gp=gpar(col="black", fontsize=16)))
ggplot(data=mergeResults) + geom_point(aes(x=Description, y=LogP), size=6) + 
  geom_bar( data = subset(mergeResults, `GroupID` == "aPD1"), aes(x=Description, y=LogP) , stat="Identity", width=0.01, color="#FFB18C", size=1) +
  geom_bar( data = subset(mergeResults, `GroupID` == "Healthy"), aes(x=Description, y=LogP) , stat="Identity", width=0.01, color="#7FAEDB", size=1) +
  coord_flip() + theme_bw() + ggtitle(title) + ylab("-Log (P)") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=18), axis.text = element_text(size=18), title = element_text(size=18)) + 
  annotation_custom(left_grob) + annotation_custom(right_grob) + scale_y_continuous(labels = abs, limits = c(-7,12))
ggplot(data=mergeResults) + geom_point(aes(x=Description, y=LogP), size=7) + 
  geom_bar( data = subset(mergeResults, `GroupID` == "aPD1"), aes(x=Description, y=LogP) , stat="Identity", width=0.03, color="#FFB18C", size=1) +
  geom_bar( data = subset(mergeResults, `GroupID` == "Healthy"), aes(x=Description, y=LogP) , stat="Identity", width=0.03, color="#7FAEDB", size=1) +
  coord_flip() + theme_bw() + ggtitle(title) + ylab("-Log (P)") + xlab(NULL) + 
  theme(axis.title.x = element_text(size=20), axis.text = element_text(size=20, color="black"), title = element_text(size=20)) + 
  # annotation_custom(left_grob) + annotation_custom(right_grob) + 
  scale_y_continuous(labels = abs, limits = c(-12,12)) + theme(panel.grid.major=element_blank()) + geom_hline(yintercept=0, size=0.5)

# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GeneOntology_PB_AvH_oW.pdf", width=5, height=8)
# write.csv(mergeResults, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3d.csv")



# ----------- PATHWAY analysis ----------- 


workingDir <- getwd()
setwd(dir = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/")

# ************ HiHi for AvH at oneWeek  ********************
pathwayPos <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/gsea_report_for_na_pos_1573243350115.xls", sep="\t")
pathwayNeg <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/gsea_report_for_na_neg_1573243350115.xls", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
plotGSEAlollipop( mergeResults, title = "Hallmark pathways", leftLabel = "Healthy", rightLabel = "aPD1")
plotGSEAlollipop( mergeResults, title = "Hallmark pathways", leftLabel = " ", rightLabel = " ", sizebyFDR = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HallmarkPathways_oW.pdf", device="pdf", width=7, height=9)
# write.csv(mergeResults, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3e.csv")



# IL-2 stat5 pathway 
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_IL2_STAT5_SIGNALING.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "IL2", title = "Hallmark: IL2-STAT5 pathway", leftLabel = "aPD1", rightLabel = "Healthy", 
                    cohortCompare =T, legendUpperRight = F)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "IL2", title = "Hallmark: IL2-STAT5 pathway", leftLabel = " ", rightLabel = " ", 
                    cohortCompare =T, legendUpperRight = F)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HallmarkPahtways_IL2STAT5_oW.pdf", device="pdf", width = 7, height=6)  
# write.csv(singlePathway, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3i.csv")

# TNF-NFkB
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "TNF", title = "Hallmark: TNF-NFkB pathway", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = F)

# IL6-STAT3
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_IL6_JAK_STAT3_SIGNALING.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "TNF", title = "Hallmark: IL6-STAT3 pathway", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = F)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "TNF", title = "Hallmark: IL6-STAT3 pathway", leftLabel = " ", rightLabel = " ", cohortCompare =T, legendUpperRight = F)

# G2m pathway 
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "G2M", title = "Hallmark: G2M checkpoint", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "G2M", title = "Hallmark: G2M checkpoint", leftLabel = " ", rightLabel = " ", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HallmarkPathways_G2M_oW.pdf", device="pdf", width=7, height=6)  
# write.csv(singlePathway, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e6b-1.csv")

# MitoticSpindle 
singlePathway <-  read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "MITOT", title = "Hallmark: Mitotic Spindle", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "MITOT", title = "Hallmark: Mitotic Spindle", leftLabel = " ", rightLabel = " ", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HallmarkPathways_MitSpindle_oW.pdf", device="pdf", width=7, height=6)    
# write.csv(singlePathway, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e6b-2.csv")

# ************ HiHi for AvH at baseline  ********************

pathwayPos <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_bL_hallmark.GseaPreranked.1573243362228/gsea_report_for_na_pos_1573243362228.xls", sep="\t")
pathwayNeg <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_bL_hallmark.GseaPreranked.1573243362228/gsea_report_for_na_neg_1573243362228.xls", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
plotGSEAlollipop( mergeResults, title = "Hallmark genesets: HiHi AvH at baseline", leftLabel = "Healthy", rightLabel = " aPD1")


# ************ PB for AvH at oneWeek  ********************
pathwayPos <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/gsea_report_for_na_pos_1574366663536.xls", sep="\t")
pathwayNeg <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/gsea_report_for_na_neg_1574366663536.xls", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
plotGSEAlollipop( mergeResults, title = "ASC at one week", leftLabel = "Healthy", rightLabel = "aPD1")
plotGSEAlollipop( mergeResults, title = "ASC at one week", leftLabel = " ", rightLabel = " ", sizebyFDR=T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/PB_HallmarkPathways_AvH_oW.pdf", device="pdf", width=6, height=4)
# write.csv(mergeResults, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3f.csv")

# G2m pathway 
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "G2M", title = "Hallmark: G2M checkpoint", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "G2M", title = "Hallmark: G2M checkpoint", leftLabel = " ", rightLabel = " ", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/PB_HallmarkPathways_G2M_oW.pdf", device="pdf", width=7, height=6)  
# write.csv(singlePathway, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e6b-3.csv")

# Mitotic Spindle
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "MITOT", title = "Hallmark: Mitotic Spindle", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "MITOT", title = "Hallmark: Mitotic Spindle", leftLabel = " ", rightLabel = " ", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/PB_HallmarkPathways_MitSpindle_oW.pdf", device="pdf", width=7, height=6)  
# write.csv(singlePathway, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e6b-4.csv")

# ************ ABC for AvH at oneWeek  ********************
pathwayPos <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/gsea_report_for_na_pos_1574696455529.xls", sep="\t")
pathwayNeg <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/gsea_report_for_na_neg_1574696455529.xls", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
plotGSEAlollipop( mergeResults, title = "ABC at one week", leftLabel = "Healthy", rightLabel = "aPD1", sizebyFDR = T) 
plotGSEAlollipop( mergeResults, title = "ABC at one week", leftLabel = " ", rightLabel = " ", sizebyFDR = T)  
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/ABC_HallmarkPathways_AvH_oW.pdf", device="pdf", width=7, height=4)
# write.csv(mergeResults, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e6a.csv")

# G2m pathway 
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "G2M", title = "Hallmark: G2M checkpoint", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "G2M", title = "Hallmark: G2M checkpoint", leftLabel = " ", rightLabel = " ", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/ABC_HallmarkPathways_G2M_oW.pdf", device="pdf", width=7, height=6)  
# write.csv(singlePathway, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e6b-5.csv")

# Mitotic Spindle
singlePathway <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "MITOT", title = "Hallmark: Mitotic Spindle", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "MITOT", title = "Hallmark: Mitotic Spindle", leftLabel = " ", rightLabel = " ", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/ABC_HallmarkPahtways_MitSpindle_oW.pdf", device="pdf", width=7, height=6)  
# write.csv(singlePathway, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e6b-6.csv")



# combined GSEA plots 

Tfh_G2M <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
Tfh_MitSpindle <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
PB_G2M <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
PB_MitSpindle <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
ABC_G2M <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
ABC_MitSpindle <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")

ggplot() + 
  geom_line(data=Tfh_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="darkgreen", size=1) + 
  geom_rug(data=Tfh_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="b", size=0.75, alpha=0.5, color="darkgreen") +
  geom_line(data=PB_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="orange", size=1) + 
  geom_rug(data=PB_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="t", size=0.75, alpha=0.5, color="orange") +
  geom_line(data=ABC_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="purple", size=1) + 
  #geom_rug(data=PB_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="t", size=0.75, alpha=0.5, color="orange") +
  theme_bw() +  ggtitle("Hallmark: Mitotic Spindle") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5))+
  geom_hline(yintercept = 0)  + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSEA_cTfh-plus-PB_Hallmark_AvH_oW_MitSpindle.pdf", device = 'pdf', width=7, height=6)

# write.csv(Tfh_MitSpindle, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3g-1.csv")
# write.csv(PB_MitSpindle, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3g-2.csv")
# write.csv(ABC_MitSpindle, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3g-3.csv")


ggplot() + 
  geom_line(data=Tfh_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="darkgreen", size=1) + 
  geom_rug(data=Tfh_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="b", size=0.75, alpha=0.5, color="darkgreen") +
  geom_line(data=PB_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="orange", size=1) + 
  geom_rug(data=PB_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="t", size=0.75, alpha=0.5, color="orange") +
  geom_line(data=ABC_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="purple", size=1) + 
  #geom_rug(data=ABC_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="t", size=0.75, alpha=0.5, color="orange") +
  theme_bw() +  ggtitle("Hallmark: G2M Checkpoint") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5))+
  geom_hline(yintercept = 0)  + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSEA_cTfh-plus-PB_Hallmark_AvH_oW_G2Mcheckpoint.pdf", device = 'pdf', width=7, height=6)
# write.csv(Tfh_G2M, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3h-1.csv")
# write.csv(PB_G2M, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3h-2.csv")
# write.csv(ABC_G2M, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/3h-3.csv")




setwd(workingDir)


# -----------  GSVA analysis ----------- 

gsets <- getGmt("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/h.all.v7.0.symbols.gmt")
HiHiv2 <- as.matrix( logDataMatrix[,grep("HiHi_oW",colnames(logDataMatrix))] )
GSVAhallmark <- gsva(HiHiv2, gsets, method="gsva")
GSVAhallmarkmeta <- makeMetaData( colnames(GSVAhallmark));  
temp <- as.data.frame(t(GSVAhallmark))
temp <- merge(x = temp, y = GSVAhallmarkmeta, by=0)
temp$Subject <- str_replace(temp$Subject,"^X","")
temp$Subject <- str_replace(temp$Subject, "[.]","-")
# row.names(temp) <- metaData[temp$Row.names,'Subject']

pheno_GSVAHiHi <- left_join(x = mergedData, y = temp, by = c("Cohort","Subject","TimeCategory"))
subsetData <- subset(pheno_GSVAHiHi, TimeCategory == "oneWeek")
univScatter(data=subsetData, yData = "HALLMARK_IL2_STAT5_SIGNALING", xData="cTfh_ICOShiCD38hi_..FreqParent", 
            fillParam = "TimeCategory", title = "IL-2 signaling", yLabel= "GSVA score - IL2-STAT5", 
            xLabel = expression(paste("ICO",S^'+', "CD3", 8^"+"," (% cTfh) at one week"))) + 
  coord_cartesian(xlim = c(0,12), ylim=c(-0.6,0.75)) + scale_x_continuous(breaks = c(seq(0,18,2)))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSVA_IL2STAT5_vs_cTfh_oW_univ.pdf", device="pdf", height=8, width=8)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","cTfh_ICOShiCD38hi_..FreqParent", "HALLMARK_IL2_STAT5_SIGNALING")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e6c-1.csv")

subsetData1 <- subset(pheno_GSVAHiHi, Cohort == "Healthy")    ; subsetData2 <- subset(pheno_GSVAHiHi, Cohort == "anti-PD-1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "anti-PD-1", yData = "HALLMARK_IL2_STAT5_SIGNALING", xData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "IL-2 signaling", yLabel= "GSVA score - IL2-STAT5", 
           xLabel = expression(paste("ICO",S^'+', "CD3", 8^"+"," (% cTfh) at one week"))) + 
  coord_cartesian(xlim = c(0,12), ylim=c(-0.6,0.75)) + scale_x_continuous(breaks = c(seq(0,18,2)))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSVA_IL2STAT5_vs_cTfh_oW_biv.pdf", device="pdf", height=8, width=8)
temp <- rbind(subsetData1, subsetData2)
# write.csv(temp[,c("dbCode","Cohort","TimeCategory","cTfh_ICOShiCD38hi_..FreqParent", "HALLMARK_IL2_STAT5_SIGNALING")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e6c-2.csv")


univScatter(data=subset(pheno_GSVAHiHi, TimeCategory == "oneWeek"), yData = "HALLMARK_MITOTIC_SPINDLE", xData="cTfh_ICOShiCD38hi_..FreqParent", 
            fillParam = "TimeCategory", title = "Mitotic Spindle", yLabel= "GSVA score - Mitotic Spindle", 
            xLabel = expression(paste("ICO",S^'+', "CD3", 8^"+"," (% cTfh) at one week"))) + 
  coord_cartesian(xlim=c(0,12),ylim = c(-1,1)) + scale_x_continuous(breaks = c(seq(0,18,2)))
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "anti-PD-1", yData = "HALLMARK_MITOTIC_SPINDLE", xData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "Mitotic Spindle", yLabel= "GSVA score - Mitotic spindle", 
           xLabel = expression(paste("ICO",S^'+', "CD3", 8^"+"," (% cTfh) at one week"))) + 
  coord_cartesian(xlim=c(0,12),ylim = c(-1,1)) + scale_x_continuous(breaks = c(seq(0,18,2)))

univScatter(data=subset(pheno_GSVAHiHi, TimeCategory == "oneWeek"), xData = "HALLMARK_MITOTIC_SPINDLE", yData="FCtfh_oW", 
            fillParam = "TimeCategory", title = "Mitotic Spindle", xLabel= "GSVA score - Mitotic Spindle", yLabel = expression(paste("ICO",S^'+', "CD3", 8^"+"," (% cTfh) at one week"))) + 
  coord_cartesian(ylim = c(0,10)) + scale_y_continuous(breaks = c(seq(0,18,2)))

bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "anti-PD-1", xData = "HALLMARK_MITOTIC_SPINDLE", yData="FCtfh_oW", 
           fillParam = "Cohort", title = "Mitotic Spindle", xLabel= "GSVA score - Mitotic spindle", yLabel = expression(paste("ICO",S^'+', "CD3", 8^"+"," (% cTfh) at one week"))) + 
  coord_cartesian(ylim = c(0,10)) + scale_y_continuous(breaks = c(seq(0,18,2)))




# -----------  link RNAseq to phenotypic analysis  ----------- 


PBv1 <- as.matrix( logDataMatrix[,grep("PB_bL",colnames(logDataMatrix))] )
GSVAhallmark <- gsva(PBv1, gsets, method="gsva")
GSVAhallmarkmeta <- makeMetaData( colnames(GSVAhallmark));  
temp <- as.data.frame(t(GSVAhallmark))
temp <- merge(x = temp, y = GSVAhallmarkmeta, by=0)
temp$Subject <- str_replace(temp$Subject,"^X","")
temp$Subject <- str_replace(temp$Subject, "[.]","-")

pheno_GSVAPB <- left_join(x = mergedData, y = temp, by = c("Cohort","Subject","TimeCategory"))

temp <- colnames(logDataMatrix)
temp <- str_replace(temp, "^X",""); temp <- str_replace(temp, "[.]","-"); temp <- str_replace(temp, "_S..",""); temp <- str_replace(temp, "_S.","")
temp <- str_replace(temp, "_bL","_baseline"); temp <- str_replace(temp, "_oW","_oneWeek");
temp2 <- logDataMatrix; colnames(temp2) <- temp; temp2 <- t(temp2)
a <- do.call(rbind.data.frame, strsplit(rownames(temp2), "_")) 
names(a) <- c("Subject", "Subset","TimeCategory")
temp2 <- as.data.frame(cbind(temp2, a))

phenoRNAHiHi <- right_join(x = subset(mergedData, Cohort != "nonPD1"), y = temp2[which(temp2$Subset == "HiHi"),], by = c("Subject","TimeCategory"), all.x = T )
phenoRNAHiHi$Cohort <- factor(phenoRNAHiHi$Cohort, levels = c("anti-PD-1","Healthy"))
phenoRNAASC <- right_join(x = mergedData, y = temp2[which(temp2$Subset == "PB"),], by = c("Subject","TimeCategory"), all.x = T )



# ----------- irAE analysis ----------- 


irAEonly <- subset(phenoRNAHiHi, !is.na(phenoRNAHiHi$irAE) & phenoRNAHiHi$irAE != "")
twoSampleBar(data=irAEonly, xData = "irAE", yData="PDCD1", fillParam = "irAE", title="PDCD1", yLabel="log counts") + 
  theme(plot.title = element_text(face='italic'), axis.title.x = element_text(size=28), axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0))+
  scale_fill_manual(values = c("grey90", "#ff9a6a"))+ coord_cartesian(ylim = c(2,7))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/irAE_PDCD1.pdf", width=4)
# write.csv(irAEonly[,c("irAE","PDCD1")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7b.csv")

twoSampleBar(data=irAEonly, xData = "irAE", yData="ICOS", fillParam = "irAE", title="ICOS", yLabel="log counts") + 
  theme(plot.title = element_text(face='italic'), axis.title.x = element_text(size=28), axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0))+
  scale_fill_manual(values = c("grey90", "#ff9a6a"))+ coord_cartesian(ylim = c(9.5,11)) 
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/irAE_ICOS.pdf", width=4)
# write.csv(irAEonly[,c("irAE","ICOS")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7c.csv")

twoSampleBar(data=irAEonly, xData = "irAE", yData="TFRC", fillParam = "irAE", title="TFRC", yLabel="log counts") + 
  theme(plot.title = element_text(face='italic'), axis.title.x = element_text(size=28), axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0))+
  scale_fill_manual(values = c("grey90", "#ff9a6a"))+ coord_cartesian(ylim = c(7,8.5))


# superseded by the heatmap that comes after the diffexp calculation

temp <- colnames(logDataMatrix)
temp <- str_replace(temp, "^X",""); temp <- str_replace(temp, "[.]","-"); temp <- str_replace(temp, "_S..",""); temp <- str_replace(temp, "_S.","")
temp <- str_replace(temp, "_bL","_baseline"); temp <- str_replace(temp, "_oW","_oneWeek");
temp2 <- logDataMatrix; colnames(temp2) <- temp;  logDataMatrix.renamed <- temp2
temp2 <- t(temp2)
a <- do.call(rbind.data.frame, strsplit(rownames(temp2), "_"))
names(a) <- c("Subject", "Subset","TimeCategory")
temp2 <- as.data.frame(cbind(temp2, a))
index <- paste0(substring(irAEonly$Label,first=0,last=9), "_HiHi_baseline")
logDataMatrix.irAE <- logDataMatrix.renamed[ , which(names(logDataMatrix.renamed) %in% index)]      # now logCounts of just the irAE-relevant Tfh samples

metaData.irAE <- data.frame(row.names=index)
metaData.irAE$irAE <- irAEonly[,c("irAE")];  metaData.irAE$irAE <- factor(metaData.irAE$irAE, levels = c("N","Y"))
probeGenes <- logDataMatrix.irAE[ c("ICOS","PDCD1","IL32","MKI67","BIRC5","XAB2","FOS","CD82","KIF2C","HELLS","CD52"),]
annotation <- metaData.irAE # [ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  irAE = c("N" ="grey90", "Y" = "#ff9a6a"))
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 24, color=inferno(100), main = "Log counts gene expression - HiHi Tfh ", show_colnames = F)

# ----  Differential expression analysis  --------

irAEdata <- dataMatrix[, which( names(logDataMatrix.renamed) %in% index)]       # since logDataMatrix.renamed is in same order as orig DataMatrix 
irAE.diffExp <- DESeqDataSetFromMatrix(countData=irAEdata, colData=metaData.irAE, design = ~irAE)
irAE.diffExp <- irAE.diffExp[idx,]   # keep hard filter from above for consistency
DESdata_irAE.diffExp <- DESeq(irAE.diffExp, parallel=TRUE)
irAE.diffExp.res <- as.data.frame(results(DESdata_irAE.diffExp, contrast = c("irAE","Y", "N" )))  # pos stats = first elem in comparison
# write.csv(irAE.diffExp.res, file = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/irAE_diffExp.csv")
diffExp<-data.frame( Up = nrow(irAE.diffExp.res[which(irAE.diffExp.res$padj < 0.05 & irAE.diffExp.res$stat>0),]), 
            Down = -nrow(irAE.diffExp.res[which(irAE.diffExp.res$padj < 0.05 & irAE.diffExp.res$stat<0),])); diffExp <- as.data.frame(t(diffExp))
diffExp$Labels <- row.names(diffExp); diffExp$Labels <- factor(diffExp$Labels, levels = c("Up","Down"))
ggplot(data = diffExp, aes(y=V1, x = Labels, fill=Labels)) + geom_bar(stat="identity", width=0.5) + theme_bw() + scale_fill_viridis_d("Direction") + 
  xlab(" ") + ylab("Number of genes") + ggtitle("irAE - Differentially\nexpressed genes")+ 
  theme(axis.text = element_text(color="black"), axis.text.x = element_text(size=24), axis.text.y = element_text(size=24), axis.title.y = element_text(size=24), 
        plot.title = element_text(size=32), legend.text = element_blank(), legend.title = element_blank() , legend.position = "none") + 
  scale_y_continuous(breaks = seq(-200,20,2)) + geom_hline(yintercept = 0, linetype = "solid", color="black")


temp <- colnames(logDataMatrix)
temp <- str_replace(temp, "^X",""); temp <- str_replace(temp, "[.]","-"); temp <- str_replace(temp, "_S..",""); temp <- str_replace(temp, "_S.","")
temp <- str_replace(temp, "_bL","_baseline"); temp <- str_replace(temp, "_oW","_oneWeek");
temp2 <- logDataMatrix; colnames(temp2) <- temp;  logDataMatrix.renamed <- temp2
temp2 <- t(temp2)
a <- do.call(rbind.data.frame, strsplit(rownames(temp2), "_")) 
names(a) <- c("Subject", "Subset","TimeCategory")
temp2 <- as.data.frame(cbind(temp2, a))
index <- paste0(substring(irAEonly$Label,first=0,last=9), "_HiHi_baseline")
logDataMatrix.irAE <- logDataMatrix.renamed[ , which(names(logDataMatrix.renamed) %in% index)]      # now logCounts of just the irAE-relevant Tfh samples

metaData.irAE <- data.frame(row.names=index)
metaData.irAE$irAE <- irAEonly[,c("irAE")];  metaData.irAE$irAE <- factor(metaData.irAE$irAE, levels = c("N","Y"))

diffExpGenes <- c(rownames(irAE.diffExp.res[which(irAE.diffExp.res$padj < 0.05 & irAE.diffExp.res$stat>0),]), 
                  rownames(irAE.diffExp.res[which(irAE.diffExp.res$padj < 0.05 & irAE.diffExp.res$stat<0),]))

newnames <- lapply(diffExpGenes, function(x) bquote(italic(.(x))))

probeGenes <- logDataMatrix.irAE[ diffExpGenes,]
annotation <- metaData.irAE # [ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  irAE = c("N" ="grey90", "Y" = "#ff9a6a"))
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 14, color=inferno(100), main = "Log counts gene expression - HiHi Tfh ", show_colnames = F,
         labels_row = as.expression(newnames),
         # filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/irAE_diffExp_heatmap.pdf"
)
dev.off()




volcanoPlot(irAE.diffExp.res, repelThresh = 0.05, title = "irAE in +/+ at baseline", leftLabel = "No irAE", rightLabel = "Yes irAE" )
volcanoPlot(irAE.diffExp.res, repelThresh = 0.05, title = "ICOS+CD38+ cTfh at baseline", leftLabel = " ", rightLabel = " " )
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/irAE_volcano.pdf")
# write.csv(irAE.diffExp.res, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7g.csv")
irAE.diffExp.res <- irAE.diffExp.res[ order(irAE.diffExp.res$stat, decreasing = T), ]; irAE.diffExp.res <-irAE.diffExp.res[ which(!is.na(irAE.diffExp.res$stat)), ]
# write.csv(irAE.diffExp.res, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/irAE.diffExp.res.csv")

pathwayPos <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/irAE_diffExp_Hallmark.GseaPreranked.1610904182341/gsea_report_for_na_pos_1610904182341.tsv", sep="\t")
pathwayNeg <- read.csv("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/Results/irAE_diffExp_Hallmark.GseaPreranked.1610904182341/gsea_report_for_na_neg_1610904182341.tsv", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
plotGSEAlollipop( mergeResults, title = "Hallmark pathways", leftLabel = " ", rightLabel = " ", sizebyFDR = T, colorRight = "#ff9a6a", colorLeft="grey90") + 
  theme(axis.text.y = element_text(angle=0, hjust=1,vjust=0.4, color="black", size=16), axis.text.x = element_text(angle=-90,hjust=1,vjust=0.4, color="black", size=16), 
        legend.text = element_text(size=16), legend.title = element_text(size=16)) + 
  scale_y_continuous(breaks=seq(-3,3,1), limits = c(-3,3), minor_breaks = NULL )
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/irAE_HallmarkPathways_baselineHiHi.pdf", device="pdf", width=7, height=10)
# write.csv(mergeResults, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/4b.csv")




