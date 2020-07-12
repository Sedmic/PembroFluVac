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
sessionInfo()
source('D:/Pembro-Fluvac/Analysis/Ranalysis_files/PembroFluSpec_PlottingFunctions.R')


## ***     read in data and make metaData  ************** 
dataMatrix <- read.table("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/PORTresults/NORMALIZED_DATA/SPREADSHEETS/FINAL_master_list_of_gene_counts_MIN.Fastqs.txt", sep="", 
                         header=T, stringsAsFactors = F)
rownames(dataMatrix) <- make.names(dataMatrix$geneSymbol, unique=TRUE)
dataMatrix$id <- dataMatrix$geneCoordinate <- dataMatrix$geneSymbol <- NULL 

makeMetaData <- function( colnVector)
{ metaData <- data.frame(row.names= colnVector);  metaData$condition <- "empty"
  metaData$Cohort <- metaData$Subject <- metaData$TimeCategory <- metaData$Subset <- "a"
  metaData$Cohort[grep("196",rownames(metaData),value=F)] <- "aPD1" ;  metaData$Cohort[grep("FS1819",rownames(metaData),value=F)] <- "Healthy"
  metaData$TimeCategory[grep("oW",rownames(metaData), value=F)] <- "oneWeek"; metaData$TimeCategory[grep("bL",rownames(metaData), value=F)] <- "baseline"; 
  metaData$Subset[grep("HiHi",rownames(metaData),value=F)] <- "HiHi_cTfh"; metaData$Subset[grep("ABC",rownames(metaData),value=F)] <- "ABC";
  metaData$Subset[grep("_PB_",rownames(metaData),value=F)] <- "PB"; metaData$Subset[grep("naiveB",rownames(metaData),value=F)] <- "navB";
  metaData$Subject <- substr(rownames(metaData), start = 1, stop = 10)
  metaData$condition <- paste0( metaData$Cohort,"_",metaData$Subset,"_",metaData$TimeCategory)
  return(metaData)      }

metaData <- makeMetaData(colnVector = colnames(dataMatrix))


#' ## ***     read in log-transformed data   **************

logDataMatrix <- read.csv( file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/logDataMatrix.csv", stringsAsFactors = F, header = 1, row.names = 1)  # created in QC script

#' ##  global heatmap of relevant/interesting genes  --------
probeGenes <- logDataMatrix[ c("CXCR5", "PRDM1",  "IFNG","MAF","CCR6","CXCR3","GATA3",
               "BTLA","TNFRSF4", "CD38","TIGIT","POU2AF1","MKI67","SH2D1A", "TOX2", "BIRC5", "TBK1", "CD19","MS4A1","FUT8","ST6GAL1","B4GALT1",
               "JCHAIN","AICDA", "CMTM6", "AK9", "ZNF770","CD3E","CD4", "PAX5", "SDC1") , ]  # some genes taken from genomicscape.com


annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression", show_colnames = F, cutree_cols = 5, cutree_rows = 4
         #, filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/SelectedGenes_allSubsets_allTimePoints_allSubjects_heatmap.pdf"
         )



# ---- glycosylation genes in ABC --------
probeGenes <- logDataMatrix[ c("FUT8","ST6GAL1","B4GALT1","AICDA") , grep("ABC",colnames(logDataMatrix)) ]
annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression - ABC ")

metaX <- makeMetaData( colnames(probeGenes))
singleGeneData <- merge(x = metaX, y = t(probeGenes[ "FUT8",]), by=0); colnames(singleGeneData)[7] <- "value"
prePostTimeGene(singleGeneData, xData = "TimeCategory" , yData = "value", fillParam = "Cohort", groupby = "Subject", title = "FUT8 in ABC", xLabel = "TimeCategory", yLabel = "log2 counts")
prePostTimeAveragedGene(singleGeneData, title = "ABC: FUT8 transcripts", xLabel = "TimeCategory", yLabel = "log2 counts")
  
singleGeneData <- merge(x = metaX, y = t(probeGenes[ "ST6GAL1",]), by=0); colnames(singleGeneData)[7] <- "value"
prePostTimeGene(singleGeneData, xData = "TimeCategory" , yData = "value", fillParam = "Cohort", groupby = "Subject", title = "ST6GAL1 in ABC", xLabel = "TimeCategory", yLabel = "log2 counts")
prePostTimeAveragedGene(singleGeneData, title = "ABC: ST6GAL1 transcripts", xLabel = "TimeCategory", yLabel = "log2 counts")

singleGeneData <- merge(x = metaX, y = t(probeGenes[ "B4GALT1",]), by=0); colnames(singleGeneData)[7] <- "value"
prePostTimeGene(singleGeneData, xData = "TimeCategory" , yData = "value", fillParam = "Cohort", groupby = "Subject", title = "B4GALT1 in ABC", xLabel = "TimeCategory", yLabel = "log2 counts")
prePostTimeAveragedGene(singleGeneData, title = "ABC: B4GALT1 transcripts", xLabel = "TimeCategory", yLabel = "log2 counts")

# ---- glycosylation genes in PB --------
annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
probeGenes <- logDataMatrix[ c("FUT8","ST6GAL1","B4GALT1") , grep("PB",colnames(logDataMatrix)) ]
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression - PB ")

metaX <- makeMetaData( colnames(probeGenes))
singleGeneData <- merge(x = metaX, y = t(probeGenes[ "FUT8",]), by=0); colnames(singleGeneData)[7] <- "value"
prePostTimeGene(singleGeneData, xData = "TimeCategory" , yData = "value", fillParam = "Cohort", groupby = "Subject", title = "FUT8 in ABC", xLabel = "TimeCategory", yLabel = "log2 counts")
prePostTimeAveragedGene(singleGeneData, title = "PB: FUT8 transcripts", xLabel = "TimeCategory", yLabel = "log2 counts")

singleGeneData <- merge(x = metaX, y = t(probeGenes[ "ST6GAL1",]), by=0); colnames(singleGeneData)[7] <- "value"
prePostTimeGene(singleGeneData, xData = "TimeCategory" , yData = "value", fillParam = "Cohort", groupby = "Subject", title = "ST6GAL1 in ABC", xLabel = "TimeCategory", yLabel = "log2 counts")
prePostTimeAveragedGene(singleGeneData, title = "PB: ST6GAL1 transcripts", xLabel = "TimeCategory", yLabel = "log2 counts")

singleGeneData <- merge(x = metaX, y = t(probeGenes[ "B4GALT1",]), by=0); colnames(singleGeneData)[7] <- "value"
prePostTimeGene(singleGeneData, xData = "TimeCategory" , yData = "value", fillParam = "Cohort", groupby = "Subject", title = "B4GALT1 in ABC", xLabel = "TimeCategory", yLabel = "log2 counts")
prePostTimeAveragedGene(singleGeneData, title = "PB: B4GALT1 transcripts", xLabel = "TimeCategory", yLabel = "log2 counts")


# ---- glycosylation genes in NaiveB  --------
annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
probeGenes <- logDataMatrix[ c("FUT8","ST6GAL1","B4GALT1") , grep("naive",colnames(logDataMatrix)) ]
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression - PB ")

metaX <- makeMetaData( colnames(probeGenes))
singleGeneData <- merge(x = metaX, y = t(probeGenes[ "FUT8",]), by=0); colnames(singleGeneData)[7] <- "value"
prePostTimeAveragedGene(singleGeneData, title = "NaiveB: FUT8 transcripts", xLabel = "TimeCategory", yLabel = "log2 counts")

singleGeneData <- merge(x = metaX, y = t(probeGenes[ "ST6GAL1",]), by=0); colnames(singleGeneData)[7] <- "value"
prePostTimeAveragedGene(singleGeneData, title = "NaiveB: ST6GAL1 transcripts", xLabel = "TimeCategory", yLabel = "log2 counts")

singleGeneData <- merge(x = metaX, y = t(probeGenes[ "B4GALT1",]), by=0); colnames(singleGeneData)[7] <- "value"
prePostTimeAveragedGene(singleGeneData, title = "NaiveB: B4GALT1 transcripts", xLabel = "TimeCategory", yLabel = "log2 counts")





## ******** Differential expression by subgroup and age category *********************

library("BiocParallel")
register(SnowParam(10))
metaData <- makeMetaData(colnVector = colnames(dataMatrix))
fullDataset <- DESeqDataSetFromMatrix(countData=dataMatrix, colData=metaData, design= ~ condition)
dds <- estimateSizeFactors(fullDataset);   idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) > 20  # hard filter for genes with at least 20 counts in 25% of samples
fullDataset <- fullDataset[idx,]
# DESdata_fullDataset <- DESeq(fullDataset, parallel=TRUE)
# approach of using full metaData with ~ condition keeps failing so will switch to simpler statistical models
# dataMatrixFiltered <- dataMatrix[idx,]


# ----------- NavB in aPD1 vs HC baseline 

navB_aPD1_v_HC_bLmeta <- metaData[ grep("naiveB_bL", rownames(metaData)),1:2]
navB_aPD1_v_HC_bLmeta$condition <- factor(navB_aPD1_v_HC_bLmeta$condition);  navB_aPD1_v_HC_bLdata <- dataMatrix[ , grep("naiveB_bL", colnames(dataMatrix) ) ]
navB_aPD1_v_HC_bL <- DESeqDataSetFromMatrix(countData=navB_aPD1_v_HC_bLdata, colData=navB_aPD1_v_HC_bLmeta, design = ~condition)
navB_aPD1_v_HC_bL <- navB_aPD1_v_HC_bL[idx,]   # keep hard filter from above for consistency
DESdata_navB_AvH_bL <- DESeq(navB_aPD1_v_HC_bL, parallel=TRUE)

navB_AvH_bL <- as.data.frame(results(DESdata_navB_AvH_bL, contrast = c("condition", "aPD1_navB_baseline", "Healthy_navB_baseline") ))  # pos stats = first elem in comparison
volcanoPlot(navB_AvH_bL, repelThresh = 0.01, title = "navB at oneWeek in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )

probeGenes <- logDataMatrix[ grep(paste(c("RANGAP1$", "PLEC$","^OLR1$","SEMA6B","KIFC3$","IRAK1$","^AES$","^ENG$","VAT1$","CD82$",
                                          "BEX4$","ALKBH4"), collapse="|"), rownames(logDataMatrix), value=F), grep("naive",colnames(logDataMatrix),value=F)]

annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression")


# ----------- NavB in aPD1 vs HC oneWeek 

navB_aPD1_v_HC_oWmeta <- metaData[ grep("naiveB_oW", rownames(metaData)),1:2]
navB_aPD1_v_HC_oWmeta$condition <- factor(navB_aPD1_v_HC_oWmeta$condition);  navB_aPD1_v_HC_oWdata <- dataMatrix[ , grep("naiveB_oW", colnames(dataMatrix) ) ]
navB_aPD1_v_HC_oW <- DESeqDataSetFromMatrix(countData=navB_aPD1_v_HC_oWdata, colData=navB_aPD1_v_HC_oWmeta, design = ~condition)
navB_aPD1_v_HC_oW <- navB_aPD1_v_HC_oW[idx,]   # keep hard filter from above for consistency
DESdata_navB_AvH_oW <- DESeq(navB_aPD1_v_HC_oW, parallel=TRUE)

navB_AvH_oW <- as.data.frame(results(DESdata_navB_AvH_oW, contrast = c("condition", "aPD1_navB_oneWeek", "Healthy_navB_oneWeek") ))  # pos stats = first elem in comparison
volcanoPlot(navB_AvH_oW, repelThresh = 0.01, title = "navB at oneWeek in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )

probeGenes <- logDataMatrix[ grep(paste(c("EGR3$", "CCNL1$","GZMB$","C3AR1","FCRL6$","QRICH2$","^LCK$","^CD5$","OAS1$","RICTOR$",
                                          "UBR4"), collapse="|"), rownames(logDataMatrix), value=F), grep("naive",colnames(logDataMatrix),value=F)]

annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression")


# ----------- NavB in HC at Baseline vs oneWeek 
navB_HC_bL_v_oWmeta <- metaData[ grep("Healthy_navB", metaData$condition),c(1:2,4)]
navB_HC_bL_v_oWmeta$condition <- factor(navB_HC_bL_v_oWmeta$condition);  navB_HC_bL_v_oWdata <- dataMatrix[ , grep("naiveB", colnames(dataMatrix) ) ]; navB_HC_bL_v_oWdata <- navB_HC_bL_v_oWdata[, grep("FS",colnames(navB_HC_bL_v_oWdata))]
navB_HC_bL_v_oW <- DESeqDataSetFromMatrix(countData=navB_HC_bL_v_oWdata, colData=navB_HC_bL_v_oWmeta, design = ~condition + Subject)   # control for Subject
navB_HC_bL_v_oW <- navB_HC_bL_v_oW[idx,]   # keep hard filter from above for consistency
DESdata_navB_H_bLvoW <- DESeq(navB_HC_bL_v_oW, parallel=TRUE)

navB_H_bL_v_oW <- as.data.frame(results(DESdata_navB_H_bLvoW, contrast = c("condition", "Healthy_navB_oneWeek", "Healthy_navB_baseline") ))  # pos stats = first elem in comparison
# write.csv(navB_H_bL_v_oW, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/HiHi_AvH_oW.csv")
volcanoPlot(navB_H_bL_v_oW, repelThresh = 0.05, title = "NavB in Healthy at bL vs oW", leftLabel = "Baseline", rightLabel = "oneWeek" )

probeGenes <- logDataMatrix[ grep(paste(c("KIAA1683", "GPER1","NR4A3","PMEPA1","CD83","^ID1$","^FOS$","^JUN$","SIGLEC1$","^IER2$",
                                          "^FAM83D$","NR4A2","^ID3$","KLHDC8B"), collapse="|"), rownames(logDataMatrix), value=F), grep("naive",colnames(logDataMatrix),value=F)]

annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression")


# ----------- NavB in aPD1 at Baseline vs oneWeek 
navB_A_bL_voWmeta <- metaData[ grep("aPD1_navB", metaData$condition),c(1:2,4)]
navB_A_bL_voWmeta$condition <- factor(navB_A_bL_voWmeta$condition);  navB_A_bL_voWdata <- dataMatrix[ , grep("naiveB", colnames(dataMatrix) ) ]; navB_A_bL_voWdata <- navB_A_bL_voWdata[, grep("19616",colnames(navB_A_bL_voWdata))]
navB_A_bL_voW <- DESeqDataSetFromMatrix(countData=navB_A_bL_voWdata, colData=navB_A_bL_voWmeta, design = ~condition)   # cannot control for Subject
navB_A_bL_voW <- navB_A_bL_voW[idx,]   # keep hard filter from above for consistency
DESdata_navB_A_bLvoW <- DESeq(navB_A_bL_voW, parallel=TRUE)

navB_A_bL_v_oW <- as.data.frame(results(DESdata_navB_A_bLvoW, contrast = c("condition", "aPD1_navB_oneWeek", "aPD1_navB_baseline") ))  # pos stats = first elem in comparison
volcanoPlot(navB_A_bL_v_oW, repelThresh = 0.05, title = "NavB in aPD1 at bL vs oW", leftLabel = "Baseline", rightLabel = "oneWeek" )

probeGenes <- logDataMatrix[ grep(paste(c("KIAA1683", "GPER1","NR4A3","PMEPA1","CD83","^ID1$","^FOS$","^JUN$","SIGLEC1$","^IER2$",
                                          "^FAM83D$","NR4A2","^ID3$","KLHDC8B"), collapse="|"), rownames(logDataMatrix), value=F), grep("naive",colnames(logDataMatrix),value=F)]

annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression")




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

HiHi_AvH_bL <- as.data.frame(results(DESdata_HiHi_AvH_bL, contrast = c("condition", "aPD1_HiHi_cTfh_baseline", "Healthy_HiHi_cTfh_baseline") ))  # pos stats = first elem in comparison
# write.csv(HiHi_AvH_bL, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/HiHi_AvH_bL.csv")
volcanoPlot(HiHi_AvH_bL, repelThresh = 0.10, title = "HiHi at oneWeek in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )

probeGenes <- logDataMatrix[ grep(paste(c("TACC3$", "TXNIP","BIN2$","EREG$","GNG4$","^IL1R1$","^FYB$","IL27RA","FOXP3","FGFR1$",
                                          "CTLA4$","RASA3$","TUBB1$","^BCL2$"), collapse="|"), rownames(logDataMatrix), value=F), grep("HiHi",colnames(logDataMatrix),value=F)]

annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression")



# ----------- HiHi in aPD1 vs HC   oneWeek 

HiHi_aPD1_v_HC_oWmeta <- metaData[ grep("HiHi_oW", rownames(metaData)),1:2]
HiHi_aPD1_v_HC_oWmeta$condition <- factor(HiHi_aPD1_v_HC_oWmeta$condition);  HiHi_aPD1_v_HC_oWdata <- dataMatrix[ , grep("HiHi_oW", colnames(dataMatrix) ) ]
HiHi_aPD1_v_HC_oW <- DESeqDataSetFromMatrix(countData=HiHi_aPD1_v_HC_oWdata, colData=HiHi_aPD1_v_HC_oWmeta, design = ~condition)
HiHi_aPD1_v_HC_oW <- HiHi_aPD1_v_HC_oW[idx,]   # keep hard filter from above for consistency
DESdata_HiHi_AvH_oW <- DESeq(HiHi_aPD1_v_HC_oW, parallel=TRUE)

HiHi_AvH_oW <- as.data.frame(results(DESdata_HiHi_AvH_oW, contrast = c("condition", "aPD1_HiHi_cTfh_oneWeek", "Healthy_HiHi_cTfh_oneWeek") ))  # pos stats = first elem in comparison
# write.csv(HiHi_AvH_oW, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/HiHi_AvH_oW.csv")
volcanoPlot(HiHi_AvH_oW, repelThresh = 0.15, title = "HiHi at oneWeek in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/HiHi_H-v-aPD1_d7_volcano.pdf", device="pdf")

probeGenes <- logDataMatrix[ grep(paste(c("MKI67$", "HELLS","GGCT","PMM1","PCSK7","^GK$","^OASL$","TNFRSF1B","OASL","UGDH$",
                                          "TMEM71","RALBP1","ORAI2","CKAP2L", "PRDM1$"), collapse="|"), rownames(logDataMatrix), value=F), grep("HiHi",colnames(logDataMatrix),value=F)]

annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression", border_color = NA
#        , filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HvA_diffExp_heatmap.pdf", device="pdf"
         )
dev.off()

metaX <- makeMetaData(colnames(probeGenes)); 
probeGenes <- logDataMatrix[ grep(paste(c("BCL6$","IL21","IL10", "PRDM1$","TOX2"), collapse="|"), rownames(logDataMatrix), value=F), grep("HiHi",colnames(logDataMatrix),value=F)]

data <- merge(metaX, y=t(probeGenes), by=0) ; data <- data[data$TimeCategory == "oneWeek", ]
data$Cohort <- factor(data$Cohort, levels=c("Healthy","aPD1"))
twoSampleBox(data[, c("TimeCategory","Cohort","BCL6")], xData = "Cohort", yData = "BCL6", fillParam = "Cohort", title = "BCL6 in +/+ cTfh", yLabel = "log2 counts")
twoSampleBox(data[, c("TimeCategory","Cohort","PRDM1")], xData = "Cohort", yData = "PRDM1", fillParam = "Cohort", title = "PRDM1 in +/+ cTfh", yLabel = "log2 counts")
twoSampleBox(data[, c("TimeCategory","Cohort","IL10")], xData = "Cohort", yData = "IL10", fillParam = "Cohort", title = "IL10 in +/+ cTfh", yLabel = "log2 counts")
twoSampleBox(data[, c("TimeCategory","Cohort","TOX2")], xData = "Cohort", yData = "TOX2", fillParam = "Cohort", title = "TOX2 in +/+ cTfh", yLabel = "log2 counts")

data$BCL6_to_PRDM1_ratio <- data$BCL6 / data$PRDM1
twoSampleBox(data[, c("TimeCategory","Cohort","BCL6_to_PRDM1_ratio")], xData = "Cohort", yData = "BCL6_to_PRDM1_ratio", fillParam = "Cohort", 
             title = "BCL6_to_PRDM1_ratio in +/+ cTfh", yLabel = "Bcl6/PRDM1 counts")


# ----------- PB in aPD1 vs HC   baseline  

PB_aPD1_v_HC_bLmeta <- metaData[ grep("PB_bL", rownames(metaData)),1:2]
PB_aPD1_v_HC_bLmeta$condition <- factor(PB_aPD1_v_HC_bLmeta$condition);  PB_aPD1_v_HC_bLdata <- dataMatrix[ , grep("PB_bL", colnames(dataMatrix) ) ]
PB_aPD1_v_HC_bL <- DESeqDataSetFromMatrix(countData=PB_aPD1_v_HC_bLdata, colData=PB_aPD1_v_HC_bLmeta, design = ~condition)
PB_aPD1_v_HC_bL <- PB_aPD1_v_HC_bL[idx,]   # keep hard filter from above for consistency
DESdata_PB_AvH_bL <- DESeq(PB_aPD1_v_HC_bL, parallel=TRUE)

PB_AvH_bL <- as.data.frame(results(DESdata_PB_AvH_bL, contrast = c("condition", "aPD1_PB_baseline", "Healthy_PB_baseline") ))  # pos stats = first elem in comparison
# write.csv(PB_AvH_bL, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/PB_AvH_bL.csv")
volcanoPlot(PB_AvH_bL, repelThresh = 0.15, title = "PB at baseline in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )

probeGenes <- logDataMatrix[ grep(paste(c("ZNF44$", "^SDPR$","CH25H$","CXCR1$","IFIT2$","PEG13$","^FOS$","UNC119$","DOCK9$","^SNN$",
                                          "S1PR5","IGHV2.70.1"), collapse="|"), rownames(logDataMatrix), value=F), grep("PB",colnames(logDataMatrix),value=F)]

annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression")

# ----------- PB in aPD1 vs HC   oneWeek 

PB_aPD1_v_HC_oWmeta <- metaData[ grep("PB_oW", rownames(metaData)),1:2]
PB_aPD1_v_HC_oWmeta$condition <- factor(PB_aPD1_v_HC_oWmeta$condition);  PB_aPD1_v_HC_oWdata <- dataMatrix[ , grep("PB_oW", colnames(dataMatrix) ) ]
PB_aPD1_v_HC_oW <- DESeqDataSetFromMatrix(countData=PB_aPD1_v_HC_oWdata, colData=PB_aPD1_v_HC_oWmeta, design = ~condition)
PB_aPD1_v_HC_oW <- PB_aPD1_v_HC_oW[idx,]   # keep hard filter from above for consistency
DESdata_PB_AvH_oW <- DESeq(PB_aPD1_v_HC_oW, parallel=TRUE)

PB_AvH_oW <- as.data.frame(results(DESdata_PB_AvH_oW, contrast = c("condition", "aPD1_PB_oneWeek", "Healthy_PB_oneWeek") ))  # pos stats = first elem in comparison
# write.csv(PB_AvH_oW, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/PB_AvH_oW.csv")
volcanoPlot(PB_AvH_oW, repelThresh = 0.15, title = "PB at oneWeek in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )

probeGenes <- logDataMatrix[ grep(paste(c("CD6$", "SLBP$","GNAI2$","STK17A","PNRC1$","DYRK2$","^DUSP1$","KLF3$","B3GNT2$","HEBP2$",
                                          "PNRC1","SNX20","IGLV9.49", "KIF2C","AURKB"), collapse="|"), rownames(logDataMatrix), value=F), grep("PB",colnames(logDataMatrix),value=F)]

annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression")



# ----------- ABC in aPD1 vs HC   oneWeek 

ABC_aPD1_v_HC_oWmeta <- metaData[ grep("ABC_oW", rownames(metaData)),1:2]
ABC_aPD1_v_HC_oWmeta$condition <- factor(ABC_aPD1_v_HC_oWmeta$condition);  ABC_aPD1_v_HC_oWdata <- dataMatrix[ , grep("ABC_oW", colnames(dataMatrix) ) ]
ABC_aPD1_v_HC_oW <- DESeqDataSetFromMatrix(countData=ABC_aPD1_v_HC_oWdata, colData=ABC_aPD1_v_HC_oWmeta, design = ~condition)
ABC_aPD1_v_HC_oW <- ABC_aPD1_v_HC_oW[idx,]   # keep hard filter from above for consistency
DESdata_ABC_AvH_oW <- DESeq(ABC_aPD1_v_HC_oW, parallel=TRUE)

ABC_AvH_oW <- as.data.frame(results(DESdata_ABC_AvH_oW, contrast = c("condition", "aPD1_ABC_oneWeek", "Healthy_ABC_oneWeek") ))  # pos stats = first elem in comparison
# write.csv(ABC_AvH_oW, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/ABC_AvH_oW.csv")
volcanoPlot(ABC_AvH_oW, repelThresh = 0.05, title = "ABC at oneWeek in aPD1 vs Healthy", leftLabel = "Healthy", rightLabel = "aPD1" )

probeGenes <- logDataMatrix[ grep(paste(c("RPP40$", "ZSCAN22$","NTHL1$","RARRES3","MRPL10$","MGEA5$","^CCNL2$","NKTR$","COG6$","DDX39A$",
                                          "GOLGB1","SREK1","IGKV2.29"), collapse="|"), rownames(logDataMatrix), value=F), grep("ABC",colnames(logDataMatrix),value=F)]

annotation <- metaData[ , -grep(paste(c("condition","Subject"),collapse="|"),colnames(metaData),value=F)]
ann_colors = list(  Cohort = c("Healthy" ="#7FAEDB", "aPD1" = "#FFB18C"), Subset = c("HiHi_cTfh"="#aaccbb", "ABC"="blue", "PB"="yellow", "navB"="orange"), TimeCategory = c("baseline"="grey90","oneWeek"="grey40")  )
pheatmap(probeGenes, scale="row", cluster_col=T, cluster_row=T, annotation_col = annotation, annotation_colors= ann_colors,
         fontsize_row = 12, color=inferno(100), main = "Log counts gene expression")


# ----------- ----------- 

#                                 PATHWAY analysis 

# ----------- ----------- 


workingDir <- getwd()
setwd(dir = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/")

# ************ HiHi for AvH at oneWeek  ********************
pathwayPos <- read.csv("Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/gsea_report_for_na_pos_1573243350115.xls", sep="\t")
pathwayNeg <- read.csv("Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/gsea_report_for_na_neg_1573243350115.xls", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
plotGSEAlollipop( mergeResults, title = "Hallmark: HiHi AvH at oW", leftLabel = "Healthy", rightLabel = "aPD1")
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HallmarkPathways_oW.pdf", device="pdf", width=7, height=9)

# IL-2 stat5 pathway 
singlePathway <- read.csv("Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_IL2_STAT5_SIGNALING.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "IL2", title = "IL2-STAT5: HiHi AvH at oW", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = F)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HallmarkPahtways_IL2STAT5_oW.pdf", device="pdf", width = 7, height=6)  

# TNF-NFkB
singlePathway <- read.csv("Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_TNFA_SIGNALING_VIA_NFKB.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "TNF", title = "TNF-NFkB: HiHi AvH at oW", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = F)

# G2m pathway 
singlePathway <- read.csv("Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "G2M", title = "G2M checkpoint: HiHi AvH at oW", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HallmarkPathways_G2M_oW.pdf", device="pdf", width=7, height=6)  

# MitoticSpindle 
singlePathway <-  read.csv("Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "MITOT", title = "Mitotic Spindle: HiHi AvH at oW", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/cTfh_HallmarkPathways_MitSpindle_oW.pdf", device="pdf", width=7, height=6)    


# ************ HiHi for AvH at baseline  ********************

pathwayPos <- read.csv("Results/HiHi_AvH_bL_hallmark.GseaPreranked.1573243362228/gsea_report_for_na_pos_1573243362228.xls", sep="\t")
pathwayNeg <- read.csv("Results/HiHi_AvH_bL_hallmark.GseaPreranked.1573243362228/gsea_report_for_na_neg_1573243362228.xls", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
plotGSEAlollipop( mergeResults, title = "Hallmark genesets: HiHi AvH at baseline", leftLabel = "Healthy", rightLabel = " aPD1")


# ************ PB for AvH at oneWeek  ********************
pathwayPos <- read.csv("Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/gsea_report_for_na_pos_1574366663536.xls", sep="\t")
pathwayNeg <- read.csv("Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/gsea_report_for_na_neg_1574366663536.xls", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
plotGSEAlollipop( mergeResults, title = "Hallmark: PB AvH at oW", leftLabel = "Healthy", rightLabel = "aPD1")
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/PB_HallmarkPathways_AvH_oW.pdf", device="pdf", width=6, height=7)

# G2m pathway 
singlePathway <- read.csv("Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "G2M", title = "G2M checkpoint: PB AvH at oW", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/PB_HallmarkPathways_G2M_oW.pdf", device="pdf", width=7, height=6)  

# Mitotic Spindle
singlePathway <- read.csv("Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "MITOT", title = "Mitotic Spindle: PB AvH at oW", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/PB_HallmarkPahtways_MitSpindle_oW.pdf", device="pdf", width=7, height=6)  


# ************ ABC for AvH at oneWeek  ********************
pathwayPos <- read.csv("Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/gsea_report_for_na_pos_1574696455529.xls", sep="\t")
pathwayNeg <- read.csv("Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/gsea_report_for_na_neg_1574696455529.xls", sep="\t")
mergeResults <- rbind(pathwayPos, pathwayNeg)
plotGSEAlollipop( mergeResults, title = "Hallmark: ABC AvH at oW", leftLabel = "Healthy", rightLabel = "aPD1")
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/ABC_HallmarkPathways_AvH_oW.pdf", device="pdf", width=7, height=6)

# G2m pathway 
singlePathway <- read.csv("Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "G2M", title = "G2M checkpoint: ABC AvH at oW", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/ABC_HallmarkPathways_G2M_oW.pdf", device="pdf", width=7, height=6)  

# Mitotic Spindle
singlePathway <- read.csv("Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
plotGSEAtraditional(mergeResults, singlePathway, pathwayName = "MITOT", title = "Mitotic Spindle: ABC AvH at oW", leftLabel = "aPD1", rightLabel = "Healthy", cohortCompare =T, legendUpperRight = T)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/ABC_HallmarkPahtways_MitSpindle_oW.pdf", device="pdf", width=7, height=6)  





# combined GSEA plots 

Tfh_G2M <- read.csv("Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
Tfh_MitSpindle <- read.csv("Results/HiHi_AvH_oW_hallmark.GseaPreranked.1573243350115/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
PB_G2M <- read.csv("Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
PB_MitSpindle <- read.csv("Results/PB_AvH_oW_hallmark.GseaPreranked.1574366663536/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")
ABC_G2M <- read.csv("Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/HALLMARK_G2M_CHECKPOINT.xls", sep="\t")
ABC_MitSpindle <- read.csv("Results/ABC_AvH_oW_hallmark.GseaPreranked.1574696455529/HALLMARK_MITOTIC_SPINDLE.xls", sep="\t")

ggplot() + 
  geom_line(data=Tfh_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="darkgreen", size=1) + 
  geom_rug(data=Tfh_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="b", size=0.75, alpha=0.5, color="darkgreen") +
  geom_line(data=PB_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="orange", size=1) + 
  geom_rug(data=PB_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="t", size=0.75, alpha=0.5, color="orange") +
  geom_line(data=ABC_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="purple", size=1) + 
  #geom_rug(data=ABC_G2M, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="t", size=0.75, alpha=0.5, color="orange") +
  theme_bw() +  ggtitle("Hallmark: G2M Checkpoint") + ylab("Enrichment score") + xlab("Rank in gene list") + 
  theme(axis.text = element_text(size=24,hjust = 0.5), axis.title = element_text(size=24,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5))+
  geom_hline(yintercept = 0)  
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSEA_cTfh-plus-PB_Hallmark_AvH_oW_G2Mcheckpoint.pdf", device = 'pdf', width=7, height=6)
  
ggplot() + 
    geom_line(data=Tfh_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="darkgreen", size=1) + 
    geom_rug(data=Tfh_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="b", size=0.75, alpha=0.5, color="darkgreen") +
    geom_line(data=PB_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="orange", size=1) + 
    geom_rug(data=PB_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="t", size=0.75, alpha=0.5, color="orange") +
    geom_line(data=ABC_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), color="purple", size=1) + 
    #geom_rug(data=PB_MitSpindle, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES), sides="t", size=0.75, alpha=0.5, color="orange") +
    theme_bw() +  ggtitle("Hallmark: Mitotic Spindle") + ylab("Enrichment score") + xlab("Rank in gene list") + 
    theme(axis.text = element_text(size=24,hjust = 0.5), axis.title = element_text(size=24,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5))+
    geom_hline(yintercept = 0)  
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSEA_cTfh-plus-PB_Hallmark_AvH_oW_MitSpindle.pdf", device = 'pdf', width=7, height=6)

  


setwd(workingDir)



# ----------- ----------- 

#                             INGENUITY PATHWAY analysis 

# ----------- ----------- 



workingDir <- getwd()
setwd(dir = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/IPA/")

pathwaysIPA <- read.csv("Reports/HiHi_AvH_oW_canonical.txt", sep="\t",header=T)    # works if you delete the Qiagen line at top
pathwaysIPA <- pathwaysIPA[order(pathwaysIPA$z.score, decreasing=F),]  # positive z scores mean enrichment for the POSITIVE stat from DiffExp so aPD1 here

pathwaysIPA_filter <- pathwaysIPA[ which(pathwaysIPA$z.score < -2.5 | pathwaysIPA$z.score > 1), -5]
plotIPAlollipop( pathwaysIPA_filter, title = "Ingenuity Canonical: HiHi AvH oneWeek", leftLabel = "Enrich for Healthy", rightLabel = "Enrich for aPD1")
  
upstreamIPA <- read.csv("Reports/HiHi_AvH_oW_upstream.txt", sep="\t", header=T)
upstreamIPA_filter <- upstreamIPA[ which(upstreamIPA$Activation.z.score < -3.8 | upstreamIPA$Activation.z.score > 3), -c(8:9)]
names(upstreamIPA_filter)[5] <- "z.score"
plotIPAlollipop( upstreamIPA_filter, title = "Ingenuity UpstreamPred: HiHi AvH oneWeek", leftLabel = "Predict for Healthy", rightLabel = "Predict for aPD1")


setwd(workingDir)

# ----------- ----------- 

#' ## -                           GSVA analysis 

# ----------- ----------- 


gsets <- getGmt("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/differentialExpression/GSEA/h.all.v7.0.symbols.gmt")
HiHiv2 <- as.matrix( logDataMatrix[,grep("HiHi_oW",colnames(logDataMatrix))] )
GSVAhallmark <- gsva(HiHiv2, gsets, method="gsva")
GSVAhallmarkmeta <- makeMetaData( colnames(GSVAhallmark))
temp <- as.data.frame(t(GSVAhallmark))
temp <- merge(x = temp, y = GSVAhallmarkmeta, by=0)
temp$Subject <- str_replace(temp$Subject,"^X","")
temp$Subject <- str_replace(temp$Subject, "[.]","-")
loadData <- function()
{
  mergedData <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/allMergedData.csv", stringsAsFactors = F, header = T, row.names = 1)
  mergedData$TimeCategory <- factor(mergedData$TimeCategory, levels = c("baseline", "oneWeek","late"))
  mergedData$Cohort <- factor(mergedData$Cohort, levels = c("Healthy", "aPD1","nonPD1"))
  mergedData$dummy <- "dummy"
  return(mergedData)
}

# mergedData <- loadData()

pheno_RNA <- merge(x = mergedData, y = temp, by = c("Cohort","Subject","TimeCategory"))
univScatter(data=pheno_RNA, xData = "HALLMARK_IL2_STAT5_SIGNALING", yData="cTfh_ICOShiCD38hi_..FreqParent", 
            fillParam = "TimeCategory", title = "+/+ cTfh vs GSVA IL2/STAT5 at oneWeek", xLabel= "GSVA score - IL2-STAT5", yLabel = "ICOS+CD38+ cTfh frequency") + 
  coord_cartesian(ylim = c(0,18)) + scale_y_continuous(breaks = c(seq(0,18,2)))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSVA_IL2STAT5_vs_cTfh_oW_univ.pdf", device="pdf", height=8, width=8)
subsetData1 <- subset(pheno_RNA, Cohort == "Healthy")    ; subsetData2 <- subset(pheno_RNA, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "HALLMARK_IL2_STAT5_SIGNALING", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "+/+ cTfh vs GSVA IL2/STAT5 at oneWeek", xLabel= "GSVA score - IL2-STAT5", yLabel = "ICOS+CD38+ cTfh frequency") + 
  coord_cartesian(ylim = c(0,18)) + scale_y_continuous(breaks = c(seq(0,18,2)))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSVA_IL2STAT5_vs_cTfh_oW_biv.pdf", device="pdf", height=8, width=8)

univScatter(data=pheno_RNA, xData = "HALLMARK_MITOTIC_SPINDLE", yData="cTfh_ICOShiCD38hi_..FreqParent", 
            fillParam = "TimeCategory", title = "+/+ cTfh vs GSVA Mitotic Spindle at oneWeek", xLabel= "GSVA score - Mitotic Spindle", yLabel = "ICOS+CD38+ cTfh frequency") + 
  coord_cartesian(ylim = c(0,18)) + scale_y_continuous(breaks = c(seq(0,18,2)))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSVA_MitoticSpindle_vs_cTfh_oW_univ.pdf", device="pdf", height=8, width=8)
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "HALLMARK_MITOTIC_SPINDLE", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "+/+ cTfh vs GSVA Mitotic Spindle at oneWeek", xLabel= "GSVA score - Mitotic spindle", yLabel = "ICOS+CD38+ cTfh frequency") + 
  coord_cartesian(ylim = c(0,18)) + scale_y_continuous(breaks = c(seq(0,18,2)))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSVA_MitoticSpindle_vs_cTfh_oW_biv.pdf", device="pdf", height=8, width=8)


univScatter(data=pheno_RNA, xData = "HALLMARK_E2F_TARGETS", yData="cTfh_ICOShiCD38hi_..FreqParent", 
            fillParam = "TimeCategory", title = "+/+ cTfh vs GSVA E2F_TARGETS at oneWeek", xLabel= "GSVA score - E2F_TARGETS", yLabel = "ICOS+CD38+ cTfh frequency") + 
  coord_cartesian(ylim = c(0,18)) + scale_y_continuous(breaks = c(seq(0,18,2)))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSVA_E2F_TARGETS_vs_cTfh_oW_univ.pdf", device="pdf", height=8, width=8)
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "HALLMARK_E2F_TARGETS", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "+/+ cTfh vs GSVA E2F_TARGETS at oneWeek", xLabel= "GSVA score - E2F_TARGETS", yLabel = "ICOS+CD38+ cTfh frequency") + 
  coord_cartesian(ylim = c(0,18)) + scale_y_continuous(breaks = c(seq(0,18,2)))
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/GSVA_E2F_TARGETS_vs_cTfh_oW_biv.pdf", device="pdf", height=8, width=8)



# ----------- ----------- 

#' ## -                           link RNAseq to phenotypic analysis 

# ----------- ----------- 

temp <- colnames(logDataMatrix)
temp <- str_replace(temp, "^X","")
temp <- str_replace(temp, "[.]","-")
temp <- str_replace(temp, "_S..",""); temp <- str_replace(temp, "_S.","")
temp <- str_replace(temp, "_bL","_baseline"); temp <- str_replace(temp, "_oW","_oneWeek");

loadData <- function()   # run if mergedData not loaded yet
{
  mergedData <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/allMergedData.csv", stringsAsFactors = F, header = T, row.names = 1)
  mergedData$TimeCategory <- factor(mergedData$TimeCategory, levels = c("baseline", "oneWeek","late"))
  mergedData$Cohort <- factor(mergedData$Cohort, levels = c("Healthy", "aPD1","nonPD1"))
  mergedData$dummy <- "dummy"
  return(mergedData)
}                     # mergedData <- loadData()

temp2 <- logDataMatrix; colnames(temp2) <- temp; temp2 <- t(temp2)


melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_NaiveB...FreqParent"))

mergedData
pheno_RNA <- merge(x = mergedData, y = temp, by = c("Cohort","Subject","TimeCategory"))














