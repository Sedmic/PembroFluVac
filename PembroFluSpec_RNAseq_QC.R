library("gplots")
library("ggplot2")
library("ggrepel")
library("DESeq2")
library("viridis")
library("Rtsne")
sessionInfo()
source('D:/Pembro-Fluvac/Analysis/PembroFluSpec_Ranalysis_files/PembroFluSpec_PlottingFunctions.R')


dataMatrix <- read.table("D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/PORTresults/NORMALIZED_DATA/SPREADSHEETS/FINAL_master_list_of_gene_counts_MIN.Fastqs.txt", sep="", 
                         header=T, stringsAsFactors = F)
rownames(dataMatrix) <- make.names(dataMatrix$geneSymbol, unique=TRUE)
dataMatrix$id <- dataMatrix$geneCoordinate <- dataMatrix$geneSymbol <- NULL 

colSums2(as.matrix(dataMatrix))

temp <- dataMatrix +0.001
x <- log(apply(temp, 1, sum))
hist(x, breaks=100, xlab = "log2 counts", col = "lightblue",main = "Log2 counts distribution") 

sampleDists <- dist(t(dataMatrix))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- inferno(255) 
hc <- hclust(sampleDists)
# png(filename = "Images/AllSamples_bestTranscripts_ClusterDendrogram.png", width=1500, res = 75)
plot(hc, main="Cluster dendrogram of best transcripts")
# dev.off()
# pdf(file = "Images/AllSamples_SampleDistanceMatrix_heatmap.pdf")
heatmap.2(sampleDistMatrix, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=colors, margins=c(2,10), labCol=FALSE,cexRow = 0.55)
# dev.off()


pca1 = prcomp(t(dataMatrix), scale = FALSE)
# create data frame with scores
scores = as.data.frame(pca1$x)
scores$Cohort <- "a"; scores$TimeCategory <- "a"; scores$Subset <- "a"
scores$Cohort[grep("196",rownames(scores),value=F)] <- "aPD1" ;  scores$Cohort[grep("FS1819",rownames(scores),value=F)] <- "Healthy"
scores$TimeCategory[grep("oW",rownames(scores), value=F)] <- "oneWeek"; scores$TimeCategory[grep("bL",rownames(scores), value=F)] <- "baseline"; 
scores$Subset[grep("HiHi",rownames(scores),value=F)] <- "ICOShiCD38hi cTfh"; scores$Subset[grep("ABC",rownames(scores),value=F)] <- "Activated B cells";
scores$Subset[grep("_PB_",rownames(scores),value=F)] <- "Plasmablasts"; scores$Subset[grep("naiveB",rownames(scores),value=F)] <- "Naive B cells";
# plot of observations
ggplot(data = scores, aes(x = PC1, y = PC2, color = Cohort)) + theme_bw() +
  geom_hline(yintercept = 0, colour = "gray65") +   geom_vline(xintercept = 0, colour = "gray65") +
  geom_point( alpha=1, size=3) + 
  geom_text_repel(aes(label=Subset), alpha = 1, size = 3) +
  ggtitle("PCA plot for all samples in final pool") + theme(title = element_text(size=20))

ggplot(data = scores, aes(x = PC3, y = PC4, color = Cohort)) + theme_bw() +
  geom_hline(yintercept = 0, colour = "gray65") +   geom_vline(xintercept = 0, colour = "gray65") +
  geom_point( alpha=1, size=3) + 
  geom_text_repel(aes(label=Subset), alpha = 1, size = 3) +
  ggtitle("PCA plot for all samples in final pool") + theme(title = element_text(size=20))


## ***************************     Log transformation  **************************************

metaData <- data.frame(row.names=colnames(dataMatrix));  metaData$condition <- "empty"
metaData$Cohort <- metaData$Subject <- metaData$TimeCategory <- metaData$Subset <- "a"
metaData$Cohort[grep("196",rownames(metaData),value=F)] <- "aPD1" ;  metaData$Cohort[grep("FS1819",rownames(metaData),value=F)] <- "Healthy"
metaData$TimeCategory[grep("oW",rownames(metaData), value=F)] <- "oneWeek"; metaData$TimeCategory[grep("bL",rownames(metaData), value=F)] <- "baseline"; 
metaData$Subset[grep("HiHi",rownames(metaData),value=F)] <- "HiHi_cTfh"; metaData$Subset[grep("ABC",rownames(metaData),value=F)] <- "ABC";
metaData$Subset[grep("_PB_",rownames(metaData),value=F)] <- "PB"; metaData$Subset[grep("naiveB",rownames(metaData),value=F)] <- "navB";
metaData$Subject <- substr(rownames(metaData), start = 8, stop = 10)
metaData$condition <- paste0( metaData$Cohort,"_",metaData$Subset,"_",metaData$TimeCategory)

fullDataset <- DESeqDataSetFromMatrix(countData=dataMatrix, colData=metaData, design= ~ condition)

# vsd <- varianceStabilizingTransformation(fullDataset)
# logDataMatrix <- as.data.frame(assay(vsd))

# write.csv(logDataMatrix, file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/logDataMatrix.csv")
# logDataMatrix <- read.csv( file="D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/logDataMatrix.csv", stringsAsFactors = F, header = 1)


## ***************************     tSNE projection  **************************************

set.seed(42)
tsneMap <- Rtsne(t(logDataMatrix),epoch=50,perplexity=10,theta=0,verbosity=TRUE,max_iter=5000)
tsneMap <- as.data.frame(tsneMap$Y); tsneMap$subgroup <- metaData$condition; rownames(tsneMap) <- colnames(logDataMatrix)
tsneMap$Cohort <- metaData$Cohort; tsneMap$Subset <- metaData$Subset; tsneMap$TimeCategory <- metaData$TimeCategory
tsneMap$Subset <- factor(tsneMap$Subset, levels = c("ABC", "HiHi_cTfh","navB","PB" ))
tsneMap$TimeCategory <- factor(tsneMap$TimeCategory, levels = c("baseline","oneWeek" ))
tsneMap$Cohort <- factor(tsneMap$Cohort, levels = c("Healthy","aPD1" ))
customPalette <- colorRampPalette(brewer.pal(12,"Paired"))(12)

ggplot(tsneMap, aes(x=V1, y=V2)) + 
  geom_point(size=10,pch=21,colour="black", aes(fill= Subset )) +                                     # color by subset
  xlab("") + ylab("") +   ggtitle("t-SNE - By lymphocyte subset") +  theme_light(base_size=24)  + 
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort", legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) #+ xlim(-60,60) + ylim(-55,70)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/tsne_AllSamples_colorBySubset.pdf", device="pdf")

ggplot(tsneMap, aes(x=V1, y=V2)) + 
  geom_point(size=10,pch=21,colour="black", aes(fill=TimeCategory)) +                                      # color by TimeCategory
  xlab("") + ylab("") +   ggtitle("t-SNE - By time after vaccine") +  theme_light(base_size=24)  +  
  theme(strip.background = element_blank()) + labs(fill = "Time Category", legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) #+ xlim(-60,60) + ylim(-55,70)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/tsne_AllSamples_colorByTimeCategory.pdf", device="pdf")


ggplot(tsneMap, aes(x=V1, y=V2)) + 
  geom_point(size=10,pch=21,colour="black", aes(fill=Cohort )) +                                      # color by Cohort
  xlab("") + ylab("") +   ggtitle("t-SNE - By cohort") +  theme_light(base_size=24)  + 
  theme(strip.background = element_blank()) + labs(fill = "Sample cohort", legend.key = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=8))) #+ xlim(-60,60) + ylim(-55,70)
# ggsave(filename = "D:/Pembro-Fluvac/18-19season/RNAseq/Analysis/Images/tsne_AllSamples_colorByCohort.pdf", device="pdf")



# unique(metaData$condition[order(metaData$condition)])       "aPD1_ABC_baseline", "aPD1_ABC_oneWeek",  "aPD1_HiHi_cTfh_baseline", "aPD1_HiHi_cTfh_oneWeek", "aPD1_navB_baseline", "aPD1_navB_oneWeek", "aPD1_PB_baseline",  "aPD1_PB_oneWeek", "Healthy_ABC_baseline", "Healthy_ABC_oneWeek", "Healthy_HiHi_cTfh_baseline", "Healthy_HiHi_cTfh_oneWeek","Healthy_navB_baseline", "Healthy_navB_oneWeek",  "Healthy_PB_baseline", "Healthy_PB_oneWeek"   
