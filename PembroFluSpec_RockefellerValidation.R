#+ fig.width=8, fig.height=8
#' ##--------- call libraries --
library("genefilter")
library("ggplot2")
library("grid");  library("reshape2")
library("ggrepel")
library("ggbeeswarm")
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
source('D:/Pembro-Fluvac/Analysis/Ranalysis/PembroFluSpec_PlottingFunctions.R')
sessionInfo()



#' ## ----------- BAA Year 3 analysis  --------------------
#'

BAAyear3long <- read.csv(file="D:/Pembro-Fluvac/Analysis/mergedData/BAAyear3long.csv"); names(BAAyear3long)[6] <- 'TimePoint'
BAAyear3long$Cohort <- factor(BAAyear3long$Cohort, levels=c("Young","Elderly"))
BAAyear3wide <- read.csv(file="D:/Pembro-Fluvac/Analysis/mergedData/BAAyear3wide.csv") 
BAAyear3wide$Cohort <- factor(BAAyear3wide$Cohort, levels=c("Young","Elderly"))
prePostTime(data = BAAyear3long, xData = 'TimePoint', yData = 'ICOShiCD38hi_FreqParent', fillParam = 'Cohort', groupby = 'Subject', title = "cTfh with Aging",
            xLabel = 'Age Group', yLabel = "ICOS+CD38+ (% cTfh)") + scale_y_continuous(breaks=seq(0,60,5),limits = c(0,55)) + 
  scale_fill_manual(values = c("#EABD81", "#7B4D9A"))  + ylab(expression(paste(paste("ICO",S^'+', "CD3",8^'+' ), " (% cTfh)")))

# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/BAAyear3_FCtfh.pdf")
tukey_hsd( aov(data = BAAyear3long, formula = `ICOShiCD38hi_FreqParent` ~ `Cohort`:`TimePoint`))
twoSampleBar(data = BAAyear3wide, xData = "Cohort", yData = "FCtfh", fillParam = "Cohort", title = "ICOS+CD38+ cTfh", yLabel = "Fold-change", FCplot = T) + 
  scale_y_continuous(breaks=seq(0,6,1), limits=c(0,6))


#' ## ----------- Dataset assembly  --------------------
#'
#'
#'
#'
setwd("D:/Pembro-Fluvac/Analysis")

SubjTable.Rock <- read.csv(file="D:/Pembro-Fluvac/Analysis/mergedData/Rockefeller_SubjectTable.csv")
mergedDataRock <- SubjTable.Rock
 
serology.Rock <- read.csv(file= "D:/Pembro-Fluvac/Analysis/mergedData/Rockefeller_Serology.csv", header = T)
print("Data by year and cohort for serology")
temp1 <- merge(x = mergedDataRock, y=serology.Rock, all = T, suffixes = c(".Tflow",".Serology"), by= c('Subject'))

TflowRock <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/Rockefeller_Tflow_freqParent.csv", header = T)
print("Data by year and cohort for T cell flow cytometry for Rockefeller Cohort")


temp2 <- merge(x = temp1, y=TflowRock, all=T, suffixes = c(".temp1",".TflowRock"), by=c('Label'))

TmfiRock <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/Rockefeller_Tmfi.csv")
temp3 <- merge(x = temp2, y=TmfiRock, all=T, suffixes = c(".temp2",".TmfiRock"), by=c('Label'))
mergedDataRock <- temp3  

mergedDataRock <- mergedDataRock[-which(mergedDataRock$Subject == "13-226-91"),]    # ambiguity about what time points are being studied
mergedDataRock <- mergedDataRock[-which(mergedDataRock$Subject == "13-226-67"),]    # ambiguity about what time points are being studied
mergedDataRock <- mergedDataRock[-which(mergedDataRock$Label == "13-226-51E"),]     # two draws at oneWeek, this one is at d10

table(mergedDataRock$Cohort, mergedDataRock$TimeCategory)

mergedDataRock$logHAI <- log(mergedDataRock$H1N1pdm09.HAI.titer)

#' ## ----------- calculate fold-changes for HAI, Tfh, and PB --------------------
subsetData <- subset(mergedDataRock, TimeCategory != "oneWeek" &  Cohort!="")
FC_response <- dcast(subsetData, `Subject` + `Cohort` ~`TimeCategory`, value.var = c("logHAI"))
FC_response$FChai_late <- FC_response$`late`/FC_response$`baseline`; FC_response$baseline <- FC_response$late <- NULL

FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("H1N1pdm09.HAI.titer")) 
FC_response2$FCH1N1_hai_late <- FC_response2$`late`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject",all.x=T); FC_response$baseline <- FC_response$late <- NULL
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("H3N2.HAI.titer")) 
FC_response2$FCH3N2_hai_late <- FC_response2$`late`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject",all.x=T); FC_response$baseline <- FC_response$late <- NULL
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("FluB.HAI.titer")) 
FC_response2$FCFluB_hai_late <- FC_response2$`late`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject", all.x=T); FC_response$baseline <- FC_response$late <- NULL


subsetData <- subset(mergedDataRock, TimeCategory != "late" &  Cohort!="")
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("ICOS.CD38..cTfh.freq")) 
FC_response2$FCtfh_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject",all.x=T); FC_response$baseline <- FC_response$oneWeek <- NULL
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD19..CD27.CD38..freq")) 
FC_response2$FCPB_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject", all.x=T); FC_response$baseline <- FC_response$oneWeek <- NULL
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("Plasma.CXCL13..pg.mL.")) 
FC_response2$FCCXCL13_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject", all.x=T); FC_response$baseline <- FC_response$oneWeek <- NULL
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgG1_Total.sialylated")) 
FC_response2$IgG1sial_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject", all.x=T); FC_response$baseline <- FC_response$oneWeek <- NULL

mergedDataRock <- merge(x = mergedDataRock, y= FC_response, all=T, by = c('Subject', 'Cohort'))


# write.csv(mergedDataRock, file = "D:/Pembro-Fluvac/Analysis/mergedData/allmergedDataRock.csv")
# mergedDataRock <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/allmergedDataRock.csv", stringsAsFactors = F, header = T, row.names = 1)

mergedDataRock$dummy <- "dummy"
mergedDataRockComplete <-  mergedDataRock

mergedDataRock <- mergedDataRockComplete
mergedDataRock <- mergedDataRock[ -grep(paste(c("IpiNivo"), collapse = "|"), mergedDataRock$Cohort), ]   # leave out aPD1 + TKI and IpiNivo
mergedDataRock$Cohort <- factor(mergedDataRock$Cohort, levels = c("Chemo","TKI","aPDL1","aPD1","aPD1_TKI"))          
mergedDataRock$TimeCategory <- factor(mergedDataRock$TimeCategory, levels = c("baseline", "oneWeek","late"))

saveRDS(mergedDataRock, file="D:/Pembro-Fluvac/Analysis/mergedData/mergedDataRock.RDS")
mergedDataRock <- readRDS(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergedDataRock.RDS")

#' ## ----------- demographics --------------------
#'
#'
temp <- mergedDataRock[,1:10]
temp$Cohort <- ifelse(temp$Cohort == 'aPD1' | temp$Cohort == 'aPD1_TKI', 'aPD1', 'non-aPD1')
temp <- temp[ !duplicated(temp$Subject), ]    # eliminate duplicate rows
table(temp$Cohort)
table(temp$Cohort, temp$Treatment) 
temp[which(temp$Cohort == "aPD1" ),] %>% get_summary_stats(type="common")        
temp[which(temp$Cohort == "non-aPD1" ),] %>% get_summary_stats(type="common")        
table(temp[which(temp$Cohort == "aPD1"), "Sex"]) %>% prop.table()
table(temp[which(temp$Cohort == "non-aPD1"), "Sex"]) %>% prop.table()

ggplot( data = subset(temp, Cohort == "aPD1"), aes(x=Cycle)) + geom_histogram(color="white",fill="#FFB18C") + ggtitle("Cohort 1 - aPD1 cohort") + 
  theme_bw() + theme(axis.text=element_text(color="black",size=18), axis.title=element_text(color="black",size=24), plot.title = element_text(size=24)) + 
  scale_x_continuous(breaks=seq(0,100,5)) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CyclesImmunotherapy_ROCK.pdf", device = 'pdf', height=3, width=8)




#' ## ----------------- quick look at dataset ----------------
#' 

temp <- mergedDataRockComplete
temp$Cohort <- factor(temp$Cohort, levels = c("Chemo","TKI","aPDL1","aPD1","aPD1_TKI","IpiNivo"))
data=subset(temp, Cohort != "nonPD1" & TimeCategory == "oneWeek"); xData = "Cohort"; yData = "FCtfh_oW"; fillParam = 'Cohort'; title = "ICOS+CD38+ cTfh responses"; 
yLabel = "Fold-change at one week";
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T); names(overTime)[which(names(overTime) == 'x')] <- yData
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + #scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) + 
  geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity") + 
  geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.15)) + 
  ggtitle(title) + ylab(yLabel) +  theme_bw() +
  theme(axis.text = element_text(size=28,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=32,hjust = 0.5)) 

summary(fit<-aov(FCtfh_oW ~ Cohort, data=data)); tukey_hsd(fit)

# group & factor by Cohort but will display symbols by treatment
mergedDataRock$Cohort <- ifelse(mergedDataRock$Cohort == 'aPD1' | mergedDataRock$Cohort == 'aPD1_TKI', 'aPD1', 'non-aPD1')
mergedDataRock$Cohort <- factor(mergedDataRock$Cohort, levels = c("non-aPD1","aPD1"))

mergedDataRock$TreatmentFactored <- mergedDataRock$Treatment
mergedDataRock$TreatmentFactored[grep(paste0(c("Axitinib", "Pazopanib", "Cabozantinib" ), collapse = "|"),mergedDataRock$TreatmentFactored)] <- "TKI"
mergedDataRock$TreatmentFactored[grep(paste0(c("Nivolumab$", "Pembrolizumab$"), collapse = "|"),mergedDataRock$TreatmentFactored)] <- "aPD1"
mergedDataRock$TreatmentFactored[grep(paste0(c("Nivolumab/cabozantinib", "Pembrolizumab/lenvatinib" ), collapse = "|"),mergedDataRock$TreatmentFactored)] <- "aPD1+TKI"
mergedDataRock$TreatmentFactored[grep(paste0(c("Atezolizumab" ), collapse = "|"),mergedDataRock$TreatmentFactored)] <- "aPDL1"
mergedDataRock$TreatmentFactored[grep(paste0(c("Carboplatin/Gemcitabine" ), collapse = "|"),mergedDataRock$TreatmentFactored)] <- "Chemo"
mergedDataRock$TreatmentFactored <- factor(mergedDataRock$TreatmentFactored, levels = c("Chemo", "TKI","aPDL1","aPD1","aPD1+TKI"))
mergedDataRock.nolate <- mergedDataRock[-which(mergedDataRock$TimeCategory == 'late'),]
mergedDataRock.nolate$TimeCategory <- factor(mergedDataRock.nolate$TimeCategory, levels = c("baseline","oneWeek"))



##' ## ----------- Tfh analysis --------------------
#'
#'

subsetData <- subset(mergedDataRock, TimeCategory == 'oneWeek')
data=subsetData; xData="Cohort"; yData="FCtfh_oW"; fillParam="Cohort"; title= expression(paste(paste("ICO",S^'+', "CD3",8^'+' ), " cTfh")); 
yLabel="Fold-change at one week"; position="left"
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T); names(overTime)[which(names(overTime) == 'x')] <- yData
justforttest <- data[, c(xData,yData)]; fit <- rstatix::wilcox_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p;   annotationInfo <- paste0("P = ", round(pValue, 2))
my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.8)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))  + 
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5) + 
  geom_bar(data=overTime, aes_string(x=xData, y=yData, width=0.6), position = position_dodge(), stat = "identity",color="black",size=0.1) + 
  #geom_point(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5, position = position_jitter(width=0.3, seed=58)) +   #48
  geom_quasirandom(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5) + 
  ggtitle(title) + ylab(yLabel) +  theme_bw() + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=28,hjust = 0.5),axis.text.x=element_text(angle=45,vjust=1,hjust=1), legend.text = element_text(size=14), legend.title = element_text(size=14)) + annotation_custom(my_grob)  + 
  scale_y_continuous(breaks = seq(0,12,1), limits=c(0,11))  + theme(axis.title.y = element_text(vjust=-0.5))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResponse_foldchange_ROCK.pdf", width = 5)
wilcox_test(FCtfh_oW ~ Cohort, data=subsetData)





# prePostTime(data = mergedDataRock.nolate, xData = "TimeCategory", yData = "ICOS.CD38..cTfh.freq", fillParam = "Cohort", title = "cTfh response by subject", xLabel = "TimeCategory",
            # yLabel = "ICOS+CD38+ (% cTfh)", groupby = "Subject")  + scale_y_continuous(breaks=seq(0,150,1), limits=c(0,18)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_PerSubject_ROCK.pdf", width = 6)



##' ## ----------- B cell analysis --------------------
#'
#'

subsetData <- subset(mergedDataRock, TimeCategory == 'oneWeek')
subsetData[which(subsetData$Cohort == "non-aPD1"),"FCPB_oW"]
outliers::grubbs.test(subsetData[which(subsetData$Cohort == "non-aPD1"),"FCPB_oW"])
outlier <- which(subsetData$FCPB_oW == max(subsetData[which(subsetData$Cohort == "non-aPD1"),"FCPB_oW"], na.rm = T))
subsetData <- subsetData[-outlier, ]
data=subsetData; xData="Cohort"; yData="FCPB_oW"; fillParam="Cohort"; title="Plasmablasts"; yLabel="Fold-change at one week"; position="left"
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T); names(overTime)[which(names(overTime) == 'x')] <- yData
justforttest <- data[, c(xData,yData)]; fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p;   annotationInfo <- paste0("P = ", round(pValue, 2))
my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28))) 
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.8)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))  + 
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5) + 
  geom_bar(data=overTime, aes_string(x=xData, y=yData, width=0.65), position = position_dodge(), stat = "identity",color="black",size=0.1) + 
  # geom_point(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5, position = position_jitter(width=0.3, seed=59)) +   #48
  geom_quasirandom(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5) +   
  ggtitle(title) + ylab(yLabel) +  theme_bw() + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=28,hjust = 0.5),axis.text.x=element_text(angle=45,vjust=1,hjust=1), legend.text = element_text(size=14), legend.title = element_text(size=14)) + annotation_custom(my_grob) # + 
  # scale_y_continuous(breaks = seq(0,10,1))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/PBResponse_foldchange_ROCK.pdf", width = 5)
wilcox_test(FCPB_oW ~ Cohort, data=subsetData)



oneWeek <- subset(mergedDataRock, TimeCategory == "oneWeek" )   
subsetData1 <- subset(oneWeek, Cohort == "non-aPD1")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
data1 = subsetData1; data2 = subsetData2; name1 = "non-aPD1"; name2 = "aPD1"; yData = "CD19..CD27.CD38..freq"; xData="ICOS.CD38..cTfh.freq"; 
fillParam = "Cohort"; title = "cTfh vs PB at oneWeek"; yLabel= "CD27+CD38+ PB"; xLabel = "ICOS+CD38+ cTfh"
kendall1 <- round(cor(data1[,xData], data1[,yData], method = "kendall", use = "complete.obs"), 2)
pValue1 <- cor.test(data1[,xData], data1[,yData], method="kendall")
kendall2 <- round(cor(data2[,xData], data2[,yData], method = "kendall", use = "complete.obs"), 2)
pValue2 <- cor.test(data2[,xData], data2[,yData], method="kendall")
annotationInfo1 <- paste0(name1, " Kendall t = ", kendall1,"\n","P = ", formatC(pValue1$p.value, format="e", digits=1) )
annotationInfo2 <- paste0("\n", name2, " Kendall t = ", kendall2,"\n","P = ", formatC(pValue2$p.value, format="e", digits=1) )
my_grob1 = grobTree(textGrob(annotationInfo1, x=0.05,  y=0.88, hjust=0, gp=gpar(col="#7FAEDB", fontsize=28)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.75, hjust=0, gp=gpar(col="#ff9a6a", fontsize=28)))
ggplot() + 
  geom_point(data = data1, aes_string(x=xData, y=yData, fill= fillParam), size=8, color="black", pch=21,alpha=0.75) + theme_bw() + 
  geom_point(data = data2, aes_string(x=xData, y=yData, fill=fillParam), size=8, color="black", pch=21,alpha=0.75) + theme_bw() + 
  geom_smooth(data = data1, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
  geom_smooth(data = data2, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
  scale_color_manual(values=c("#ff9a6a", "#d9eafb")) + 
  scale_fill_manual(values=c("#ff9a6a", "#d9eafb")) + 
  ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
  theme(axis.text = element_text(size=28,hjust = 0.5,color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5)) + 
  annotation_custom(my_grob1) + annotation_custom(my_grob2) + theme(legend.position = "none") #+ scale_y_continuous(limits=c(0,50)) + scale_x_continuous(breaks=seq(0,15,1))


oneWeek <- subset(mergedDataRock, TimeCategory == "oneWeek" )   
subsetData1 <- subset(oneWeek, Cohort == "non-aPD1")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
data1 = subsetData1; data2 = subsetData2; name1 = "non-aPD1"; name2 = "aPD1"; yData = "FCPB_oW"; xData="FCtfh_oW"; 
fillParam = "Cohort"; title = "Cohort 1 - One Week"; yLabel= expression(paste("Fold-change in Plasmablasts")); xLabel = expression(paste("Fold-change in ICO",S^'+',"CD3",8^'+'," cTfh")) 
kendall1 <- round(cor(data1[,xData], data1[,yData], method = "kendall", use = "complete.obs"), 2)
pValue1 <- cor.test(data1[,xData], data1[,yData], method="kendall")
kendall2 <- round(cor(data2[,xData], data2[,yData], method = "kendall", use = "complete.obs"), 2)
pValue2 <- cor.test(data2[,xData], data2[,yData], method="kendall")
annotationInfo1 <- paste0(name1, " Kendall t = ", kendall1,"\n","P = ", formatC(pValue1$p.value, format="e", digits=1) )
annotationInfo2 <- paste0("\n", name2, " Kendall t = ", kendall2,"\n","P = ", formatC(pValue2$p.value, format="e", digits=1) )
my_grob1 = grobTree(textGrob(annotationInfo1, x=0.05,  y=0.88, hjust=0, gp=gpar(col="#7FAEDB", fontsize=28)))
my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.75, hjust=0, gp=gpar(col="#ff9a6a", fontsize=28)))
ggplot() + 
  geom_point(data = data1, aes_string(x=xData, y=yData, fill= fillParam), size=8, color="black", pch=21,alpha=0.75) + theme_bw() + 
  geom_point(data = data2, aes_string(x=xData, y=yData, fill=fillParam), size=8, color="black", pch=21,alpha=0.75) + theme_bw() + 
  geom_smooth(data = data1, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
  geom_smooth(data = data2, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
  scale_color_manual(values=c("#ff9a6a", "#d9eafb")) + 
  scale_fill_manual(values=c("#ff9a6a", "#d9eafb")) + 
  ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
  theme(axis.text = element_text(size=28,hjust = 0.5,color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5)) + 
  # annotation_custom(my_grob1) + annotation_custom(my_grob2) + 
  theme(legend.position = "none") + scale_y_continuous(limits=c(0,50)) + scale_x_continuous(breaks=seq(0,15,1))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/PB_correl_cTfh_FC_ROCK.pdf", width =8, height=6)



##' ## ----------- CXCL13 analysis --------------------
#'
#'

subsetData <- subset(mergedDataRock, TimeCategory == 'oneWeek')
data=subsetData; xData="Cohort"; yData="FCCXCL13_oW"; fillParam="Cohort"; title="Plasma CXCL13"; yLabel="Fold-change at one week"; position="left"
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T); names(overTime)[which(names(overTime) == 'x')] <- yData
justforttest <- data[, c(xData,yData)]; fit <- rstatix::wilcox_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p;   annotationInfo <- paste0("P = ", round(pValue, 3))
my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.8)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))  + 
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5) + 
  geom_bar(data=overTime, aes_string(x=xData, y=yData, width=0.65), position = position_dodge(), stat = "identity",color="black",size=0.1) + 
  # geom_point(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5, position = position_jitter(width=0.3, seed=58)) +   #48
  geom_quasirandom(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5) + 
  ggtitle(title) + ylab(yLabel) +  theme_bw() + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=28,hjust = 0.5),axis.text.x=element_text(angle=45,vjust=1,hjust=1), legend.text = element_text(size=14), legend.title = element_text(size=14)) + annotation_custom(my_grob)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_foldchange_ROCK.pdf", width = 4.5)
wilcox_test(FCCXCL13_oW ~ Cohort, data=subsetData)

subsetData <- subset(mergedDataRock, TimeCategory != 'late')
data=subsetData; xData="TimeCategory"; yData="Plasma.CXCL13..pg.mL."; fillParam="Cohort"; title="CXCL13 by subject"; yLabel="Plasma CXCL13 (pg/mL)"; xLabel="TimeCategory";
groupby = "Subject"

subsetData <- subsetData[order(subsetData$Subject, subsetData$TimePoint, decreasing = F),]
ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
    geom_bar(stat = "summary", aes_string(fill=fillParam), color="black",size=0.1) + 
    geom_path(aes_string(group=groupby), color="grey40", alpha=0.5) + 
    geom_point(size = 3, color="black", fill="black", alpha=0.5, aes(shape=TreatmentFactored)) + facet_wrap(fillParam ) +   
    scale_color_manual(values=c("#d9eafb", "#ff9a6a")) + scale_fill_manual(values=c("#d9eafb", "#ff9a6a")) + 
    ggtitle(title) + ylab(yLabel) + xlab(xLabel)  + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
    theme(axis.text = element_text(size=18,hjust = 0.5, color="black"), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
          strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white"), # legend.position = "none",
          axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1,vjust=1), legend.text = element_text(size=14), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_overTime_PerSubject_ROCK.pdf", width = 6)

  


melted <- melt(mergedDataRock.nolate, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("Plasma.CXCL13..pg.mL."))
prePostTimeAveraged(melted, title = "CXCL13 responses", xLabel = NULL, yLabel = "Plasma CXCL13 (pg/mL)") + scale_y_continuous(breaks=seq(0,100,10), limits=c(0,100)) + 
  scale_fill_manual(values = c( "#abceef", "#ff9a6a"))  
summary( fit <- aov(value ~ Cohort+TimeCategory + Cohort:TimeCategory, data=melted ) );  tukey_hsd(fit)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_TimeAveraged_ROCK.pdf", device="pdf", width=5, height=7)

subsetData <- subset(melted, TimeCategory == "oneWeek")
rownames(subsetData) <- seq(1:nrow(subsetData))
twoSampleBar(data=subsetData, xData="Cohort", yData="value", fillParam="Cohort", title="CXCL13\nOne week", yLabel="plasma CXCL13 (pg/mL)") + 
  scale_y_continuous(limits = c(0,250)) + scale_fill_manual(values = c( "#abceef", "#ff9a6a"))  
aggregate( subsetData, by= list(subsetData$variable, subsetData$TimeCategory, subsetData$Cohort), FUN=mean, na.rm = T)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_oW_ROCK.pdf", device="pdf",  width = 4)



##' ## ----------- NeutAb - HAI analysis --------------------
#'
#'


data=mergedDataRock; xData="TimeCategory"; yData="H1N1pdm09.HAI.titer"; fillParam="Cohort"; title="H1N1"; yLabel="HAI titer"; xLabel="TimeCategory"; groupby = "Subject"
subsetData <- subsetData[order(subsetData$Subject, subsetData$TimePoint, decreasing = F),]
ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
  geom_bar(stat = "summary", aes_string(fill=fillParam), color="black",size=0.1) + 
  geom_path(aes_string(group=groupby), color="grey40", alpha=0.5) + 
  geom_point(size = 3, color="black", fill="black", alpha=0.5, aes(shape=TreatmentFactored)) + facet_wrap(fillParam ) +   
  scale_color_manual(values=c("#d9eafb", "#ff9a6a")) + scale_fill_manual(values=c("#d9eafb", "#ff9a6a")) + 
  ggtitle(title) + ylab(yLabel) + xlab(xLabel)  + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=22,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
        strip.text = element_text(size = 28, color="black"), strip.background = element_rect(fill="white"), # legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1,vjust=1), legend.text = element_text(size=18), legend.title = element_text(size=18))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAItiter_overTime_PerSubject_H1N1_ROCK.pdf", width=7)

data=mergedDataRock; xData="TimeCategory"; yData="H3N2.HAI.titer"; fillParam="Cohort"; title="H3N2"; yLabel="HAI titer"; xLabel="TimeCategory"; groupby = "Subject"
subsetData <- subsetData[order(subsetData$Subject, subsetData$TimePoint, decreasing = F),]
ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
  geom_bar(stat = "summary", aes_string(fill=fillParam), color="black",size=0.1) + 
  geom_path(aes_string(group=groupby), color="grey40", alpha=0.5) + 
  geom_point(size = 3, color="black", fill="black", alpha=0.5, aes(shape=TreatmentFactored)) + facet_wrap(fillParam ) +   
  scale_color_manual(values=c("#d9eafb", "#ff9a6a")) + scale_fill_manual(values=c("#d9eafb", "#ff9a6a")) + 
  ggtitle(title) + ylab(yLabel) + xlab(xLabel)  + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=22,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
        strip.text = element_text(size = 28, color="black"), strip.background = element_rect(fill="white"), # legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1,vjust=1), legend.text = element_text(size=18), legend.title = element_text(size=18))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAItiter_overTime_PerSubject_H3N2_ROCK.pdf", width=7)

data=mergedDataRock; xData="TimeCategory"; yData="FluB.HAI.titer"; fillParam="Cohort"; title="FluB"; yLabel="HAI titer"; xLabel="TimeCategory"; groupby = "Subject"
subsetData <- subsetData[order(subsetData$Subject, subsetData$TimePoint, decreasing = F),]
ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
  geom_bar(stat = "summary", aes_string(fill=fillParam), color="black",size=0.1) + 
  geom_path(aes_string(group=groupby), color="grey40", alpha=0.5) + 
  geom_point(size = 3, color="black", fill="black", alpha=0.5, aes(shape=TreatmentFactored)) + facet_wrap(fillParam ) +   
  scale_color_manual(values=c("#d9eafb", "#ff9a6a")) + scale_fill_manual(values=c("#d9eafb", "#ff9a6a")) + 
  ggtitle(title) + ylab(yLabel) + xlab(xLabel)  + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=22,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
        strip.text = element_text(size = 28, color="black"), strip.background = element_rect(fill="white"), # legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1,vjust=1), legend.text = element_text(size=18), legend.title = element_text(size=18))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAItiter_overTime_PerSubject_FluB_ROCK.pdf", width=7)



seroconv <- dcast(mergedDataRock, `Subject`+`Cohort`+`TreatmentFactored`~`TimeCategory`, value.var = c("H1N1pdm09.HAI.titer")); 
seroconv$FC <- seroconv$`late`/seroconv$`baseline`
data=seroconv; xData="Cohort"; yData="FC"; fillParam="Cohort"; title="H1N1"; yLabel="Fold-change in HAI"; 
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T); names(overTime)[which(names(overTime) == 'x')] <- yData
justforttest <- data[, c(xData,yData)]; fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p;   annotationInfo <- paste0("P = ", round(pValue, 2))
my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.8)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))  + 
  geom_bar(data=overTime, aes_string(x=xData, y=yData, width=0.65),  position = position_dodge(), stat = "identity",color="black",size=0.1) + 
  geom_quasirandom(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5) +
  ggtitle(title) + ylab(yLabel) +  theme_bw() + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=28,hjust = 0.5),axis.text.x=element_text(angle=45,vjust=1,hjust=1), legend.text = element_text(size=18), legend.title = element_text(size=18)) + annotation_custom(my_grob)  + 
  scale_y_continuous(trans='log2',breaks=c(2^(0:14)), limits = c(1,175))  
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/SeroconversionFactor_H1N1_ROCK.pdf", device = "pdf", width=5)


seroconv <- dcast(mergedDataRock, `Subject`+`Cohort`+`TreatmentFactored`~`TimeCategory`, value.var = c("H3N2.HAI.titer")); 
seroconv$FC <- seroconv$`late`/seroconv$`baseline`
data=seroconv; xData="Cohort"; yData="FC"; fillParam="Cohort"; title="H3N2"; yLabel="Fold-change in HAI"; 
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T); names(overTime)[which(names(overTime) == 'x')] <- yData
justforttest <- data[, c(xData,yData)]; fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p;   annotationInfo <- paste0("P = ", round(pValue, 2))
my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.8)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))  + 
  geom_bar(data=overTime, aes_string(x=xData, y=yData, width=0.65),  position = position_dodge(), stat = "identity",color="black",size=0.1) + 
  geom_quasirandom(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5) +
  ggtitle(title) + ylab(yLabel) +  theme_bw() + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=28,hjust = 0.5),axis.text.x=element_text(angle=45,vjust=1,hjust=1), legend.text = element_text(size=18), legend.title = element_text(size=18)) + annotation_custom(my_grob)  + 
  scale_y_continuous(trans='log2',breaks=c(2^(0:14)), limits = c(1,175))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/SeroconversionFactor_H3N2_ROCK.pdf", device = "pdf", width=5)

seroconv <- dcast(mergedDataRock, `Subject`+`Cohort`+`TreatmentFactored`~`TimeCategory`, value.var = c("FluB.HAI.titer")); 
seroconv$FC <- seroconv$`late`/seroconv$`baseline`
data=seroconv; xData="Cohort"; yData="FC"; fillParam="Cohort"; title="FluB"; yLabel="Fold-change in HAI"; 
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T); names(overTime)[which(names(overTime) == 'x')] <- yData
justforttest <- data[, c(xData,yData)]; fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p;   annotationInfo <- paste0("P = ", round(pValue, 2))
my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.8)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))  + 
  geom_bar(data=overTime, aes_string(x=xData, y=yData, width=0.65),  position = position_dodge(), stat = "identity",color="black",size=0.1) + 
  geom_quasirandom(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5) +
  ggtitle(title) + ylab(yLabel) +  theme_bw() + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=28,hjust = 0.5),axis.text.x=element_text(angle=45,vjust=1,hjust=1), legend.text = element_text(size=18), legend.title = element_text(size=18)) + annotation_custom(my_grob)  + 
  scale_y_continuous(trans='log2',breaks=c(2^(0:14)), limits = c(1,175))  
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/SeroconversionFactor_FluB_ROCK.pdf", device = "pdf", width=5)

##' ## ----------- Glycosylation analysis --------------------
#'
#'
subsetData <- subset(mergedDataRock )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgG1_Total.Galactosylation..G1.G2."))
prePostTimeAveraged(melted, title = "Galactosylation", xLabel = NULL, yLabel = "% anti-H1 IgG1") + scale_y_continuous(breaks=seq(0,100,1), limits = c(90,101))  +
  scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_galactosylation_overTime_ROCK.pdf", width=4)

melted %>%  anova_test(value ~ TimeCategory+Cohort )                        # two-way anova without interaction term 
Anova(fit <- lm( value ~ Cohort * TimeCategory, data = melted)); tukey_hsd(fit)


subsetData <- mergedDataRock # subset(mergedDataRock, TimeCategory != 'late')
data=subsetData; xData="TimeCategory"; yData="IgG1_Total.Galactosylation..G1.G2."; fillParam="Cohort"; title="Total galactosylation"; yLabel="% anti-H1 IgG1"; xLabel="TimeCategory";
groupby = "Subject"
subsetData <- subsetData[order(subsetData$Subject, subsetData$TimePoint, decreasing = F),]
ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
  geom_bar(stat = "summary", aes_string(fill=fillParam),color="black",size=0.1) + 
  geom_path(aes_string(group=groupby), color="grey60", alpha=0.5) + 
  geom_point(size = 3, color="black", fill="black",alpha=0.5, aes(shape=TreatmentFactored)) + facet_wrap(fillParam ) +   
  scale_color_manual(values=c("#d9eafb", "#ff9a6a")) + scale_fill_manual(values=c("#d9eafb", "#ff9a6a")) + 
  ggtitle(title) + ylab(yLabel) + xlab(xLabel)  + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=18,hjust = 0.5, color="black"), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
        strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white"), # legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1,vjust=1), legend.text = element_text(size=14), legend.title = element_text(size=14)) + coord_cartesian(ylim=c(80,100))

subsetData <- subset(mergedDataRock, TimeCategory == 'baseline')
shapiro_test(data = subsetData, vars = c("IgG1_Total.Galactosylation..G1.G2."))          # not normal distribution so will use nonparametric
data=subsetData; xData="Cohort"; yData="IgG1_Total.Galactosylation..G1.G2."; fillParam="Cohort"; title="Total galactosylation\nbaseline"; yLabel="% anti-H1 IgG1"; position="left"
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T); names(overTime)[which(names(overTime) == 'x')] <- yData
justforttest <- data[, c(xData,yData)]; fit <- rstatix::wilcox_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p;   annotationInfo <- paste0("P = ", round(pValue, 3))
my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=24))) 
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.8)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))  + 
  geom_bar(data=overTime, aes_string(x=xData, y=yData, width=0.65),  position = position_dodge(), stat = "identity",color="black",size=0.1) + 
  geom_quasirandom(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5) + 
  ggtitle(title) + ylab(yLabel) +  theme_bw() + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=28,hjust = 0.5),axis.text.x=element_text(angle=45,vjust=1,hjust=1), legend.text = element_text(size=18), legend.title = element_text(size=18)) + annotation_custom(my_grob)  + 
  scale_y_continuous(breaks = seq(0,100,1)) +   coord_cartesian(ylim = c(88,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_totGalactosylation_baseline_ROCK.pdf", width = 5)
wilcox_test(IgG1_Total.sialylated ~ Cohort, data=subsetData)



subsetData <- subset(mergedDataRock )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgG1_Total.sialylated"))
prePostTimeAveraged(melted, title = "Sialylation", xLabel = NULL, yLabel = "% anti-H1 IgG1") + scale_y_continuous(breaks=seq(0,99,2), limits = c(6,12)) + 
  scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_sialylation_overTime_ROCK.pdf", width=4)

melted %>%  anova_test(value ~ TimeCategory+Cohort )                        # two-way anova without interaction term 
Anova(fit <- lm( value ~ Cohort * TimeCategory, data = melted)); tukey_hsd(fit)


subsetData <- subset(mergedDataRock, TimeCategory != 'late')
data=subsetData; xData="TimeCategory"; yData="IgG1_Total.sialylated"; fillParam="Cohort"; title="Sialylation"; yLabel="% anti-H1 IgG1"; xLabel="TimeCategory";groupby = "Subject"
subsetData <- subsetData[order(subsetData$Subject, subsetData$TimePoint, decreasing = F),]
ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
  geom_bar(stat = "summary", aes_string(fill=fillParam),color="black",size=0.1) + 
  geom_path(aes_string(group=groupby), color="grey60") + 
  geom_point(size = 3, color="grey60", aes(shape=TreatmentFactored)) + facet_wrap(fillParam ) +   
  scale_color_manual(values=c("#d9eafb", "#ff9a6a")) + scale_fill_manual(values=c("#d9eafb", "#ff9a6a")) + 
  ggtitle(title) + ylab(yLabel) + xlab(xLabel)  + scale_shape_manual(values=c(0:2,5:6)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=18,hjust = 0.5, color="black"), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
        strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white"), # legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1,vjust=1), legend.text = element_text(size=18), legend.title = element_text(size=18))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Sialylation_overTime_PerSubject_ROCK.pdf", width = 6)


subsetData <- subset(mergedDataRock, TimeCategory == 'baseline')
shapiro_test(data = subsetData, vars = c("IgG1_Total.sialylated"))          # not normal distribution so will use nonparametric
twoSampleBarRock(data=subsetData, xData="Cohort", yData="IgG1_Total.sialylated", fillParam="Cohort", title="Sialylation\nbaseline", yLabel="% anti-H1 IgG1") + 
  scale_y_continuous(breaks = seq(0,100,2),limits=c(0,22)) + theme( legend.text = element_text(size=14), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Sialylation_baseline_ROCK.pdf", width = 4.5)
wilcox_test(IgG1_Total.sialylated ~ Cohort, data=subsetData)


##' ## ----------- Affinity analysis --------------------
#'
#'

subsetData <- subset(mergedDataRock, TimeCategory != 'late')
data=subsetData; xData="TimeCategory"; yData="Lo.Hi.HA.affinity"; fillParam="Cohort"; title="Affinity"; yLabel="Lo-Hi Affinity - H1N1pdm09"; xLabel="TimeCategory"; groupby = "Subject"
subsetData <- subsetData[order(subsetData$Subject, subsetData$TimePoint, decreasing = F),]
ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
  geom_bar(stat = "summary", aes_string(fill=fillParam),color="black",size=0.1) + 
  geom_path(aes_string(group=groupby), color="grey40", alpha=0.5) + 
  geom_point(size = 3, color="black", fill="black", alpha=0.5, aes(shape=TreatmentFactored)) + facet_wrap(fillParam ) +   
  scale_color_manual(values=c("#d9eafb", "#ff9a6a")) + scale_fill_manual(values=c("#d9eafb", "#ff9a6a")) + 
  ggtitle(title) + ylab(yLabel) + xlab(xLabel)  + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=18,hjust = 0.5, color="black"), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
        strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white"), # legend.position = "none",
        axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1,vjust=1), legend.text = element_text(size=14), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/LoHiAffinity_overTime_PerSubject_ROCK.pdf")


subsetData <- subset(mergedDataRock, TimeCategory == 'baseline')
twoSampleBarRock(data=subsetData, xData="Cohort", yData="Lo.Hi.HA.affinity", fillParam="Cohort", title="anti-HA Affinity", yLabel="Lo-Hi Affinity - H1N1pdm09") + 
  scale_y_continuous(breaks = seq(0,100,0.2),limits=c(0,1)) + theme(legend.text = element_text(size=14), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/LoHiAffinity_baseline_ROCK.pdf", width = 5)
wilcox_test(FCCXCL13_oW ~ Cohort, data=subsetData)


subsetData <- subset(mergedDataRock, TimeCategory != "oneWeek" &  Cohort!="")
FC_response <- dcast(subsetData, `Subject` + `Cohort`+ `TreatmentFactored` ~`TimeCategory` , value.var = c("Lo.Hi.HA.affinity"))
FC_response$FCaffinity_late <- FC_response$`late`/FC_response$`baseline`; FC_response$baseline <- FC_response$late <- NULL

data=FC_response; xData="Cohort"; yData="FCaffinity_late"; fillParam="Cohort"; title="Lo-Hi Affinity"; yLabel="Fold-change (late / baseline)"; position="left"
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T); names(overTime)[which(names(overTime) == 'x')] <- yData
justforttest <- data[, c(xData,yData)]; fit <- rstatix::wilcox_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p;   annotationInfo <- paste0("P = ", round(pValue, 2))
my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28))) 
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.8)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))  + 
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5) + 
  #geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity") + 
  geom_violin(draw_quantiles = 0.5, size=0.1) + 
  geom_point(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5, position = position_jitter(width=0.3, seed=58)) +   #48
  ggtitle(title) + ylab(yLabel) +  theme_bw() + scale_shape_manual(values=c(21:25)) + labs(col="Cohort",shape="Treatment") + 
  theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=28,hjust = 0.5),axis.text.x=element_text(angle=45,vjust=1,hjust=1)) + annotation_custom(my_grob) + 
scale_y_continuous(breaks = seq(0,10,0.25), limits = c(0.75,1.75))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/LoHiAffinity_foldchange_ROCK.pdf", width = 4.5)
wilcox_test(FCaffinity_late ~ Cohort, data=FC_response)



##' ## ----------- irAE analysis --------------------
#'
#'

subsetData <- mergedDataRock.nolate[which(!is.na(mergedDataRock.nolate$irAE) & mergedDataRock.nolate$irAE != "" ), ]
subsetData <- subset(subsetData, TimeCategory=="oneWeek" & Cohort=="aPD1")
twoSampleBarRock(data = subsetData, yData="FCtfh_oW", yLabel = "Fold-change at one week", title="ICOS+CD38+ cTfh", xData = "irAE",fillParam = "irAE",position = 'right') + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))

twoSampleBarRock(data = subsetData, yData="ICOS.CD38..cTfh.medFI.aIgG4", yLabel = "medianFI aIgG4 - baseline", title="ICOS+CD38+ cTfh", 
                 xData = "irAE",fillParam = "irAE",position = 'left') + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28),axis.text.x = element_text(angle=0, hjust=0.5,vjust=0,face="plain")) + 
  coord_cartesian(ylim = c(0,800))

twoSampleBarRock(data = subsetData, yData="ICOS.CD38..cTfh.medFI.ICOS", yLabel = "medianFI ICOS - baseline", title="ICOS+CD38+ cTfh", 
                 xData = "irAE",fillParam = "irAE",position = 'left') + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28),axis.text.x = element_text(angle=0, hjust=0.5,vjust=0,face="plain")) + 
  coord_cartesian(ylim = c(0,20000))

twoSampleBarRock(data = subsetData, yData="FChai_late", yLabel = "Fold-change in HAI titers", title="HAI titers", 
                 xData = "irAE",fillParam = "irAE",position = 'left') + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28),axis.text.x = element_text(angle=0, hjust=0.5,vjust=0,face="plain")) + 
  coord_cartesian(ylim = c(0,3))


# ---------- exclude any combination treatment --------------------

subsetData <- mergedDataRock.nolate[which(!is.na(mergedDataRock.nolate$irAE) & mergedDataRock.nolate$irAE != "" ), ]
subsetData <- subset(subsetData, Treatment == "Nivolumab" | Treatment == "Pembrolizumab")
subsetData <- subset(subsetData, TimeCategory=="oneWeek" & Cohort=="aPD1")

twoSampleBarRock(data = subsetData, yData="ICOS.CD38..cTfh.medFI.aIgG4", yLabel = "medianFI aIgG4 - baseline", title="ICOS+CD38+ cTfh", 
                 xData = "irAE",fillParam = "irAE",position = 'left') + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28),axis.text.x = element_text(angle=0, hjust=0.5,vjust=0,face="plain")) + 
  coord_cartesian(ylim = c(0,800)) + theme( legend.text = element_text(size=14), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_baselineTfh_aIgG4_ROCK.pdf", width=4)

twoSampleBarRock(data = subsetData, yData="ICOS.CD38..cTfh.medFI.ICOS", yLabel = "medianFI ICOS - baseline", title="ICOS+CD38+ cTfh", 
                 xData = "irAE",fillParam = "irAE",position = 'left') + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28),axis.text.x = element_text(angle=0, hjust=0.5,vjust=0,face="plain")) + 
  coord_cartesian(ylim = c(0,20000)) + theme( legend.text = element_text(size=14), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_baselineTfh_ICOS_ROCK.pdf", width=4)

twoSampleBarRock(data = subsetData, yData="FCtfh_oW", yLabel = "Fold-change at one week", title="ICOS+CD38+ cTfh", 
                 xData = "irAE",fillParam = "irAE",position = 'left') + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))+ theme( legend.text = element_text(size=14), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_FCtfh_ROCK.pdf", width=5)

twoSampleBarRock(data = subsetData, yData="FChai_late", yLabel = "Fold-change in HAI titers", title="HAI titers", 
                 xData = "irAE",fillParam = "irAE",position = 'left') + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28),axis.text.x = element_text(angle=0, hjust=0.5,vjust=0,face="plain")) + 
  coord_cartesian(ylim = c(0,3))+ theme( legend.text = element_text(size=14), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_HAItiters_ROCK.pdf", width=5)






