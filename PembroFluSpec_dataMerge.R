#+ fig.width=8, fig.height=8
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
source('D:/Pembro-Fluvac/Analysis/Ranalysis/PembroFluSpec_PlottingFunctions.R')
sessionInfo()

#' ## ----------- merge & qc data from spreadsheets  --------------------
#'
setwd("D:/Pembro-Fluvac/Analysis")
mergedData <- Tflow <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_freqParent_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
print("Data by year and cohort for T cell flow cytometry")
table(mergedData$Cohort, mergedData$TimeCategory, mergedData$Year)
fit <- lm(data=mergedData[ , - which( colnames(mergedData) == "Label" | colnames(mergedData) == "Subject" | colnames(mergedData) == "Cohort" | 
                                        colnames(mergedData) == "TimeCategory" | colnames(mergedData) == "TimePoint")], formula = TflowBatch  ~ .)
a <- as.data.frame(summary(fit)$coefficients)
a <- rownames(a[which(a$`Pr(>|t|)` < 0.05), ] )
if (length (a) > 0) { colnames(mergedData) [ grep( paste(a, collapse="|"), colnames(mergedData)) ] <- paste0( colnames(mergedData)[ grep( paste(a, collapse="|"), colnames(mergedData)) ], ".batchEffect" ) }

serology <- read.csv(file= "D:/Pembro-Fluvac/Analysis/mergedData/SerologyMeasurements.csv", stringsAsFactors = F, header = T)
print("Data by year and cohort for serology")
table(serology$Cohort, serology$TimeCategory, serology$Year)
temp1 <- merge(x = mergedData, y=serology, all = T, suffixes = c(".Tflow",".Serology"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year'))

serologyKd <- read.csv(file="D:/Pembro-Fluvac/Analysis/mergedData/Serology_Kd.csv", header=T)
print("Data by year and cohort for 1/Kd measurements")
table(serologyKd$Cohort, serologyKd$TimeCategory, serologyKd$Year)
temp1 <- merge(x = temp1, y=serologyKd, all = T, suffixes = c(".Tflow",".SerologyKd"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year','Visit'))


Bflow <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_freqParent_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
print("Data by year and cohort for B cell flow cytometry")
table(Bflow$Cohort, Bflow$TimeCategory, Bflow$Year)
BflowBATCH <- Bflow[, -c( grep( colnames(Bflow), pattern=paste0(c("Naive","Lymphocytes", "CD4","CD19hi_..FreqParent","CD19hi_CD11c","CD19hi_CD16","CD19hi_CD23","CD1hi9_CD25",
                                                             "CD19hi_CD32","CD19hi_CD38lo...FreqParent","CD19hi_CD80","CD19hi_CD86","CD19hi_IgM","CD19hi_CD138",
                                                             "CD19hi_CXCR4","CD19hi_CD71","CD19hi_Ki67..CXCR4","CD19hi_Ki67..Ki67","CD38lo...FreqParent",
                                                             "CD19hi_Tbet"), collapse ="|")) )]    # cut out columns to test for batch fx, else variables > observations
fit <- lm(data=BflowBATCH[ , - which( colnames(BflowBATCH) == "Label" | colnames(BflowBATCH) == "Subject" | colnames(BflowBATCH) == "Cohort" | 
                                        colnames(BflowBATCH) == "TimeCategory" | colnames(BflowBATCH) == "TimePoint")], 
          formula = BflowBatch  ~ .)
a <- as.data.frame(summary(fit)$coefficients)
a <- rownames(a[which(a$`Pr(>|t|)` < 0.05), ] )
if (length (a) > 0) { colnames(Bflow) [ grep( paste(a, collapse="|"), colnames(Bflow)) ] <- paste0( colnames(Bflow)[ grep( paste(a, collapse="|"), colnames(Bflow)) ], ".batchEffect" )}
temp2 <- merge(x = temp1, y=Bflow, all = T, suffixes = c(".Tflow",".Bflow"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year'))


Tmfi <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_medianFI_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
Tmfi <- Tmfi[, -c(grep( colnames(Tmfi), pattern=paste0(c("Dead","medfi.CD20","medfi.CD4","Freq","ICOSlo"), collapse ="|")))]  # cut out columns to test for batch fx, else variables > observations
fit <- lm(data=Tmfi[ , - which( colnames(Tmfi) == "Label" | colnames(Tmfi) == "Subject" | colnames(Tmfi) == "Cohort" | 
                                   colnames(Tmfi) == "TimeCategory" | colnames(Tmfi) == "TimePoint")], 
          formula = TmfiBatch  ~ .)
a <- as.data.frame(summary(fit)$coefficients)
a <- rownames(a[which(a$`Pr(>|t|)` < 0.05), ] )
if (length (a) > 0) { colnames(Tmfi) [ grep( paste(a, collapse="|"), colnames(Tmfi)) ] <- paste0( colnames(Tmfi)[ grep( paste(a, collapse="|"), colnames(Tmfi)) ], ".batchEffect" )}

temp3 <- merge(x = temp2, y=Tmfi, all = T, suffixes = c(".one",".Tmedfi"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year'))

demog <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/clinicalMetadata.csv", stringsAsFactors = F, header=T); 
demog$TimePoint <- demog$Visit <- NULL
demog$DateFirstIRAE <- strptime(demog$DateFirstIRAE, format = "%Y%m%d")
demog$DateFirstFlowFile <- strptime(demog$DateFirstFlowFile, format = "%Y%m%d")

temp4 <- merge(x = temp3, y=demog, all=T, suffixes = c(".two",".demog"), by = c('Subject', 'TimeCategory', 'Year','Cohort'))


Bmfi <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_Bmfi_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
fit <- lm(data=Bmfi[ , - which( colnames(Bmfi) == "Label" | colnames(Bmfi) == "Subject" | colnames(Bmfi) == "Cohort" | 
                                    colnames(Bmfi) == "TimeCategory" | colnames(Bmfi) == "TimePoint")], 
          formula = BmfiBatch  ~ .)
a <- as.data.frame(summary(fit)$coefficients)
a <- rownames(a[which(a$`Pr(>|t|)` < 0.05), ] )
if (length (a) > 0) { colnames(Bmfi) [ grep( paste(a, collapse="|"), colnames(Bmfi)) ] <- paste0( colnames(Bmfi)[ grep( paste(a, collapse="|"), colnames(Bmfi)) ], ".batchEffect" )}
temp5 <- merge(x = temp4, y=Bmfi, all = T, suffixes = c(".one",".TBpmedfi"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year'))


TBpmfi <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_TmfiFromB_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
fit <- lm(data=TBpmfi[ , - which( colnames(TBpmfi) == "Label" | colnames(TBpmfi) == "Subject" | colnames(TBpmfi) == "Cohort" | 
                                  colnames(TBpmfi) == "TimeCategory" | colnames(TBpmfi) == "TimePoint")], 
          formula = TfromBBatch  ~ .)
a <- as.data.frame(summary(fit)$coefficients)
a <- rownames(a[which(a$`Pr(>|t|)` < 0.05), ] )
if (length (a) > 0) { colnames(TBpmfi) [ grep( paste(a, collapse="|"), colnames(TBpmfi)) ] <- paste0( colnames(TBpmfi)[ grep( paste(a, collapse="|"), colnames(TBpmfi)) ], ".batchEffect" )}
temp6 <- merge(x = temp5, y=TBpmfi, all = T, suffixes = c(".one",".TBpmedfi"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year'))


mergedData <- temp6

index <- demog[which(demog$Current.Immunotherapy == 'Ipi/Nivo'),"Subject"]
mergedData$Cohort[which(mergedData$Subject %in% index)] <- "nonPD1"      # exclude combo therapy


rm(Bflow); rm(temp1); rm(temp2); rm(temp3); rm(temp4); rm(temp5); rm(temp6)

##' ## ----------- calculate fold-changes for HAI, Tfh, and PB --------------------

mergedData$logHAI <- log(mergedData$H1N1pdm09.HAI.titer)
mergedData$X5hi_PD1hi_ICOSCD38_FreqX5 <- mergedData$cTfh_PD1..ICOS.CD38....FreqParent * mergedData$cTfh_PD1....FreqParent / 100
subsetData <- subset(mergedData, TimeCategory != "oneWeek" & Cohort != "nonPD1" )
FC_response <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("logHAI"))
FC_response$FChai_late <- FC_response$`late`/FC_response$`baseline`
FC_response2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("H1N1pdm09.HAI.titer"))
FC_response2$FCH1N1_hai_late <- FC_response2$`late`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject")
FC_response2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("H3N2.HAI.titer"))
FC_response2$FCH3N2_hai_late <- FC_response2$`late`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
  FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject")
FC_response2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("FluB.HAI.titer"))
FC_response2$FCFluB_hai_late <- FC_response2$`late`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
  FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject")

subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" )
FC_response2 <- dcast( subset(subsetData, cTfh_ICOShiCD38hi_..FreqParent > 0), `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")) 
FC_response2$FCtfh_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject", all.x=T)
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD27hiCD38hi_..FreqParent")) 
FC_response2$FCPB_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject", all.x=T)
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("Plasma.CXCL13..pg.mL.")) 
FC_response2$FCCXCL13_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject", all.x=T)
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgG1_Total.sialylated")) 
FC_response2$IgG1sial_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject", all.x=T)
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("X5hi_PD1hi_ICOSCD38_FreqX5")) 
FC_response2$X5PD1_ICOSCD38_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
   FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject", all.x=T)


FC_response <- FC_response[, c("Subject","Cohort","FChai_late","FCtfh_oW", "FCPB_oW", "FCCXCL13_oW", "IgG1sial_oW","X5PD1_ICOSCD38_oW","FCH1N1_hai_late","FCH3N2_hai_late", "FCFluB_hai_late" )]
mergedData <- merge(x = mergedData, y= FC_response, all=T, by = c('Subject', 'Cohort'))

mergedData$ASC_ABCratio <- mergedData$IgDlo_CD71hi_CD20loCD71hi...FreqParent / mergedData$IgDlo_CD71hi_ActBCells...FreqParent

mergedData$TimeCategory <- factor(mergedData$TimeCategory, levels = c("baseline", "oneWeek","late"))
mergedData$Cohort[which(mergedData$Cohort == "aPD1")] <- "anti-PD-1"
demog$Cohort[which(demog$Cohort == "aPD1")] <- "anti-PD-1"
mergedData$Cohort <- factor(mergedData$Cohort, levels = c("Healthy", "anti-PD-1","nonPD1"))
mergedData$dummy <- "dummy"


##' -------------------- Double coding -----------------------


dbCode <- read.csv(file="D:/Pembro-Fluvac/Analysis/doubleCoding.csv")
dbCode <- dbCode[!duplicated(dbCode$Subject),]
dbCode$Concat <- NULL

temp <- demog[,1:5]
temp <- merge(x=temp, y=dbCode, all=T, by="Subject")

totNumAPD1 <- length(grep(pattern="anti-PD-1", x = temp$Cohort))          # 32 people in the aPD1 arm
totNumHealthy <- length(grep(pattern="Healthy", x = temp$Cohort))    # 27 people in the Healthy arm
notlabeled <- temp[which(is.na(temp$dbCode)),]                       
   notlabeled_aPD1 <- nrow(notlabeled[grep(x=notlabeled$Subject, pattern="19616"),])      # 22 in aPD1 were not labeled
   notlabeled_HC <- nrow(notlabeled[grep(x=notlabeled$Subject, pattern="FS"),])           # 13 in HC were not labeled
labeled <- temp[which(!is.na(temp$dbCode)),]                         
   labeled_aPD1 <- nrow(labeled[grep(x=labeled$Subject, pattern="19616"),])      # 22 in aPD1 were not labeled
   labeled_HC <- nrow(labeled[grep(x=labeled$Subject, pattern="FS"),])           # 13 in HC were not labeled   

aPD1_sequence <- (labeled_aPD1+1):(totNumAPD1)
HC_sequence <- (labeled_HC+1):(totNumHealthy)

notlabeled[grep(pattern = "19616", x=notlabeled$Subject),"dbCode"] <- paste("anti-PD1-treated adult", aPD1_sequence)
notlabeled[grep(pattern = "FS", x=notlabeled$Subject),"dbCode"] <- paste("Healthy adult", HC_sequence)

fullLabel <- rbind(labeled,notlabeled)
temp <- merge(x = temp, y=fullLabel[,c("Subject","dbCode")], all=T, by="Subject")
temp$dbCode <- temp$dbCode.y
temp$dbCode.x <- temp$dbCode.y <- NULL
mergedData.x <- merge(x=mergedData, y=temp[,c("Subject","dbCode", "GSMid", "GEOlabel")], all.x=T, by="Subject")
mergedData <- mergedData.x
rm(temp, mergedData.x, fullLabel, notlabeled,labeled,aPD1_sequence, HC_sequence)

##' ## ----------- Save mergedData as RDS  --------------------

# saveRDS(mergedData, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergedData.RDS")
# saveRDS(demog, file="D:/Pembro-Fluvac/Analysis/mergedData/demog.RDS")
# saveRDS(serology, file="D:/Pembro-Fluvac/Analysis/mergedData/serology.RDS")
# write.csv(mergedData, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergedData.csv")


##' ## ----------- Dataset overview --------------------


dataAvail <- mergedData[ , - grep(paste(c("TimePoint", "dummy", "Visit","Label"), collapse = "|"), colnames(mergedData))]
temp <- reshape(dataAvail, idvar = c('Subject', 'Cohort','Year'), timevar = c('TimeCategory'),direction = 'wide')
saveCohort <- temp[which(temp$Cohort != "nonPD1"),"Cohort"]
saveYear <- temp[which(temp$Cohort != "nonPD1"),"Year"]
simplifiedDataAvail <- temp[which(temp$Cohort != "nonPD1"), c('Subject',  
                                                              'cTfh_ICOShiCD38hi_CCR6hi...FreqParent.baseline' ,          # Tflow_freq
                                                              'Plasma.CXCL13..pg.mL..baseline' ,                          # Serology
                                                              'IgDlo_CD71hi_ActBCells.Tbet....FreqParent.baseline' ,      # Bflow_freq
                                                              'CD19hi_CD27hiCD38hi_CD20lo.PB..._medfi.Foxp3..baseline' ,  # Tflow_medFI
                                                              'Sex.baseline' ,                                            # Demographics
                                                              'CD20loCD71hi...medfi.Blimp1..baseline')]                   # Bflow_medFI

rownames(simplifiedDataAvail) <- simplifiedDataAvail$Subject; simplifiedDataAvail$Subject <- NULL;
simplifiedDataAvail <- as.data.frame(!is.na(simplifiedDataAvail));  simplifiedDataAvail <- simplifiedDataAvail * 2 
simplifiedDataAvail <- cbind(simplifiedDataAvail, Cohort = (as.numeric(saveCohort)-1)*2, Year = as.numeric(saveYear)-1)
names(simplifiedDataAvail) <- c("Tflow_freq","Serology","Bflow_freq", "Tflow_medFI","Demographics","Bflow_medFI","Cohort", "Year")
pheatmap(as.data.frame(t(simplifiedDataAvail)), cluster_col = F,color=gray.colors(100, start=1,end=0), fontsize_row = 18, fontsize_col = 6, main = "Data available by subject"
         # , filename = "D:/Pembro-Fluvac/Analysis/Images/dataAvailable.pdf"
)
# dev.off()

