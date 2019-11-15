library("grid")
library("ggplot2")
library("gplots")
library("viridis")
library("reshape2")
library("gridExtra")
library("RColorBrewer")
library("scales")


source('D:/Pembro-Fluvac/Analysis/PembroFluSpec_Ranalysis_files/PembroFluSpec_PlottingFunctions.R')


## ------------------------------------------------- merge & qc data from spreadsheets  ----------------------------------------------------------

setwd("D:/Pembro-Fluvac/Analysis")
mergedData <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_freqParent_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
print("Data by year and cohort for T cell flow cytometry")
table(mergedData$Cohort, mergedData$TimeCategory, mergedData$Year)
fit <- lm(data=mergedData[ , - which( colnames(mergedData) == "Label" | colnames(mergedData) == "Subject" | colnames(mergedData) == "Cohort" | 
                                        colnames(mergedData) == "TimeCategory" | colnames(mergedData) == "TimePoint")], 
          formula = TflowBatch  ~ .)
a <- as.data.frame(summary(fit)$coefficients)
a <- rownames(a[which(a$`Pr(>|t|)` < 0.05), ] )
if (length (a) > 0) { colnames(mergedData) [ grep( paste(a, collapse="|"), colnames(mergedData)) ] <- paste0( colnames(mergedData)[ grep( paste(a, collapse="|"), colnames(mergedData)) ], ".batchEffect" ) }

serology <- read.csv(file= "D:/Pembro-Fluvac/Analysis/mergedData/SerologyMeasurements.csv", stringsAsFactors = F, header = T)
print("Data by year and cohort for serology")
table(serology$Cohort, serology$TimeCategory, serology$Year)
temp1 <- merge(x = mergedData, y=serology, all = T, suffixes = c(".Tflow",".Serology"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year'))


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

demog <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/clinicalMetadata.csv", stringsAsFactors = F, header=T); demog$TimePoint <- demog$Visit <- NULL
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

mergedData$TimeCategory <- factor(mergedData$TimeCategory, levels = c("baseline", "oneWeek","late"))
mergedData$Cohort <- factor(mergedData$Cohort, levels = c("Healthy", "aPD1","nonPD1"))
mergedData$dummy <- "dummy"

rm(Bflow); rm(temp1); rm(temp2); rm(temp3); rm(temp4); rm(temp5); rm(temp6)

# write.csv(mergedData, file = "D:/Pembro-Fluvac/Analysis/mergedData/allMergedData.csv")
   
## ------------------------------------------------- Cohort description ----------------------------------------------------------

table(mergedData$Cohort, mergedData$Sex, mergedData$Year)/3  # divide by 3 because demog data repeated 3x per subject
hist(mergedData[which(mergedData$Cohort == "aPD1" & mergedData$TimeCategory == "baseline"), 'Cycle.of.Immunotherapy'], breaks = 50, col = "grey90",main = "Cycle of IT")   
median(mergedData[which(mergedData$Cohort == "aPD1" & mergedData$TimeCategory == "baseline"), 'Cycle.of.Immunotherapy'], na.rm = T)  



## ------------------------------------------------- Global overview of phenotypic differences ----------------------------------------------------------

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_NaiveB...FreqParent"))
prePostTimeAveraged(melted, title = "Naive B frequencies averaged", xLabel = NULL, yLabel = "NaiveB (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, CD19hi_NaiveB...FreqParent ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_Ki67....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+Ki67+ frequencies averaged", xLabel = NULL, yLabel = "CD19+Ki67+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, CD19hi_Ki67....FreqParent ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_CD71....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+CD71+ frequencies averaged", xLabel = NULL, yLabel = "CD19+CD71+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, CD19hi_CD71....FreqParent ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_CD80....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+CD80+ frequencies averaged", xLabel = NULL, yLabel = "CD19+CD80+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, CD19hi_CD80....FreqParent ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_CD86....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+CD86+ frequencies averaged", xLabel = NULL, yLabel = "CD19+CD86+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, CD19hi_CD86....FreqParent ~ Cohort + Year + TimeCategory))


subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_IgM....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+IgM+ frequencies averaged", xLabel = NULL, yLabel = "CD19+IgM+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_CD138....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+CD138+ frequencies averaged", xLabel = NULL, yLabel = "CD19+CD138+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_Tbet....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+Tbet+ frequencies averaged", xLabel = NULL, yLabel = "CD19+Tbet+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))





subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("Live_CD3..CD4..ICOShi...FreqParent"))
prePostTimeAveraged(melted, title = "CD4+ICOS+ frequencies averaged", xLabel = NULL, yLabel = "CD4+ICOS+ (freq CD4)") #+ scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, Live_CD3..CD4..ICOShi...FreqParent ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" & Year == 3 )  # given strong influence of Year in the linear model 
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("Live_CD3..CD4..CTLA4hi...FreqParent"))
prePostTimeAveraged(melted, title = "CD4+CTLA4+ frequencies averaged", xLabel = NULL, yLabel = "CD4+CTLA4+ (freq CD4)") #+ scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, Live_CD3..CD4..CTLA4hi...FreqParent ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("Live_CD3..CD4..CD38hi...FreqParent"))
prePostTimeAveraged(melted, title = "CD4+CD38hi frequencies averaged", xLabel = NULL, yLabel = "CD4+CD38hi (freq CD4)") #+ scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, Live_CD3..CD4..CD38hi...FreqParent ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("Live_CD3..CD4..Foxp3hi...FreqParent"))
prePostTimeAveraged(melted, title = "CD4+Foxp3hi frequencies averaged", xLabel = NULL, yLabel = "CD4+Foxp3hi (freq CD4)") #+ scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, Live_CD3..CD4..Foxp3hi...FreqParent ~ Cohort + Year + TimeCategory))
summary(lm( data = mergedData, Live_CD3..CD4..Foxp3hi...FreqParent ~ Cohort + Year + TimeCategory + Live_CD3..CD4..CTLA4hi...FreqParent + Live_CD3..CD4..ICOShi...FreqParent))



subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("Live_CD3..CD4..Ki67hi...FreqParent"))
prePostTimeAveraged(melted, title = "CD4+Ki67hi frequencies averaged", xLabel = NULL, yLabel = "CD4+Ki67hi (freq CD4)") #+ scale_y_continuous(breaks=seq(0,99,0.25))
summary(lm( data = mergedData, Live_CD3..CD4..Ki67hi...FreqParent ~ Cohort + Year + TimeCategory))


## ------------------------------------------------- Global overview of phenotypic differences ----------------------------------------------------------

# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$cTfh_ICOShiCD38hi_..FreqParent <- subsetData$cTfh_ICOShiCD38hi_..FreqParent * subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /10000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_..FreqParent"))
prePostTimeAveraged(melted, title = "+/+ cTfh frequencies averaged", xLabel = NULL, yLabel = "+/+ cTfh (freq CD4)") 


# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$cTfh_ICOShiCD38hi_CTLA4hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_CTLA4hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent * 
  subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /1000000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CTLA4hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ CTLA4hi frequencies averaged", xLabel = NULL, yLabel = "+/+ CTLA4hi (freq CD4)") 
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_CTLA4hi...FreqParent ~ Cohort + Year + TimeCategory))


# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" & Year == 3)
subsetData$cTfh_ICOShiCD38hi_Ki67hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_Ki67hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent * 
  subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /1000000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_Ki67hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ Ki67hi frequencies averaged", xLabel = NULL, yLabel = "+/+ Ki67hi (freq CD4)") 
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_Ki67hi...FreqParent ~ Cohort + Year + TimeCategory))  # mergedData because of influence of Year



# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$cTfh_ICOShiCD38hi_CXCR3hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_CXCR3hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent * 
  subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /1000000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CXCR3hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ CXCR3hi frequencies averaged", xLabel = NULL, yLabel = "+/+ CXCR3hi (freq CD4)") 
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_CXCR3hi...FreqParent ~ Cohort + Year + TimeCategory))


# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$cTfh_ICOShiCD38hi_CD25hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_CD25hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent * 
  subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /1000000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CD25hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ CD25hi frequencies averaged", xLabel = NULL, yLabel = "+/+ CD25hi (freq CD4)") 
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_CD25hi...FreqParent ~ Cohort + Year + TimeCategory))





subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_NaiveB...FreqParent"))
prePostTimeAveraged(melted, title = "Naive B frequencies averaged", xLabel = NULL, yLabel = "NaiveB (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))





## ------------------------------------------------- PB and Tfh just baseline frequencies ----------------------------------------------------------


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBox(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: cTfh +/+ at d0", yLabel="ICOS+CD38+ cTfh frequency at d0")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBox(data=subsetData, xData="Cohort", yData="CD27hiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: PB at d0", yLabel="CD19+IgD-CD27++CD38++ frequency at d0")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBox(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_ActBCells...FreqParent", fillParam="Cohort", title="All yrs: ABC at d0", yLabel="ABC (freq of CD19+IgD-CD71+) at d0")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBox(data=subsetData, xData="Cohort", yData="CD19hi_NaiveB...FreqParent", fillParam="Cohort", title="All yrs: Naive B at d0", yLabel="CD19+IgDhi frequency at d0")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBox(data=subsetData, xData="Cohort", yData="Live_CD3..CD4..ICOShi...FreqParent", fillParam="Cohort", title="All yrs: CD4+ICOS+ at d0", yLabel="CD4+ICOS+ frequency at d0")





## ------------------------------------------------- PB and Tfh just day 7 frequencies ----------------------------------------------------------

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek" & Year == "3" & cTfh_ICOShiCD38hi_..FreqParent > 0)
twoSampleBox(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="Yr3: ICOS+CD38+ cTfh at d7", yLabel="ICOS+CD38+ cTfh frequency at d7") + scale_y_continuous(breaks=seq(1:12), limits=c(0.5,11))

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek" & Year == "3")
twoSampleBox(data=subsetData, xData="Cohort", yData="CD19_CD27.CD38....FreqParent", fillParam="Cohort", title="Yr3: Plasmablast at d7", yLabel="CD19+CD27++CD38++ frequency at d7") + scale_y_continuous(breaks=seq(0,99,1))

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 0)
twoSampleBox(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: ICOS+CD38+ cTfh at d7", yLabel="ICOS+CD38+ cTfh frequency at d7") + scale_y_continuous(breaks=seq(0,99,1))

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBox(data=subsetData, xData="Cohort", yData="CD27hiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: Plasmablast at d7", yLabel="CD19+CD27++CD38++ frequency at d7")

# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
subsetData$CD27hiCD38hi_..FreqParent <-  subsetData$CD27hiCD38hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 100
twoSampleBox(data=subsetData, xData="Cohort", yData="CD27hiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: PB (freq CD19) at d7", yLabel="CD19+CD27++CD38++ frequency at d7")



subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBox(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_CD20loCD71hi...FreqParent", fillParam="Cohort", title="All yrs: Plasmablast at d7", yLabel="IgDloCD71+CD20lo frequency at d7")

# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 10000
twoSampleBox(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_CD20loCD71hi...FreqParent", fillParam="Cohort", title="All yrs: Plasmablast at d7", yLabel="IgDloCD71+CD20lo frequency at d7")



subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBox(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_ActBCells...FreqParent", fillParam="Cohort", title="All yrs: ABC at d7", yLabel="ABC frequency at d7")
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
subsetData$IgDlo_CD71hi_ActBCells...FreqParent <-  subsetData$IgDlo_CD71hi_ActBCells...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 10000
twoSampleBox(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_ActBCells...FreqParent", fillParam="Cohort", title="All yrs: ABC at d7", yLabel="ABC frequency at d7")









## ------------------------------------------------- Tfh and PB Fold-change at day 7  ----------------------------------------------------------

subsetData <- subset(mergedData, TimeCategory != "late" & Year == "3" & Cohort != "nonPD1")
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="Yr3: cTfh response at one week", yLabel="ICOS+CD38+ cTfh fold-change") + scale_y_continuous(breaks=seq(0,99,1))

subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & cTfh_ICOShiCD38hi_..FreqParent > 0)
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="All yrs: cTfh response at one week", yLabel="ICOS+CD38+ cTfh fold-change") + scale_y_continuous(breaks=seq(0,99,1))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResponse_foldchange_AllYrs.pdf", device="pdf")


subsetData <- subset(mergedData, TimeCategory != "late"  & Cohort != "nonPD1")
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD27hiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="AllYrs: PB response at one week", yLabel="CD19+CD27++CD38++ (freq. IgD-) fold-change") + scale_y_continuous(breaks=seq(0,99,1))

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Year == "3" & Cohort != "nonPD1")
subsetData$CD27hiCD38hi_..FreqParent <-  subsetData$CD27hiCD38hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent /100
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD27hiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="Yr3: PB response (freq CD19) at one week", yLabel="CD19+CD27++CD38++ fold-change") + scale_y_continuous(breaks=seq(0,99,1))





subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent > 0)
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="All yrs: PB response at one week", yLabel="IgDloCD71+CD20lo fold-change") + scale_y_continuous(breaks=seq(0,99,1))

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent > 0)
subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="All yrs: PB response (freq CD19) at one week", yLabel="IgDloCD71+CD20lo fold-change") + scale_y_continuous(breaks=seq(0,99,1))


subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent > 0)
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_ActBCells...FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="AllYrs: ABC response FC at one week", yLabel="ABC fold-change")

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent > 0)
subsetData$IgDlo_CD71hi_ActBCells...FreqParent <-  subsetData$IgDlo_CD71hi_ActBCells...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_ActBCells...FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="AllYrs: ABC response (freq CD19) FC at one week", yLabel="ABC fold-change")






subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1")
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC_Tfh <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
FC_Tfhresponse2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD27hiCD38hi_..FreqParent")); 
FC_Tfhresponse$PB <- FC_Tfhresponse2$`oneWeek`/FC_Tfhresponse2$`baseline`
bivScatter(data1 = FC_Tfhresponse[which(FC_Tfhresponse$Cohort == "Healthy"),], data2 = FC_Tfhresponse[ which(FC_Tfhresponse$Cohort == "aPD1") ,], 
           name1 = "Healthy", name2 = "aPD1", xData = "FC_Tfh", yData = "PB", fillParam = "Cohort", title = "FC PB (freqParent) vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "CD19+CD27++CD38++ PB fold-change") + scale_y_continuous(limits=c(0,12))

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1"  )
subsetData$CD27hiCD38hi_..FreqParent <-  subsetData$CD27hiCD38hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC_Tfh <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
FC_Tfhresponse2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD27hiCD38hi_..FreqParent")); 
FC_Tfhresponse$PB <- FC_Tfhresponse2$`oneWeek`/FC_Tfhresponse2$`baseline`
bivScatter(data1 = FC_Tfhresponse[which(FC_Tfhresponse$Cohort == "Healthy"),], data2 = FC_Tfhresponse[ which(FC_Tfhresponse$Cohort == "aPD1") ,], 
           name1 = "Healthy", name2 = "aPD1", xData = "FC_Tfh", yData = "PB", fillParam = "Cohort", title = "FC PB (as freq CD19) vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "CD19+CD27++CD38++ PB fold-change") + scale_y_continuous(limits=c(0,12))



  
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1"  & IgDlo_CD71hi_CD20loCD71hi...FreqParent > 0 )
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC_Tfh <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
FC_Tfhresponse2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent")); 
FC_Tfhresponse$PB <- FC_Tfhresponse2$`oneWeek`/FC_Tfhresponse2$`baseline`
bivScatter(data1 = FC_Tfhresponse[which(FC_Tfhresponse$Cohort == "Healthy"),], data2 = FC_Tfhresponse[ which(FC_Tfhresponse$Cohort == "aPD1") ,], 
           name1 = "Healthy", name2 = "aPD1", xData = "FC_Tfh", yData = "PB", fillParam = "Cohort", title = "FC PB vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "IgDloCD71+CD20lo PB fold-change")

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent > 0 )
subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC_Tfh <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
FC_Tfhresponse2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent")); 
FC_Tfhresponse$PB <- FC_Tfhresponse2$`oneWeek`/FC_Tfhresponse2$`baseline`
bivScatter(data1 = FC_Tfhresponse[which(FC_Tfhresponse$Cohort == "Healthy"),], data2 = FC_Tfhresponse[ which(FC_Tfhresponse$Cohort == "aPD1") ,], 
           name1 = "Healthy", name2 = "aPD1", xData = "FC_Tfh", yData = "PB", fillParam = "Cohort", title = "FC PB (as freq CD19) vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "IgDloCD71+CD20lo PB fold-change") + scale_y_continuous(limits=c(0,12))





subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1"  & IgDlo_CD71hi_ActBCells...FreqParent > 0 )
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC_Tfh <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
FC_Tfhresponse2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_ActBCells...FreqParent")); 
FC_Tfhresponse$ABC <- FC_Tfhresponse2$`oneWeek`/FC_Tfhresponse2$`baseline`
bivScatter(data1 = FC_Tfhresponse[which(FC_Tfhresponse$Cohort == "Healthy"),], data2 = FC_Tfhresponse[ which(FC_Tfhresponse$Cohort == "aPD1") ,], 
           name1 = "Healthy", name2 = "aPD1", xData = "FC_Tfh", yData = "ABC", fillParam = "Cohort", title = "FC ABC vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "IgDloCD71+CD20hi ABC fold-change") + scale_y_continuous(limits=c(0,4))

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent > 0 )
subsetData$IgDlo_CD71hi_ActBCells...FreqParent <-  subsetData$IgDlo_CD71hi_ActBCells...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC_Tfh <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
FC_Tfhresponse2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_ActBCells...FreqParent")); 
FC_Tfhresponse$ABC <- FC_Tfhresponse2$`oneWeek`/FC_Tfhresponse2$`baseline`
bivScatter(data1 = FC_Tfhresponse[which(FC_Tfhresponse$Cohort == "Healthy"),], data2 = FC_Tfhresponse[ which(FC_Tfhresponse$Cohort == "aPD1") ,], 
           name1 = "Healthy", name2 = "aPD1", xData = "FC_Tfh", yData = "ABC", fillParam = "Cohort", title = "FC ABC (as freq CD19) vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "IgDloCD71+CD20hi ABC fold-change") + scale_y_continuous(limits=c(0,12))



## ------------------------------------------------- PB vs cTfh response correlations ----------------------------------------------------------

## ***************    CD27++CD38++     ********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
subsetData1 <- subset(oneWeek, Cohort == "Healthy" & Year == "3")    
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & Year == "3" )    
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD27hiCD38hi_..FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "Yr3: cTfh vs PB response at oneWeek", xLabel= "CD27++CD38++ (freq IgD-)", yLabel = "ICOS+CD38+ cTfh frequency") + scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1))
univScatter(data = oneWeek, xData = "CD27hiCD38hi_..FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "AllYrs: cTfh vs CD27++CD38++ at oneWeek", xLabel= "CD27++CD38++ (freq IgD-)", yLabel = "ICOS+CD38+ cTfh frequency")
plot(lm(subsetData1$CD27hiCD38hi_..FreqParent ~  subsetData1$cTfh_ICOShiCD38hi_..FreqParent)) 
# Row 14 in subsetData1 has Cook's distance >0.5 and is an outlier with high leverage, will exclude
bivScatter(data1 = subsetData1[-14,], data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD27hiCD38hi_..FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "Yr3: cTfh vs PB response at oneWeek", xLabel= "Plasmablast (freq CD19+IgD-)", yLabel = "ICOS+CD38+ cTfh frequency") + 
  scale_x_continuous(breaks=seq(0,99,1)) + scale_y_continuous(breaks=seq(0,99,2))
# ************** different denominator *********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$CD27hiCD38hi_..FreqParent <- oneWeek$CD27hiCD38hi_..FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000
subsetData1 <- subset(oneWeek, Cohort == "Healthy" )    
subsetData2 <- subset(oneWeek, Cohort == "aPD1"  & cTfh_ICOShiCD38hi_..FreqParent > 1)    #  *****   OUTLIER REMOVED
#  plot(lm(subsetData1$CD27hiCD38hi_..FreqParent ~  subsetData1$cTfh_ICOShiCD38hi_..FreqParent)) 
#  by review of Cook's distance, row 26 in subsetData1 has the highest Cook's distance and will thus exclude
bivScatter(data1 = subsetData1[-26,], data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD27hiCD38hi_..FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "AllYrs: cTfh vs CD27++CD38++ at oneWeek", xLabel= "CD27++CD38++ (freq CD19)", yLabel = "ICOS+CD38+ cTfh frequency") + 
  scale_x_continuous(breaks=seq(0,99,0.05)) + scale_y_continuous(breaks=seq(0,99,2))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_PB-27-38-d7_freqCD19.pdf", device="pdf")
univScatter(data = oneWeek, xData = "CD27hiCD38hi_..FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "AllYrs: cTfh vs CD27++CD38++  at oneWeek", xLabel= "CD27++CD38++ (freq CD19)", yLabel = "ICOS+CD38+ cTfh frequency")
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_PB-27-38-d7_freqCD19_univ.pdf", device="pdf")



## ***************    CD71+ CD20loPB    ********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
subsetData1 <- subset(oneWeek, Cohort == "Healthy" & Year == "3")   
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & Year == "3" )    
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "Yr3: cTfh vs CD20lo response at oneWeek", xLabel= "CD71+CD20lo frequency", yLabel = "ICOS+CD38+ cTfh frequency") + 
  scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1), limits = c(0,15))
univScatter(data = oneWeek, xData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "Yr3: cTfh vs CD20lo at oneWeek", xLabel= "CD20lo (freq CD71+)", yLabel = "ICOS+CD38+ cTfh frequency")
# ************** different denominator *********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  oneWeek$IgDlo_CD71hi_CD20loCD71hi...FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000
subsetData1 <- subset(oneWeek, Cohort == "Healthy" & Year == "3" )   
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & Year == "3" )    
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "Yr3: cTfh vs CD20lo response at oneWeek", xLabel= "PB CD71+CD20lo (freq CD19)", yLabel = "ICOS+CD38+ cTfh frequency")#  + scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1), limits = c(0,15))
univScatter(data = oneWeek, xData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "Yr3: cTfh vs CD20lo at oneWeek", xLabel= "PB CD71+CD20lo (freq CD19)", yLabel = "ICOS+CD38+ cTfh frequency")
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_PB-71hi-20lo-d7_freqCD19_univ.pdf", device="pdf")


## ***************    CD71+ CD20hi ABC    ********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
subsetData1 <- subset(oneWeek, Cohort == "Healthy" & Year == "3" )    
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & Year == "3" ) 
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDlo_CD71hi_ActBCells...FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "Yr3: cTfh vs ABC at oneWeek", xLabel= "ABC CD20hi (freq CD71+)", yLabel = "ICOS+CD38+ cTfh frequency") + 
  scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1), limits = c(0,15))
univScatter(data = oneWeek, xData = "IgDlo_CD71hi_ActBCells...FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "Yr3: cTfh vs CD20hi at oneWeek", xLabel= "ABC CD20hi (freq CD71+)", yLabel = "ICOS+CD38+ cTfh frequency")

# ************** different denominator *********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$IgDlo_CD71hi_ActBCells...FreqParent <-  oneWeek$IgDlo_CD71hi_ActBCells...FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000
subsetData1 <- subset(oneWeek, Cohort == "Healthy" & Year == "3"  )   
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & Year == "3" )      
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDlo_CD71hi_ActBCells...FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "cTfh vs ABC at oneWeek", xLabel= "ABC CD20hi (freq CD19)", yLabel = "ICOS+CD38+ cTfh frequency") #+ scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1), limits=c(0,12))
univScatter(data = oneWeek, xData = "IgDlo_CD71hi_ActBCells...FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "Yr3: cTfh vs CD20hi at oneWeek", xLabel= "ABC CD20hi (freq CD19)", yLabel = "ICOS+CD38+ cTfh frequency")
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ABC-71hi-20hi-d7_freqCD19_univ.pdf", device="pdf")




oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
subsetData1 <- subset(oneWeek, Cohort == "Healthy" )    
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & cTfh_ICOShiCD38hi_..FreqParent > 1)    #  *****   OUTLIER REMOVED
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD27hiCD38hi_..FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "AllYrs: cTfh vs PB response at oneWeek", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency") + scale_x_continuous(breaks=seq(0,99,1)) + scale_y_continuous(breaks=seq(0,99,1))
# ************** different denominator *********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$CD27hiCD38hi_..FreqParent <-  oneWeek$CD27hiCD38hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent /100
subsetData1 <- subset(oneWeek, Cohort == "Healthy" & Year == "3" )    
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & Year == "3" & cTfh_ICOShiCD38hi_..FreqParent > 1)    #  *****   OUTLIER REMOVED
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD27hiCD38hi_..FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "cTfh vs PB response at oneWeek", xLabel= "Plasmablast frequency (freq CD19)", yLabel = "ICOS+CD38+ cTfh frequency") #+ scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1))




oneWeek <- subset(mergedData, TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
subsetData1 <- subset(oneWeek, Cohort == "Healthy")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDlo_CD71hi_ActBCells...FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "AllYrs: cTfh vs ActivB response at oneWeek", xLabel= "CD20hi (freq CD71+)", yLabel = "ICOS+CD38+ cTfh frequency")
univScatter(data = oneWeek, xData = "IgDlo_CD71hi_ActBCells...FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "AllYrs: cTfh vs CD20hi response at oneWeek", xLabel= "CD20hi (freq CD71+)", yLabel = "ICOS+CD38+ cTfh frequency")



cTfhCorrel <- data.frame(CellType = 0, pHC = 0, paPD1 = 0)
z <- array()
data <- subset(oneWeek, Cohort == "Healthy" & Year == "3" )
z[1] <- "CD71+CD20hi_ActivBcells"
z[2] <- round(cor(data[ ,"IgDlo_CD71hi_ActBCells...FreqParent"], data[ ,"cTfh_ICOShiCD38hi_..FreqParent"], method = "pearson", use = "complete.obs"), 2)
z[3] <- round(cor(data[ ,"IgDlo_CD71hi_ActBCells...FreqParent"], data[ ,"cTfh_ICOShiCD38hi_..FreqParent"], method = "pearson", use = "complete.obs"), 2)


## ------------------------------------------------- Averaged frequencies over time  ----------------------------------------------------------


subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "Cohort", title = "cTfh response over time", xLabel = "TimeCategory",
                 yLabel = "ICOS+CD38+ cTfh frequency", groupby = "Subject") + scale_y_continuous(breaks=seq(0,150,5))     # limits = c(0,45)

subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "CD19_CD27.CD38....FreqParent", fillParam = "Cohort", title = "PB response over time", xLabel = "TimeCategory",
                 yLabel = "Plasmablast frequency", groupby = "Subject")  + scale_y_continuous(breaks=seq(0,150,5))     # limits = c(0,45)

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_..FreqParent"))
prePostTimeAveraged(melted, title = "cTfh responses averaged", xLabel = NULL, yLabel = "Average ICOS+CD38+ cTfh frequency") + scale_y_continuous(breaks=seq(0,99,0.5))






subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_..FreqParent"))
prePostTimeAveraged(melted, title = "CD19+ frequencies averaged", xLabel = NULL, yLabel = "CD19+ (freq of CD14-CD4- live)") #+ scale_y_continuous(breaks=seq(0,99,0.25))



subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD27hiCD38hi_..FreqParent"))
prePostTimeAveraged(melted, title = "Plasmablast responses averaged", xLabel = NULL, yLabel = "CD27++CD38++ Plasmablast frequency") + scale_y_continuous(breaks=seq(0,99,0.25))
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$CD27hiCD38hi_..FreqParent <-  subsetData$CD27hiCD38hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent /100
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD27hiCD38hi_..FreqParent"))
prePostTimeAveraged(melted, title = "Plasmablast responses averaged", xLabel = NULL, yLabel = "CD27++CD38++ Plasmablast frequency") + scale_y_continuous(breaks=seq(0,99,0.25))



subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent"))
prePostTimeAveraged(melted, title = "IgDloCD71+CD20lo responses averaged", xLabel = NULL, yLabel = "Average IgDloCD71+CD20lo frequency") + scale_y_continuous(breaks=seq(0,99,3))
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 10000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent"))
prePostTimeAveraged(melted, title = "IgDloCD71+CD20lo responses averaged", xLabel = NULL, yLabel = "Average IgDloCD71+CD20lo frequency") #+ scale_y_continuous(breaks=seq(0,99,3))




subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgDlo_CD71hi_ActBCells...FreqParent"))
prePostTimeAveraged(melted, title = "IgDloCD71+CD20+ ABC responses averaged", xLabel = NULL, yLabel = "Average IgDloCD71+CD20+ frequency") + scale_y_continuous(breaks=seq(0,99,3))
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$IgDlo_CD71hi_ActBCells...FreqParent <-  subsetData$IgDlo_CD71hi_ActBCells...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 10000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgDlo_CD71hi_ActBCells...FreqParent"))
prePostTimeAveraged(melted, title = "IgDloCD71+CD20+ ABC responses averaged", xLabel = NULL, yLabel = "Average IgDloCD71+CD20+ frequency") #+ scale_y_continuous(breaks=seq(0,99,3))



## ------------------------------------------------- Vaccine response determinants and biomarkers ----------------------------------------------------------



subsetData <- subset(mergedData, TimeCategory != "oneWeek" & Cohort != "nonPD1" )
FC_response <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("HAI.titer..H1N1pdm09."))
FC_response$FChai <- FC_response$`late`/FC_response$`baseline`
FC_response2 <- dcast( subset(mergedData, TimeCategory != "late" & cTfh_ICOShiCD38hi_..FreqParent > 0), `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")) 
FC_response2$FCtfh <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject")
bivScatter(data1 = FC_response[which(FC_response$Cohort == "Healthy"),], data2 = FC_response[ which(FC_response$Cohort == "aPD1") ,], 
           name1 = "Healthy", name2 = "aPD1", xData = "FCtfh", yData = "FChai", fillParam = "Cohort", title = "FC HAI vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "H1N1pdm09 titer fold-change") + scale_y_continuous(limits=c(1,128), trans = "log2")



subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("Plasma.CXCL13..pg.mL."))
prePostTimeAveraged(melted, title = "CXCL13 responses averaged", xLabel = NULL, yLabel = "Average plasma CXCL13 (pg/mL)") + scale_y_continuous(breaks=seq(0,99,1))
summary( fit <- aov(value ~ Cohort*TimeCategory, data=melted ) )
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_TimeAveraged.pdf", device="pdf")

subsetData <- subset(melted, TimeCategory == "oneWeek")
rownames(subsetData) <- seq(1:nrow(subsetData))
twoSampleBox(data=subsetData, xData="Cohort", yData="value", fillParam="Cohort", title="CXCL13 at oneWeek", yLabel="plasma CXCL13 (pg/mL)")
# plot(fit <- lm (value ~ Cohort, subsetData))    # row 38 removed for high Cook's distance
twoSampleBox(data=subsetData[-38,], xData="Cohort", yData="value", fillParam="Cohort", title="CXCL13 at oneWeek", yLabel="plasma CXCL13 (pg/mL)")
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_oW.pdf", device="pdf")


subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("HAI.titer..H1N1pdm09."))
prePostTimeAveraged(melted, title = "HAI responses averaged", xLabel = NULL, yLabel = "Average H1N1pdm09 HAI titer") + scale_y_continuous(trans = 'log', breaks=c(2^(2:14)))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAIresponsesTimeAveraged.pdf", device = "pdf")
summary( fit <- aov(value ~ Cohort + TimeCategory, data=melted ) )

twoSampleBox(data = subset(melted, TimeCategory == "late"), xData = "Cohort", yData = "value", fillParam = "Cohort", title = "late HAI titer", yLabel = "H1N1pdm09 HAI titer") +
  scale_y_continuous(trans = "log2")
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAItiter_late_byCohort.pdf", device="pdf")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBox(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="ICOS+CD38+ cTfh at baseline", yLabel="ICOS+CD38+ cTfh freq")
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh_baseline_freq_byCohort.pdf", device="pdf")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late")
twoSampleBox(data=subsetData, xData="Cohort", yData="CD27hiCD38hi_..FreqParent", fillParam="Cohort", title="PB response at Late", yLabel="Plasmablast freq - Late")


## ------------------------------------------------- Scatterplots of HAI titer vs CXCL13  ----------------------------------------------------------


subsetData <- subset(mergedData, Cohort == "aPD1" & TimeCategory == "oneWeek" )
univScatter(subsetData, xData = "HAI.titer..H1N1pdm09.", yData="Plasma.CXCL13..pg.mL.", fillParam = "Cohort", 
                 title = "All yrs aPD1: CXCL13 vs HAI at one week", xLabel= "HAI titer - H1N1pdm09", yLabel = "Plasma CXCL13 (pg/mL)")


subsetData <- subset(mergedData, Cohort != "nonPD1" )
crossTimeCategory <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c('Plasma.CXCL13..pg.mL.' , 'HAI.titer..H1N1pdm09.') )
crossTimeCategory2 <- dcast(crossTimeCategory, Subject + Cohort ~ TimeCategory + variable)
univScatter(crossTimeCategory2, xData = "baseline_Plasma.CXCL13..pg.mL.", yData="late_HAI.titer..H1N1pdm09.", fillParam = NULL, 
            title = "All yrs: Late HAI vs baseline CXCL13", xLabel= "baseline Plasma CXCL13 (pg/mL)", yLabel = "Late HAI titer - H1N1pdm09")   # nonstd appearance due to lack of fillParam
univScatter(crossTimeCategory2, xData = "oneWeek_Plasma.CXCL13..pg.mL.", yData="late_HAI.titer..H1N1pdm09.", fillParam = NULL, 
            title = "All yrs: Late HAI vs d7 CXCL13", xLabel= "oneWeek Plasma CXCL13 (pg/mL)", yLabel = "Late HAI titer - H1N1pdm09")   # nonstd appearance due to lack of fillParam



 
subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "HAI.titer..H1N1pdm09.", fillParam = "Cohort", title = "HAI titers over time", xLabel = "TimeCategory",
            yLabel = "HAI titer", groupby = "Subject"); a + scale_y_continuous(trans='log2')

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBox(data=subsetData, xData="Cohort", yData="HAI.titer..H1N1pdm09.", fillParam="Cohort", title="All yrs: HAI titers at d0", yLabel="H1N1pdm09 titer") + scale_y_continuous(trans='log2', breaks=c(2^(2:14) ) )
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBox(data=subsetData, xData="Cohort", yData="HAI.titer..H1N1pdm09.", fillParam="Cohort", title="All yrs: HAI titers at d7", yLabel="H1N1pdm09 titer") + scale_y_continuous(trans='log2', breaks=c(2^(2:14) ) )
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late")
twoSampleBox(data=subsetData, xData="Cohort", yData="HAI.titer..H1N1pdm09.", fillParam="Cohort", title="All yrs: HAI titers at d21-28", yLabel="H1N1pdm09 titer") + scale_y_continuous(trans='log2' , breaks=c(2^(2:14) ) )

subsetData <- subset(mergedData, Cohort != "nonPD1" )
seroconv <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("HAI.titer..H1N1pdm09.")); 
seroconv$FC <- seroconv$`late`/seroconv$`baseline`
twoSampleBox(data=seroconv, xData="Cohort", yData="FC", fillParam="Cohort", title="All yrs: seroconversion", yLabel="Seroconversion factor") + scale_y_continuous(breaks=seq(0,99,4))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/SeroconversionFactor.pdf", device = "pdf")



subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "Plasma.CXCL13..pg.mL.", fillParam = "Cohort", title = "CXCL13 over itme", xLabel = "TimeCategory",
            yLabel = "HAI titer", groupby = "Subject") + scale_y_continuous(breaks=seq(0,150,5))     # limits = c(0,45)


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBox(data=subsetData, xData="Cohort", yData="Plasma.CXCL13..pg.mL.", fillParam="Cohort", title="All yrs: Plasma CXCL13 at d0", yLabel="Plasma CXCL13 (pg/mL)")
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBox(data=subsetData, xData="Cohort", yData="Plasma.CXCL13..pg.mL.", fillParam="Cohort", title="All yrs: Plasma CXCL13 at d7", yLabel="Plasma CXCL13 (pg/mL)")
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late")
twoSampleBox(data=subsetData, xData="Cohort", yData="Plasma.CXCL13..pg.mL.", fillParam="Cohort", title="All yrs: Plasma CXCL13 at d21-28", yLabel="Plasma CXCL13 (pg/mL)")

summary(lm( Plasma.CXCL13..pg.mL. ~ Cohort + Year + TimeCategory, data=mergedData))





## ------------------------------------------------- Working with the medianFI for different parameters  ----------------------------------------------------------


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBox(data=subsetData, xData="Cohort", yData="Ki67hiCD38hicTfh_medfi.Blimp1.", fillParam="Cohort", title="All yrs: Blimp1 in +/+ at d7", yLabel="Blimp1 median FI")


## ------------------------------------------------- additional analyses  ----------------------------------------------------------


# subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
# a <- colnames(subsetData[ , sapply(subsetData, class) == 'character'])
# a <- cor(subsetData[ , - grep(paste(c(a, "Cohort"), collapse="|"), colnames(subsetData))], use="pairwise.complete.obs")
# a[order(a[,20],decreasing = F),20]


oneWeek <- subset(mergedData, TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
subsetData1 <- subset(oneWeek, Cohort == "Healthy")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDlo_CD71hi_CD20loCD71hi.CD32hi...FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "AllYrs: cTfh vs PB CD32hi at oneWeek", xLabel= "CD20lo - CD32hi frequency", yLabel = "ICOS+CD38+ cTfh frequency")
univScatter(data = oneWeek, xData = "IgDlo_CD71hi_CD20loCD71hi.CD32hi...FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "AllYrs: cTfh vs PB CD32hi at oneWeek", xLabel= "CD20lo - CD32hi frequency", yLabel = "ICOS+CD38+ cTfh frequency")
  

oneWeek <- subset(mergedData, TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
subsetData1 <- subset(oneWeek, Cohort == "Healthy")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDlo_CD71hi_CD20loCD71hi.CD86....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "AllYrs: cTfh vs PB CD86 resp at oneWeek", xLabel= "CD20lo - CD86hi frequency", yLabel = "ICOS+CD38+ cTfh frequency")
univScatter(data = oneWeek, xData = "IgDlo_CD71hi_CD20loCD71hi.CD86....FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "AllYrs: cTfh vs CD20lo CD86+ at oneWeek", xLabel= "CD20lo - CD86hi frequency", yLabel = "ICOS+CD38+ cTfh frequency")


# as predicted by the elastic net regression
oneWeek <- subset(mergedData, TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
subsetData1 <- subset(oneWeek, Cohort == "Healthy")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDloCD71hi..medfi.CD86.", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "AllYrs: cTfh vs CD86 resp at oneWeek", xLabel= "IgDloCD71hi - CD86 medFI", yLabel = "ICOS+CD38+ cTfh frequency")
univScatter(data = oneWeek, xData = "IgDloCD71hi..medfi.CD86.", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "TimeCategory", title = "AllYrs: cTfh vs CD86 resp at oneWeek", xLabel= "IgDloCD71hi - CD86 medFI", yLabel = "ICOS+CD38+ cTfh frequency")












