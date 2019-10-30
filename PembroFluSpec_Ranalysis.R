library("grid")
library("ggplot2")
library("gplots")
library("viridis")
library("reshape2")
# library("ggpubr")
# library("ggcorrplot")
library("gridExtra")
library("RColorBrewer")
library("scales")
library("MASS")

source('D:/Pembro-Fluvac/Analysis/PembroFluSpec_Ranalysis_files/PembroFluSpec_PlottingFunctions.R')


## ------------------------------------------------- merge & qc data from spreadsheets  ----------------------------------------------------------

setwd("D:/Pembro-Fluvac/Analysis")
mergedData <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_freqParent_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
fit <- lm(data=mergedData[ , - which( colnames(mergedData) == "Label" | colnames(mergedData) == "Subject" | colnames(mergedData) == "Cohort" | 
                                        colnames(mergedData) == "TimeCategory" | colnames(mergedData) == "TimePoint")], 
          formula = TflowBatch  ~ .)
a <- as.data.frame(summary(fit)$coefficients)
a <- rownames(a[which(a$`Pr(>|t|)` < 0.05), ] )
if (length (a) > 0) { colnames(mergedData) [ grep( paste(a, collapse="|"), colnames(mergedData)) ] <- paste0( colnames(mergedData)[ grep( paste(a, collapse="|"), colnames(mergedData)) ], ".batchEffect" ) }

serology <- read.csv(file= "D:/Pembro-Fluvac/Analysis/mergedData/SerologyMeasurements.csv", stringsAsFactors = F, header = T)

temp1 <- merge(x = mergedData, y=serology, all = T, suffixes = c(".Tflow",".Serology"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year'))


Bflow <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_freqParent_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
Bflow <- Bflow[, -c( grep( colnames(Bflow), pattern=paste0(c("Naive","Lymphocytes", "CD4","CD19hi_..FreqParent","CD19hi_CD11c","CD19hi_CD16","CD19hi_CD23","CD1hi9_CD25",
                                                             "CD19hi_CD32","CD19hi_CD38lo...FreqParent","CD19hi_CD80","CD19hi_CD86","CD19hi_IgM","CD19hi_CD138",
                                                             "CD19hi_CXCR4","CD19hi_CD71","CD19hi_Ki67..CXCR4","CD19hi_Ki67..Ki67","CD38lo...FreqParent",
                                                             "CD19hi_Tbet"), collapse ="|")) )]    # cut out columns to test for batch fx, else variables > observations
fit <- lm(data=Bflow[ , - which( colnames(Bflow) == "Label" | colnames(Bflow) == "Subject" | colnames(Bflow) == "Cohort" | 
                                        colnames(Bflow) == "TimeCategory" | colnames(Bflow) == "TimePoint")], 
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
mergedData <- temp3 


mergedData$TimeCategory <- factor(mergedData$TimeCategory, levels = c("baseline", "oneWeek","late"))
mergedData$Cohort <- factor(mergedData$Cohort, levels = c("Healthy", "aPD1","nonPD1"))
# rm(Tflow); rm(Bflow); rm(temp1); rm(temp2); rm(temp3)







## ------------------------------------------------- PB and Tfh just day 7 frequencies ----------------------------------------------------------

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek" & Year == "3" & cTfh_ICOShiCD38hi_..FreqParent > 0)
a <- twoSampleBox(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="Yr3: ICOS+CD38+ cTfh at d7", yLabel="ICOS+CD38+ cTfh frequency at d7")
a + scale_y_continuous(breaks=seq(1:12), limits=c(0.5,11))

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek" & Year == "3")
a <- twoSampleBox(data=subsetData, xData="Cohort", yData="CD19_CD27.CD38....FreqParent", fillParam="Cohort", title="Yr3: Plasmablast at d7", yLabel="CD19+CD27++CD38++ frequency at d7")
a + scale_y_continuous(breaks=seq(0,99,1))

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 0)
a <- twoSampleBox(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: ICOS+CD38+ cTfh at d7", yLabel="ICOS+CD38+ cTfh frequency at d7")
a + scale_y_continuous(breaks=seq(0,99,1))

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
a <- twoSampleBox(data=subsetData, xData="Cohort", yData="CD19_CD27.CD38....FreqParent", fillParam="Cohort", title="All yrs: Plasmablast at d7", yLabel="CD19+CD27++CD38++ frequency at d7")
a # + scale_y_continuous(breaks=seq(0,99,1))


## ------------------------------------------------- Tfh and PB Fold-change at day 7  ----------------------------------------------------------

subsetData <- subset(mergedData, TimeCategory != "late" & Year == "3" & Cohort != "nonPD1")
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
a <- twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="Yr3: cTfh response at one week", yLabel="ICOS+CD38+ cTfh fold-change")
a + scale_y_continuous(breaks=seq(0,99,1))

subsetData <- subset(mergedData, TimeCategory != "late" & Year == "3" & Cohort != "nonPD1")
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD19_CD27.CD38....FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
a <- twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="Yr3: Plasmablast response at one week", yLabel="CD19+CD27++CD38++ fold-change")
a + scale_y_continuous(breaks=seq(0,99,1))

subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & cTfh_ICOShiCD38hi_..FreqParent > 0)
FC_Tfhresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`oneWeek`/FC_Tfhresponse$`baseline`
a <- twoSampleBox(data=FC_Tfhresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="All yrs: cTfh response at one week", yLabel="ICOS+CD38+ cTfh fold-change")
a + scale_y_continuous(breaks=seq(0,99,1))



## ------------------------------------------------- PB vs cTfh response correlations ----------------------------------------------------------

subsetData <- subset(mergedData, Cohort == "aPD1" & Year == "3" )
subsubsetData <- subset(subsetData, TimeCategory == "oneWeek")
# xData = "CD19_CD27.CD38....FreqParent"; yData="cTfh_ICOShiCD38hi_..FreqParent"
# z <- summary(rlm(subsubsetData[, yData] ~ subsubsetData[, xData]))
# y <- summary(lm(subsubsetData[, yData] ~ subsubsetData[, xData]))
a <- univScatter(subsubsetData, xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam = "Cohort", 
                      title = "Yr3 aPD1: cTfh vs PB at one week", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")
a + scale_x_continuous(breaks=seq(0,99,0.25)) + scale_y_continuous(breaks=seq(0,99,1))

subsetData <- subset(mergedData, Cohort == "Healthy" & Year == "3" & cTfh_ICOShiCD38hi_..FreqParent > 1)    #  *****   OUTLIER REMOVED
subsubsetData <- subset(subsetData, TimeCategory == "oneWeek")
a <- univScatter(subsubsetData, xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam = "Cohort", 
                      title = "Yr3 Healthy: cTfh vs PB at one week", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")
a + scale_x_continuous(breaks=seq(0,99,0.25)) + scale_y_continuous(breaks=seq(0,99,1))

oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
subsetData1 <- subset(oneWeek, Cohort == "Healthy" & Year == "3" & cTfh_ICOShiCD38hi_..FreqParent > 1)    #  *****   OUTLIER REMOVED
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & Year == "3" & cTfh_ICOShiCD38hi_..FreqParent > 1)    #  *****   OUTLIER REMOVED
a <- bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "Yr3: cTfh vs PB response at oneWeek", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")
a + scale_x_continuous(breaks=seq(0,99,0.25)) + scale_y_continuous(breaks=seq(0,99,1))



subsetData <- subset(mergedData, Cohort == "aPD1" & TimeCategory == "oneWeek" )
univScatter(subsetData, xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam = "Cohort", 
            title = "All yrs aPD1: cTfh vs PB at one week", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")


subsetData <- subset(mergedData, Cohort == "Healthy" & TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 1)   #  *****   OUTLIER REMOVED
a <- univScatter(subsetData, xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam = "Cohort", 
            title = "All yrs Healthy: cTfh vs PB at one week", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")
a + scale_x_continuous(breaks=seq(0,99,0.25)) + scale_y_continuous(breaks=seq(0,99,1))



subsetData <- subset(mergedData, Cohort != "nonPD1" )
a <- prePostTime(subsetData, xData = "TimeCategory", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "Cohort", title = "cTfh response over time", xLabel = "TimeCategory",
                 yLabel = "ICOS+CD38+ cTfh frequency", groupby = "Subject"); a + scale_y_continuous(breaks=seq(0,150,5))     # limits = c(0,45)

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_..FreqParent"))
a <- prePostTimeAveraged(melted, title = "cTfh responses averaged", xLabel = NULL, yLabel = "Average ICOS+CD38+ cTfh frequency"); a + scale_y_continuous(breaks=seq(0,99,0.5))
  
subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19_CD27.CD38....FreqParent"))
a <- prePostTimeAveraged(melted, title = "Plasmablast responses averaged", xLabel = NULL, yLabel = "Average Plasmablast frequency"); a + scale_y_continuous(breaks=seq(0,99,0.25))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("Plasma.CXCL13..pg.mL."))
a <- prePostTimeAveraged(melted, title = "CXCL13 responses averaged", xLabel = NULL, yLabel = "Average plasma CXCL13 (pg/mL)"); a + scale_y_continuous(breaks=seq(0,99,1))



subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late")
twoSampleBox(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="cTfh response at Late", yLabel="ICOS+CD38+ cTfh freq - Late")

subsetData <- subset(mergedData, Cohort != "nonPD1" )
a <- prePostTime(subsetData, xData = "TimeCategory", yData = "CD19_CD27.CD38....FreqParent", fillParam = "Cohort", title = "PB response over time", xLabel = "TimeCategory",
                 yLabel = "Plasmablast frequency", groupby = "Subject"); a + scale_y_continuous(breaks=seq(0,150,5))     # limits = c(0,45)
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late")
twoSampleBox(data=subsetData, xData="Cohort", yData="CD19_CD27.CD38....FreqParent", fillParam="Cohort", title="PB response at Late", yLabel="Plasmablast freq - Late")


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
a <- prePostTime(subsetData, xData = "TimeCategory", yData = "HAI.titer..H1N1pdm09.", fillParam = "Cohort", title = "HAI titers over time", xLabel = "TimeCategory",
            yLabel = "HAI titer", groupby = "Subject"); a + scale_y_continuous(trans='log2')

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
a <- twoSampleBox(data=subsetData, xData="Cohort", yData="HAI.titer..H1N1pdm09.", fillParam="Cohort", title="All yrs: HAI titers at d0", yLabel="H1N1pdm09 titer"); a + scale_y_continuous(trans='log2', breaks=c(2^(2:14) ) )
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
a <- twoSampleBox(data=subsetData, xData="Cohort", yData="HAI.titer..H1N1pdm09.", fillParam="Cohort", title="All yrs: HAI titers at d7", yLabel="H1N1pdm09 titer"); a + scale_y_continuous(trans='log2', breaks=c(2^(2:14) ) )
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late")
a <- twoSampleBox(data=subsetData, xData="Cohort", yData="HAI.titer..H1N1pdm09.", fillParam="Cohort", title="All yrs: HAI titers at d21-28", yLabel="H1N1pdm09 titer"); a + scale_y_continuous(trans='log2' , breaks=c(2^(2:14) ) )

subsetData <- subset(mergedData, Cohort != "nonPD1" )
seroconv <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("HAI.titer..H1N1pdm09.")); 
seroconv$FC <- seroconv$`late`/seroconv$`baseline`
a <- twoSampleBox(data=seroconv, xData="Cohort", yData="FC", fillParam="Cohort", title="All yrs: seroconversion", yLabel="Seroconversion factor")
a + scale_y_continuous(breaks=seq(0,99,4))




subsetData <- subset(mergedData, Cohort != "nonPD1" )
a <- prePostTime(subsetData, xData = "TimeCategory", yData = "Plasma.CXCL13..pg.mL.", fillParam = "Cohort", title = "CXCL13 over itme", xLabel = "TimeCategory",
            yLabel = "HAI titer", groupby = "Subject"); a + scale_y_continuous(breaks=seq(0,150,5))     # limits = c(0,45)


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBox(data=subsetData, xData="Cohort", yData="Plasma.CXCL13..pg.mL.", fillParam="Cohort", title="All yrs: Plasma CXCL13 at d0", yLabel="Plasma CXCL13 (pg/mL)")
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBox(data=subsetData, xData="Cohort", yData="Plasma.CXCL13..pg.mL.", fillParam="Cohort", title="All yrs: Plasma CXCL13 at d7", yLabel="Plasma CXCL13 (pg/mL)")
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late")
twoSampleBox(data=subsetData, xData="Cohort", yData="Plasma.CXCL13..pg.mL.", fillParam="Cohort", title="All yrs: Plasma CXCL13 at d21-28", yLabel="Plasma CXCL13 (pg/mL)")




## ------------------------------------------------- Working with the medianFI for different parameters  ----------------------------------------------------------














