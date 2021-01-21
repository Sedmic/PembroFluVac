#+ fig.width=6, fig.height=6
#' ##--------- call libraries --
library("grid")
library("ggplot2")
library("gplots")
library("viridis")
library("reshape2")
library("gridExtra")
library("RColorBrewer")
library("scales")
library("stringr")
library("rstatix")
library("dplyr")
library("tidyr")
library("gridExtra")
source('D:/Pembro-Fluvac/Analysis/Ranalysis/PembroFluSpec_PlottingFunctions.R')

#' ## ----------- merge & qc data from spreadsheets  --------------------
#'
#'
setwd("D:/Pembro-Fluvac/Analysis")
mergedData <- Tflow <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_freqParent_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
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

#' ## ----------- calculate fold-changes for HAI, Tfh, and PB --------------------

mergedData$logHAI <- log(mergedData$H1N1pdm09.HAI.titer)
subsetData <- subset(mergedData, TimeCategory != "oneWeek" & Cohort != "nonPD1" )
FC_response <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("logHAI"))
FC_response$FChai_late <- FC_response$`late`/FC_response$`baseline`

subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" )
FC_response2 <- dcast( subset(subsetData, cTfh_ICOShiCD38hi_..FreqParent > 0), `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")) 
FC_response2$FCtfh_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject")
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD27hiCD38hi_..FreqParent")) 
FC_response2$FCPB_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject")
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("Plasma.CXCL13..pg.mL.")) 
FC_response2$FCCXCL13_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject")
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgG1_Total.sialylated")) 
FC_response2$IgG1sial_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject")

FC_response <- FC_response[, c("Subject","Cohort","FChai_late","FCtfh_oW", "FCPB_oW", "FCCXCL13_oW", "IgG1sial_oW")]
mergedData <- merge(x = mergedData, y= FC_response, all=T, by = c('Subject', 'Cohort'))

mergedData$ASC_ABCratio <- mergedData$IgDlo_CD71hi_CD20loCD71hi...FreqParent / mergedData$IgDlo_CD71hi_ActBCells...FreqParent

# write.csv(mergedData, file = "D:/Pembro-Fluvac/Analysis/mergedData/allMergedData.csv")


# mergedData <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/allMergedData.csv", stringsAsFactors = F, header = T, row.names = 1)

mergedData$TimeCategory <- factor(mergedData$TimeCategory, levels = c("baseline", "oneWeek","late"))
mergedData$Cohort <- factor(mergedData$Cohort, levels = c("Healthy", "aPD1","nonPD1"))
mergedData$dummy <- "dummy"

#' ## ----------- Dataset overview --------------------

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
         ); dev.off()

#' ## ----------- Cohort description --------------------
#'
#'
 
table(mergedData$Cohort, mergedData$Sex, mergedData$Year)  
table(serology$Cohort, serology$TimeCategory, serology$Year, !is.na(serology$H1N1pdm09.HAI.titer))
ggplot( data = subset(mergedData, Cohort == "aPD1" & TimeCategory == "baseline"), aes(x=Cycle.of.Immunotherapy)) + geom_histogram(color="white",fill="#FFB18C") + 
  theme_bw() + theme(axis.text=element_text(color="black",size=18), axis.title=element_text(color="black",size=24)) + scale_x_continuous(breaks=seq(0,30,5)) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CyclesImmunotherapy.pdf", device = 'pdf', height=3)

median(mergedData[which(mergedData$Cohort == "aPD1" & mergedData$TimeCategory == "baseline"), 'Cycle.of.Immunotherapy'], na.rm = T)  
paste("Median Age aPD1: ", median(mergedData[which(mergedData$Cohort == "aPD1" & mergedData$TimeCategory == "baseline"), 'Age'], na.rm = T)  )
paste("Median Age Healthy: ", median(mergedData[which(mergedData$Cohort == "Healthy" & mergedData$TimeCategory == "baseline"), 'Age'], na.rm = T)   )
paste("Range Age aPD1: ", range(mergedData[which(mergedData$Cohort == "aPD1" & mergedData$TimeCategory == "baseline"), 'Age'], na.rm = T)  )
paste("Range Age Healthy: ", range(mergedData[which(mergedData$Cohort == "Healthy" & mergedData$TimeCategory == "baseline"), 'Age'], na.rm = T)   )



#' ## ----------- Global overview of phenotypic differences --------------------
#'
#'

temp <- as.data.frame( cor(mergedData[ ,sapply(mergedData, is.numeric)], mergedData$Cycle.of.Immunotherapy, use="pairwise.complete.obs"))
temp$names <- rownames(temp); rownames(temp) <- NULL; 
temp[order(temp$V1,decreasing = T),][1:20,];   temp[order(temp$V1,decreasing = F),][1:20,]; 

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline");  
res <- lapply(subsetData[,c(5:80,82:161,163:292,300:450)], function(x) t.test(x ~ subsetData$Cohort, var.equal = TRUE))
res <- sapply(res, "[[", "p.value")



subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_NaiveB...FreqParent"))
prePostTimeAveraged(melted, title = "Naive B frequencies averaged", xLabel = NULL, yLabel = "NaiveB (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))
twoSampleBar(data=subset(subsetData, TimeCategory == "baseline"), xData="Cohort", yData="CD19hi_NaiveB...FreqParent", fillParam="Cohort", title="NaiveB in CD19 at baseline", 
             yLabel="IgD+CD27- (% CD19)", position="left")  
univScatter(data = subset(subsetData, TimeCategory == "baseline"), yData = "CD19hi_NaiveB...FreqParent", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Cycle of IT vs NaiveB at baseline", yLabel= "IgD+CD27- (% CD19)", xLabel = "Cycle of immunotherapy") #+ scale_y_continuous(limits= c(0,12))
summary(lm( data = mergedData, CD19hi_NaiveB...FreqParent ~ Cohort + Year + TimeCategory ))
summary(lm( data = subset(subsetData, TimeCategory == "baseline"), CD19hi_NaiveB...FreqParent ~ Cohort + Year + Age))


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
twoSampleBar(data=subset(subsetData, TimeCategory == "baseline"), xData="Cohort", yData="CD19hi_CD80....FreqParent", fillParam="Cohort", title="CD80 in CD19 at baseline", 
             yLabel="CD80+ (% CD19)", position="right")  
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CD19CD80_atBaseline_barplot.pdf")
univScatter(data = subset(subsetData, TimeCategory == "baseline"), yData = "CD19hi_CD80....FreqParent", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Cycle of IT vs CD19+CD80+ at baseline", yLabel= "CD80+ (% CD19)", xLabel = "Cycle of immunotherapy") + scale_y_continuous(limits= c(0,8))
summary(lm( data = subsetData, CD19hi_CD80....FreqParent ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_CD86....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+CD86+ frequencies averaged", xLabel = NULL, yLabel = "CD19+CD86+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))
twoSampleBar(data=subset(subsetData, TimeCategory == "baseline"), xData="Cohort", yData="CD19hi_CD86....FreqParent", fillParam="Cohort", title="CD86 in CD19 at baseline", 
             yLabel="CD86+ (% CD19)", position="left") 
summary(lm( data = subsetData, CD19hi_CD86....FreqParent ~ Cohort + Year + TimeCategory))
summary(lm( data = subset(subsetData, TimeCategory == "baseline"), logCD19CD86_FreqParent ~ Cohort ))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_IgM....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+IgM+ frequencies averaged", xLabel = NULL, yLabel = "CD19+IgM+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_CD138....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+CD138+ frequencies averaged", xLabel = NULL, yLabel = "CD19+CD138+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_Tbet....FreqParent"))
prePostTimeAveraged(melted, title = "CD19+Tbet+ frequencies averaged", xLabel = NULL, yLabel = "CD19+Tbet+ (freq CD19)") # + scale_y_continuous(breaks=seq(0,99,0.25))


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="CD27hiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: PB at d0", yLabel="CD19+IgD-CD27++CD38++ frequency at d0")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline" & IgDlo_CD71hi_ActBCells...FreqParent > 0 & IgDlo_CD71hi_CD20loCD71hi...FreqParent >0)
twoSampleBar(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_ActBCells...FreqParent", fillParam="Cohort", title="baseline ABC", 
             yLabel="CD20+ ABC (% CD71+)") + scale_y_continuous(limits=c(0,100), breaks=seq(0,100,10))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/ABC_atBaseline_barPlot.pdf", width=4)
twoSampleBar(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_CD20loCD71hi...FreqParent", fillParam="Cohort", title="baseline ASC", 
             yLabel="CD20- ASC (% CD71+)")+ scale_y_continuous(limits=c(0,100), breaks=seq(0,100,10))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/ASC_atBaseline_barPlot.pdf", width=4)

twoSampleBar(data=subsetData, xData="Cohort", yData="ASC_ABCratio", fillParam="Cohort", title="baseline ASC/ABC ratio", 
             yLabel="ASC/ABC ratio") + scale_y_continuous(limits=c(0,8), breaks=seq(0,100,1))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/ASC_ABCratio_atBaseline_barPlot.pdf", width=4)



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



#' ## ----------- Averaged frequencies over time  --------------------
#'
#'

subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(data = subsetData, xData = "TimeCategory", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "Cohort", title = "cTfh response", xLabel = "TimeCategory",
            yLabel = "ICOS+CD38+ (% cTfh)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,150,2), limits=c(0,18))     
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_PerSubject.pdf")
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_..FreqParent"))
melted <- melted[which(melted$TimeCategory != 'late'),]
prePostTimeAveraged(data = melted, title = "cTfh", xLabel = NULL, yLabel = "ICOS+CD38+ (% cTfh)") + scale_y_continuous(breaks=seq(1,10,0.5), limits=c(3.0, 7))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_freqcTfh.pdf", width=4)
summary( fit <- lm(formula = value ~ Cohort * TimeCategory, data = melted) );  tukey_hsd(fit)       # testing using unpaired two-way ANOVA with Tukey's posthoc
melted %>% group_by(Cohort) %>% kruskal_test(formula = value ~ TimeCategory)                        # testing per cohort using Kruskall-Wallis without adjustment
index <- melted[melted$TimeCategory=="oneWeek", 'Subject']
melted.paired <- melted[which(melted$Subject %in% index), ]
prePostTimeAveraged(data = melted.paired, title = "cTfh", xLabel = NULL, yLabel = "ICOS+CD38+ (% cTfh)") + scale_y_continuous(breaks=seq(1,10,0.5), limits=c(3.0, 7))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_freqcTfh_paired.pdf", width=4)
# write.csv(melted.paired, file = "cTfh_responses_freqParent.csv")                                  # paired two-way ANOVA results calculated in Prism

#************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$cTfh_ICOShiCD38hi_..FreqParent <- subsetData$cTfh_ICOShiCD38hi_..FreqParent * subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /10000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_..FreqParent"))
melted <- melted[which(melted$TimeCategory != 'late'),]
prePostTimeAveraged(melted, title = "cTfh", xLabel = NULL, yLabel = "ICOS+CD38+ cTfh (% CD4)") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0.25,1))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_freqCD4.pdf", width=4)
summary( fit <- lm(formula = value ~ Cohort * TimeCategory, data = melted) ); tukey_hsd(fit)        # testing using unpaired two-way ANOVA with Tukey's posthoc
melted %>% group_by(Cohort) %>% kruskal_test(formula = value ~ TimeCategory)                        # testing per cohort using Kruskall-Wallis without adjustment

index <- melted[melted$TimeCategory=="oneWeek", 'Subject']
melted.paired <- melted[which(melted$Subject %in% index), ]
prePostTimeAveraged(data = melted.paired, title = "cTfh", xLabel = NULL, yLabel = "ICOS+CD38+ cTfh (% CD4+)") + 
  scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0.25,1))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_freqCD4_paired.pdf", width=4)
# write.csv(melted.paired, file = "cTfh_responses_freqCD4.csv")








subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD19hi_..FreqParent"))
prePostTimeAveraged(melted, title = "CD19+ frequencies averaged", xLabel = NULL, yLabel = "CD19+ (freq of CD14-CD4- live)") #+ scale_y_continuous(breaks=seq(0,99,0.25))


subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "CD19_CD27.CD38....FreqParent", fillParam = "Cohort", title = "PB response over time", xLabel = "TimeCategory",
            yLabel = "Plasmablast frequency", groupby = "Subject")  + scale_y_continuous(breaks=seq(0,150,5))     # limits = c(0,45)



subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD27hiCD38hi_..FreqParent"))
prePostTimeAveraged(melted, title = "Plasmablast responses averaged", xLabel = NULL, yLabel = "CD27++CD38++ Plasmablast (% nonnaive B)") + scale_y_continuous(breaks=seq(0,99,0.25))
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$CD27hiCD38hi_..FreqParent <-  subsetData$CD27hiCD38hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent /100
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("CD27hiCD38hi_..FreqParent"))
prePostTimeAveraged(melted, title = "Plasmablast responses averaged", xLabel = NULL, yLabel = "CD27++CD38++ Plasmablast (% CD19)") + scale_y_continuous(breaks=seq(0,99,0.25))




subsetData <- subset(mergedData, Cohort != "nonPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent != 0 )   # exclude zero values due to presumptive technical artifact
prePostTime(subsetData, xData = "TimeCategory", yData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", fillParam = "Cohort", title = "ASC response by subject", 
            xLabel = "TimeCategory", yLabel = "CD20- ASC (% IgD-CD71+)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,150,5))     
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/ASCResp_overTime_PerSubject.pdf", width = 8)
subsetData <- subset(mergedData, Cohort != "nonPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent != 0 )   # exclude zero values due to presumptive technical artifact
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent"))
prePostTimeAveraged(data = subset(melted, TimeCategory != "late"), title = "ASC", xLabel = NULL, yLabel = "CD20- ASC (% IgD-CD71+)") + 
  scale_y_continuous(breaks=seq(0,99,5), limits = c(30,70))
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ASCResp_overTime_freqCD71.pdf", width=4)
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent != 0 )   # exclude zero values due to presumptive technical artifact
subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 10000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent"))
prePostTimeAveraged(data = subset(melted, TimeCategory != "late"), title = "ASC", xLabel = NULL, yLabel = "CD20- ASC (% CD19)") + 
  scale_y_continuous(breaks=seq(0,99,0.25))
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ASCResp_overTime_freqCD19.pdf", width=4)




subsetData <- subset(mergedData, Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent != 0 )   # exclude zero values due to presumptive technical artifact
prePostTime(subsetData, xData = "TimeCategory", yData = "IgDlo_CD71hi_ActBCells...FreqParent", fillParam = "Cohort", title = "ABC response by subject", 
            xLabel = "TimeCategory", yLabel = "CD20+ ABC (% IgD-CD71+)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,150,5)) 
subsetData <- subset(mergedData, Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent != 0 )   # exclude zero values due to presumptive technical artifact
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgDlo_CD71hi_ActBCells...FreqParent"))
prePostTimeAveraged(data = subset(melted, TimeCategory != "late"), title = "ABC", xLabel = NULL, 
                    yLabel = "CD20+ ABC (% CD71+)") + scale_y_continuous(breaks=seq(0,99,5), limits = c(20,60))
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ABCResp_overTime_freqCD71.pdf", width=4)
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent != 0 )   # exclude zero values due to presumptive technical artifact
subsetData$IgDlo_CD71hi_ActBCells...FreqParent <-  subsetData$IgDlo_CD71hi_ActBCells...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 10000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("IgDlo_CD71hi_ActBCells...FreqParent"))
prePostTimeAveraged(data = subset(melted, TimeCategory != "late"), title = "ABC", xLabel = NULL, 
                    yLabel = "CD20+ ABC (% CD19)")  + scale_y_continuous(breaks=seq(0,99,0.2), limits = c(0.3,1.5))
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ABCResp_overTime_freqCD19.pdf", width=4)





# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$cTfh_ICOShiCD38hi_CTLA4hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_CTLA4hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent * 
  subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /1000000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CTLA4hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ CTLA4hi frequencies averaged", xLabel = NULL, yLabel = "+/+ CTLA4hi (freq CD4)") 
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_CTLA4hi...FreqParent ~ Cohort + Year + TimeCategory))


subsetData <- subset(mergedData, Cohort != "nonPD1")
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_Ki67hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ Ki67hi frequencies averaged", xLabel = NULL, yLabel = "Ki67hi (% cTfh)") 
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_Ki67hi...FreqParent ~ Cohort  + TimeCategory + Year))    # Year is a batch effect
tukey_hsd(aov(data = subsetData, formula = cTfh_ICOShiCD38hi_Ki67hi...FreqParent ~ Cohort:TimeCategory))  # proportion of Parent
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBar(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_Ki67hi...FreqParent", fillParam="Cohort", title="Ki67 in ICOS+CD38+ cTfh at oW", 
             yLabel="+/+ Ki67hi (% cTfh) ")
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1")
subsetData$cTfh_ICOShiCD38hi_Ki67hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_Ki67hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent * 
  subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /1000000
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_Ki67hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ Ki67hi frequencies averaged", xLabel = NULL, yLabel = "+/+ Ki67hi (freq CD4)") 
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_Ki67hi...FreqParent ~ Cohort  + TimeCategory + Year))    # Year is a batch effect
tukey_hsd(aov(data = subsetData, formula = cTfh_ICOShiCD38hi_Ki67hi...FreqParent ~ Cohort:TimeCategory))  # proportion of CD4 not FreqParent

subsetData$cTfh_ICOShiCD38hi_Ki67hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_Ki67hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent * 
  subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /1000000
twoSampleBar(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_Ki67hi...FreqParent", fillParam="Cohort", title="Ki67 in ICOS+CD38+ cTfh at oW", 
             yLabel="+/+ Ki67hi (freq CD4) ")   # this is really just reflective of the +/+ frequency difference between cohorts
twoSampleBar(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_cTfh..._medfi.Ki67.", fillParam="Cohort", title="Ki67 in ICOS+CD38+ cTfh at oW", 
             yLabel="+/+ for medianFI of Ki67")  



subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_cTfh..._medfi.CD25."))
prePostTimeAveraged(melted, title = "+/+ medFI CD25", xLabel = NULL, yLabel = "medianFI CD25 in ICOS+CD38+ cTfh") # + scale_y_continuous(breaks=seq(0,99,0.25))
twoSampleBar(data=subset(subsetData, TimeCategory == "baseline"), xData="Cohort", yData="cTfh_ICOShiCD38hi_cTfh..._medfi.CD25.", fillParam="Cohort", 
             title="medianFI CD25 in ICOS+CD38+ cTfh at baseline", yLabel="medianFI CD25 in ICOS+CD38+ cTfh", position="left")  
univScatter(data = subset(subsetData, TimeCategory == "baseline"), yData = "cTfh_ICOShiCD38hi_cTfh..._medfi.CD25.", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Cycle of IT vs +/+ medFI CD25", yLabel= "medianFI CD25 in ICOS+CD38+ cTfh", xLabel = "Cycle of immunotherapy") #+ scale_y_continuous(limits= c(0,12))
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_cTfh..._medfi.CD25. ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CD25hi...FreqParent"))
prePostTimeAveraged(melted, title = "CD25+ ICOS+CD38+ cTfh", xLabel = NULL, yLabel = "CD25+ (% ICOS+CD38+ cTfh)") # + scale_y_continuous(breaks=seq(0,99,0.25))
twoSampleBar(data=subset(subsetData, TimeCategory == "baseline"), xData="Cohort", yData="cTfh_ICOShiCD38hi_CD25hi...FreqParent", fillParam="Cohort", 
             title="CD25+ ICOS+CD38+ cTfh at baseline", yLabel="CD25+ (% ICOS+CD38+ cTfh)", position="left")  
univScatter(data = subset(subsetData, TimeCategory == "baseline"), yData = "cTfh_ICOShiCD38hi_CD25hi...FreqParent", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Cycle of IT vs +/+ CD25+", yLabel= "CD25+ (% ICOS+CD38+ cTfh)", xLabel = "Cycle of immunotherapy") #+ scale_y_continuous(limits= c(0,12))
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_CD25hi...FreqParent ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_Foxp3hi...FreqParent"))
prePostTimeAveraged(melted, title = "Foxp3+ ICOS+CD38+ cTfh", xLabel = NULL, yLabel = "Foxp3+ (% ICOS+CD38+ cTfh)") # + scale_y_continuous(breaks=seq(0,99,0.25))
twoSampleBar(data=subset(subsetData, TimeCategory == "baseline"), xData="Cohort", yData="cTfh_ICOShiCD38hi_Foxp3hi...FreqParent", fillParam="Cohort", 
             title="Foxp3+ ICOS+CD38+ cTfh at baseline", yLabel="Foxp3+ (% ICOS+CD38+ cTfh)", position="left")  
twoSampleBar(data=subset(subsetData, TimeCategory == "oneWeek"), xData="Cohort", yData="cTfh_ICOShiCD38hi_Foxp3hi...FreqParent", fillParam="Cohort", 
             title="Foxp3+ ICOS+CD38+ cTfh at baseline", yLabel="Foxp3+ (% ICOS+CD38+ cTfh)", position="left")  
univScatter(data = subset(subsetData, TimeCategory == "baseline"), yData = "cTfh_ICOShiCD38hi_Foxp3hi...FreqParent", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Cycle of IT vs +/+ Foxp3+", yLabel= "Foxp3+ (% ICOS+CD38+ cTfh)", xLabel = "Cycle of immunotherapy") #+ scale_y_continuous(limits= c(0,12))
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_Foxp3hi...FreqParent ~ Cohort + Year + TimeCategory))


subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CD25.Foxp3....FreqParent"))
prePostTimeAveraged(melted, title = "CD25+Foxp3+ ICOS+CD38+ cTfh", xLabel = NULL, yLabel = "CD25+Foxp3+ (% ICOS+CD38+ cTfh)") # + scale_y_continuous(breaks=seq(0,99,0.25))
twoSampleBar(data=subset(subsetData, TimeCategory == "baseline"), xData="Cohort", yData="cTfh_ICOShiCD38hi_CD25.Foxp3....FreqParent", fillParam="Cohort", 
             title="CD25+Foxp3+ ICOS+CD38+ cTfh at baseline", yLabel="CD25+Foxp3+ (% ICOS+CD38+ cTfh)", position="left")  
univScatter(data = subset(subsetData, TimeCategory == "baseline"), yData = "cTfh_ICOShiCD38hi_CD25.Foxp3....FreqParent", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Cycle of IT vs +/+ CD25+Foxp3+", yLabel= "CD25+Foxp3+ (% ICOS+CD38+ cTfh)", xLabel = "Cycle of immunotherapy") #+ scale_y_continuous(limits= c(0,12))
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_CD25hi...FreqParent ~ Cohort + Year + TimeCategory))


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





#' ## ----------- Cellular subset frequencies at each time point  --------------------
#'
#'

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="cTfh +/+ at d0", yLabel="ICOS+CD38+ cTfh frequency at d0")

# subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 0)
# twoSampleBar(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="ICOS+CD38+ cTfh at d7", yLabel="ICOS+CD38+ cTfh (% CXCR5+)") + 
#   coord_cartesian(ylim = c(0.5,11)) + scale_y_continuous(breaks=seq(1:12))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResponse_justDay7_AllYrs.pdf", device='pdf', width=7, height=7)

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBar(data=subsetData, xData="Cohort", yData="CD27hiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: Plasmablast at d7", yLabel="CD19+CD27++CD38++ frequency at d7")
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
subsetData$CD27hiCD38hi_..FreqParent <-  subsetData$CD27hiCD38hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 100
twoSampleBar(data=subsetData, xData="Cohort", yData="CD27hiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: PB (freq CD19) at d7", yLabel="CD19+CD27++CD38++ frequency at d7")



subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBar(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_CD20loCD71hi...FreqParent", fillParam="Cohort", title="All yrs: Plasmablast at d7", yLabel="IgDloCD71+CD20lo frequency at d7")
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 10000
twoSampleBar(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_CD20loCD71hi...FreqParent", fillParam="Cohort", title="All yrs: Plasmablast at d7", yLabel="IgDloCD71+CD20lo frequency at d7")



subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBar(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_ActBCells...FreqParent", fillParam="Cohort", title="All yrs: ABC at d7", yLabel="ABC frequency at d7")
# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
subsetData$IgDlo_CD71hi_ActBCells...FreqParent <-  subsetData$IgDlo_CD71hi_ActBCells...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent / 10000
twoSampleBar(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_ActBCells...FreqParent", fillParam="Cohort", title="All yrs: ABC at d7", yLabel="ABC frequency at d7")


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
temp <- list()
plotElements <- c("cTfh_ICOShiCD38hi_cTfh..._medfi.SLAM.", "cTfh_ICOShiCD38hi_cTfh..._medfi.CD127.", "cTfh_ICOShiCD38hi_cTfh..._medfi.CXCR5.", "cTfh_ICOShiCD38hi_cTfh..._medfi.Ki67.",
                  "cTfh_ICOShiCD38hi_cTfh..._medfi.PD1.", "cTfh_ICOShiCD38hi_cTfh..._medfi.aIgG4..batchEffect", "cTfh_ICOShiCD38hi_cTfh..._medfi.CD25.", "cTfh_ICOShiCD38hi_cTfh..._medfi.CXCR3.",
                  "cTfh_ICOShiCD38hi_cTfh..._medfi.ICOS.", "cTfh_ICOShiCD38hi_cTfh..._medfi.CD27..batchEffect", "cTfh_ICOShiCD38hi_cTfh..._medfi.CTLA4..batchEffect",   "cTfh_ICOShiCD38hi_cTfh..._medfi.CD38." ,
                  "cTfh_ICOShiCD38hi_cTfh..._medfi.CCR6.",  "cTfh_ICOShiCD38hi_cTfh..._medfi.CD28.", "cTfh_ICOShiCD38hi_cTfh..._medfi.CXCR4.", "cTfh_ICOShiCD38hi_cTfh..._medfi.Foxp3." )
plotNames <- c("SLAM", "CD127","CXCR5","Ki67",
               "PD-1", "aIgG4", "CD25","CXCR3", 
               "ICOS", "CD27", "CTLA4", "CD38",
               "CCR6", "CD28", "CXCR4", "Foxp3")
# for (i in 1:length(plotElements))  { ggsave(twoSampleBar(data=subsetData, xData="Cohort", yData=plotElements[i], fillParam="Cohort", title=plotNames[i], yLabel="medianFI"), filename = paste0("D:/Pembro-Fluvac/Analysis/Images/cTfhPhenotype_",plotNames[i], ".pdf"), width=4, height=5) }

#' ## ----------- Tfh and PB Fold-change at day 7  --------------------
#'
#'

subsetData <- mergedData[,c("cTfh_ICOShiCD38hi_..FreqParent","Subject","TimeCategory")]
temp <- str_split(subsetData$Subject,"-", simplify=T); subsetData <- as.data.frame(cbind(subsetData,temp))
names(subsetData) <- c("cTfhfreq","FullLabel", "TimeCategory","Cohort","Subject")
data_wide <-  pivot_wider(data = subsetData, names_from = c('TimeCategory') , values_from = c("Subject", "cTfhfreq"))
ggplot(data_wide, aes(x=cTfhfreq_baseline, y=cTfhfreq_oneWeek, label=FullLabel)) + geom_point() + geom_label(size=3) + theme_bw() 
cor.test(x = data_wide$cTfhfreq_baseline, y=data_wide$cTfhfreq_oneWeek, use="complete.obs", method='kendall')


subsetData <- subset(mergedData, TimeCategory == "oneWeek" & Cohort != "nonPD1" & cTfh_ICOShiCD38hi_..FreqParent > 0)
twoSampleBar(data=subsetData, xData="Cohort", yData="FCtfh_oW", fillParam="Cohort", title="cTfh", yLabel="Fold-change at one week", FCplot=T, ttest = F) + 
  scale_y_continuous(breaks=seq(0,99,1)) + theme(axis.title.y = element_text(vjust=1.9))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResponse_foldchange_AllYrs.pdf", device="pdf", width=3.5)
subsetData %>% group_by(Cohort) %>% shapiro_test(FCtfh_oW)
wilcox_test(subsetData, FCtfh_oW ~ Cohort)

# aggregate( FC_Tfhresponse, by= list( FC_Tfhresponse$Cohort), FUN=mean, na.rm = T)               # orphan calculation, not sure what this goes to?

subsetData <- subset(mergedData, TimeCategory != "late"  & Cohort != "nonPD1")
FC_PBresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD27hiCD38hi_..FreqParent")); 
FC_PBresponse$FC <- FC_PBresponse$`oneWeek`/FC_PBresponse$`baseline`
twoSampleBar(data=FC_PBresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="Plasmablasts", yLabel="Fold-change at one week (% IgD-)",
             FCplot=T) +   scale_y_continuous(breaks=seq(0,99,1))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/PBresponse_foldchange_AllYrs.pdf", device="pdf", width=4.1, height=7)

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Year == "3" & Cohort != "nonPD1")
subsetData$CD27hiCD38hi_..FreqParent <-  subsetData$CD27hiCD38hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent /100
FC_PBresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD27hiCD38hi_..FreqParent")); 
FC_PBresponse$FC <- FC_PBresponse$`oneWeek`/FC_PBresponse$`baseline`
twoSampleBox(data=FC_PBresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="PB response (freq CD19) at one week", yLabel="CD19+CD27++CD38++ fold-change") + 
  scale_y_continuous(breaks=seq(0,99,1))





subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent > 0)
FC_ASCresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent")); 
FC_ASCresponse$FC <- FC_ASCresponse$`oneWeek`/FC_ASCresponse$`baseline`
twoSampleBar(data=FC_ASCresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="ASC", yLabel="ASC fold-change (% CD71+)", FCplot = T) + 
  scale_y_continuous(breaks=seq(0,99,1), limits = c(0,3)) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/ASCresponse_foldchange_AllYrs.pdf", device="pdf", width=4, height=7)

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent > 0)
subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  subsetData$IgDlo_CD71hi_CD20loCD71hi...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent
FC_ASCresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_CD20loCD71hi...FreqParent")); 
FC_ASCresponse$FC <- FC_ASCresponse$`oneWeek`/FC_ASCresponse$`baseline`
twoSampleBar(data=FC_ASCresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="ASC (freq CD19)", yLabel="ASC fold-change (% CD19)", FCplot=T) + 
  scale_y_continuous(breaks=seq(0,99,1))


subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent > 0)
FC_ABCresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_ActBCells...FreqParent")); 
FC_ABCresponse$FC <- FC_ABCresponse$`oneWeek`/FC_ABCresponse$`baseline`
twoSampleBar(data=FC_ABCresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="ABC", yLabel="ABC fold-change (% CD71+)", FCplot=T) + 
  scale_y_continuous(limits = c(0,4))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/ABCresponse_foldchange_AllYrs.pdf", device="pdf", width=4, height=7)

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent > 0)
subsetData$IgDlo_CD71hi_ActBCells...FreqParent <-  subsetData$IgDlo_CD71hi_ActBCells...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent
FC_ABCresponse <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_ActBCells...FreqParent")); 
FC_ABCresponse$FC <- FC_ABCresponse$`oneWeek`/FC_ABCresponse$`baseline`
twoSampleBar(data=FC_ABCresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="ABC (freq CD19)", yLabel="ABC fold-change (% CD19)")






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
FC_response <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_response$FC_Tfh <- FC_response$`oneWeek`/FC_response$`baseline`
FC_response2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_ActBCells...FreqParent")); 
FC_response$ABC <- FC_response2$`oneWeek`/FC_response2$`baseline`
bivScatter(data1 = FC_response[which(FC_response$Cohort == "Healthy"),], data2 = FC_response[ which(FC_response$Cohort == "aPD1") ,], 
           name1 = "Healthy", name2 = "aPD1", xData = "FC_Tfh", yData = "ABC", fillParam = "Cohort", title = "FC ABC vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "IgDloCD71+CD20hi ABC fold-change") + scale_y_continuous(limits=c(0,4))

# ************** different denominator *********************
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent > 0 )
subsetData$IgDlo_CD71hi_ActBCells...FreqParent <-  subsetData$IgDlo_CD71hi_ActBCells...FreqParent * subsetData$IgDlo_CD71hi_..FreqParent * subsetData$CD19hi_NotNaiveB_freqParent
FC_response <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("cTfh_ICOShiCD38hi_..FreqParent")); 
FC_response$FC_Tfh <- FC_response$`oneWeek`/FC_response$`baseline`
FC_response2 <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgDlo_CD71hi_ActBCells...FreqParent")); 
FC_response$ABC <- FC_response2$`oneWeek`/FC_response2$`baseline`
bivScatter(data1 = FC_response[which(FC_response$Cohort == "Healthy"),], data2 = FC_response[ which(FC_response$Cohort == "aPD1") ,], 
           name1 = "Healthy", name2 = "aPD1", xData = "FC_Tfh", yData = "ABC", fillParam = "Cohort", title = "FC ABC (as freq CD19) vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "IgDloCD71+CD20hi ABC fold-change") + scale_y_continuous(limits=c(0,12))



#' ## ----------- PB vs cTfh response correlations --------------------
#'
#'
## ***************    CD27++CD38++     ********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
univScatter(data = oneWeek, yData = "CD27hiCD38hi_..FreqParent", xData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "+/+ cTfh vs CD27++CD38++ at oneWeek", yLabel= "CD27++CD38++ (freq IgD-)", xLabel = "ICOS+CD38+ cTfh frequency")
subsetData1 <- subset(oneWeek, Cohort == "Healthy" )    
subsetData2 <- subset(oneWeek, Cohort == "aPD1")    
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", yData = "CD27hiCD38hi_..FreqParent", xData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "cTfh vs PB response at oW", yLabel= "CD27++CD38++ (% IgD-CD19+)", xLabel = "ICOS+CD38+ (% cTfh)") + 
  scale_x_continuous(breaks=seq(0,99,1),limits=c(0,12)) + scale_y_continuous(breaks=seq(0,99,2),limits=c(0,15))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_PB-27-38-d7_freqCD71.pdf", device="pdf", width=9, height=7)

oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
univScatter(data = oneWeek, yData = "FCPB_oW", xData = "FCtfh_oW", fillParam = "TimeCategory", 
            title = "+/+ cTfh vs CD27++CD38++ at oneWeek", yLabel= "CD27++CD38++ (freq IgD-)", xLabel = "ICOS+CD38+ cTfh frequency")
subsetData1 <- subset(oneWeek, Cohort == "Healthy" )    
subsetData2 <- subset(oneWeek, Cohort == "aPD1")    
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", yData = "FCPB_oW", xData="FCtfh_oW", 
           fillParam = "Cohort", title = "cTfh vs PB response", yLabel= "Fold-change CD27++CD38++ PB", xLabel = "Fold-change ICOS+CD38+ cTfh") + 
  scale_x_continuous(breaks=seq(0,99,1),limits=c(0,7)) + scale_y_continuous(breaks=seq(0,99,2),limits=c(0,15))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-vs-PB_foldchange.pdf", device="pdf", width=9, height=7)



## ***************    CD71+ CD20loASC    ********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
subsetData1 <- subset(oneWeek, Cohort == "Healthy" )   
subsetData2 <- subset(oneWeek, Cohort == "aPD1"  )    
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", yData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", xData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "Cellular response at one week", yLabel= "ASC (% CD71)", xLabel = "ICOS+CD38+ (% cTfh)") + 
  coord_cartesian(xlim=c(0,12),ylim=c(0,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ASC-71hi-20lo-d7_freqCD71_biv.pdf", device="pdf")
univScatter(data = oneWeek, yData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", xData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "+/+ cTfh vs ASC at one week", yLabel= "ASC (% CD71+)", xLabel = "ICOS+CD38+ (% cTfh)") + coord_cartesian(xlim=c(0,12),ylim=c(0,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ASC-71hi-20lo-d7_freqCD71_univ.pdf", device="pdf")

# ************** different denominator *********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  oneWeek$IgDlo_CD71hi_CD20loCD71hi...FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000
subsetData1 <- subset(oneWeek, Cohort == "Healthy" )   
subsetData2 <- subset(oneWeek, Cohort == "aPD1"  )    
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", yData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", xData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "cTfh vs CD20lo response at oneWeek", yLabel= "ASC (% CD19)", xLabel = "ICOS+CD38+ (% cTfh)") + coord_cartesian(xlim = c(0,12), ylim=c(0,3))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ASC-71hi-20lo-d7_freqCD19_biv.pdf", device="pdf")
univScatter(data = oneWeek, yData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", xData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "+/+ cTfh vs ASC at oneWeek", yLabel= "ASC (% CD19)", xLabel = "ICOS+CD38+ (% cTfh)") + coord_cartesian(xlim = c(0,12), ylim=c(0,3))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ASC-71hi-20lo-d7_freqCD19_univ.pdf", device="pdf")
  

## ***************    CD71+ CD20hi ABC    ********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
subsetData1 <- subset(oneWeek, Cohort == "Healthy"  )    
subsetData2 <- subset(oneWeek, Cohort == "aPD1" ) 
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", yData = "IgDlo_CD71hi_ActBCells...FreqParent", xData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = " cTfh vs ABC at oneWeek", yLabel= "ABC (% CD71+)", xLabel = "ICOS+CD38+ (% cTfh)") + 
  scale_x_continuous(breaks=seq(0,99,1), limits=c(0,10)) + scale_y_continuous(breaks=seq(0,100,25), limits = c(0,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ABC-71hi-20hi-d7_freqCD71_biv.pdf", device="pdf")
univScatter(data = oneWeek, yData = "IgDlo_CD71hi_ActBCells...FreqParent", xData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "+/+ cTfh vs ABC at oneWeek", yLabel= "ABC (% CD71+)", xLabel = "ICOS+CD38+ (% cTfh)")  + 
  scale_x_continuous(breaks=seq(0,99,1), limits=c(0,10)) + scale_y_continuous(breaks=seq(0,100,25), limits = c(0,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ABC-71hi-20hi-d7_freqCD71_univ.pdf", device="pdf")

# ************** different denominator *********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$IgDlo_CD71hi_ActBCells...FreqParent <-  oneWeek$IgDlo_CD71hi_ActBCells...FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000
subsetData1 <- subset(oneWeek, Cohort == "Healthy"   )   
subsetData2 <- subset(oneWeek, Cohort == "aPD1"  )      
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", yData = "IgDlo_CD71hi_ActBCells...FreqParent", xData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "cTfh vs ABC at oneWeek", yLabel= "ABC (% CD19)", xLabel = "ICOS+CD38+ (% cTfh)") + coord_cartesian(xlim=c(0,12),ylim=c(0,3))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ABC-71hi-20hi-d7_freqCD19_biv.pdf", device="pdf")
univScatter(data = oneWeek, yData = "IgDlo_CD71hi_ActBCells...FreqParent", xData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "+/+ cTfh vs CD20hi at oneWeek", yLabel= "ABC (% CD19)", xLabel = "ICOS+CD38+ (% cTfh)") + coord_cartesian(xlim=c(0,12),ylim=c(0,3))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ABC-71hi-20hi-d7_freqCD19_univ.pdf", device="pdf")


oneWeek$CD27hiCD38hi_..FreqParent <-  oneWeek$CD27hiCD38hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent /100
oneWeek$IgDlo_CD71hi_ActBCells...FreqParent <-  oneWeek$IgDlo_CD71hi_ActBCells...FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000

cTfhCorrel <- function( data, subset )
{
  z <- vector()
  z[1] <- subset
  z[2] <- round(cor(data[ which(data$Cohort == "Healthy"), subset], data[ which(data$Cohort == "Healthy") ,"cTfh_ICOShiCD38hi_..FreqParent"], 
                    method = "pearson", use = "complete.obs"), 2)
  z[3] <- round(cor(data[ which(data$Cohort == "aPD1"), subset], data[ which(data$Cohort == "aPD1"),"cTfh_ICOShiCD38hi_..FreqParent"], 
                    method = "pearson", use = "complete.obs"), 2)
  return(z)
}
result <- data.frame(CellType = 0, Healthy = 0, aPD1 = 0)
result[1,] <- cTfhCorrel(oneWeek, subset = "IgDlo_CD71hi_ActBCells...FreqParent")
# fit <- lm(oneWeek$CD27hiCD38hi_..FreqParent ~  oneWeek$cTfh_ICOShiCD38hi_..FreqParent); outlier <- which(cooks.distance(fit) == max(cooks.distance(fit)))       # one row identified with high Cook's distance
# result[2,] <- cTfhCorrel(oneWeek[ -as.numeric(names(outlier)), ], subset = "CD27hiCD38hi_..FreqParent")
result[2,] <- cTfhCorrel(oneWeek, subset = "IgDlo_CD71hi_CD20loCD71hi...FreqParent")
result$CellType <- str_replace(result$CellType, "IgDlo_CD71hi_ActBCells...FreqParent", "ABC \n(% CD19)")
#result$CellType <- str_replace(result$CellType, "CD27hiCD38hi_..FreqParent", "CD27+CD38+ PB")
result$CellType <- str_replace(result$CellType, "IgDlo_CD71hi_CD20loCD71hi...FreqParent", "ASC \n(% CD19)")
result <- melt(result, id.vars = "CellType", measure.vars = c("Healthy","aPD1")); result$value <- as.numeric(result$value)
result$CellType <- factor(result$CellType, levels = c("ASC \n(% CD19)",  "ABC \n(% CD19)"))
ggplot(result, aes(x=CellType, fill=variable, color="bl")) + geom_bar(aes(y = value), stat='identity', position="dodge" ) + theme_bw()  + 
  scale_fill_manual(values=c("#7FAEDB", "#FFB18C")) + scale_color_manual(values="black") + ggtitle("Correlation with \nTfh responses at one week") + 
  ylab("Pearson r") + 
  theme(axis.text.x = element_text(angle = 0), axis.title = element_text(size=16,hjust = 0.5), plot.title = element_text(size=18,hjust = 0.5), 
        axis.title.x = element_blank(), axis.text = element_text(size=16, color="black"), legend.position = "none") + 
  scale_y_continuous(breaks = seq(-1,1,0.1), limits=c(-0.2,0.6))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_various-d7_pearsons_freqCD19.pdf", device = "pdf", width=4, height=6)


#' ## ----------- Vaccine response determinants and biomarkers --------------------
#'
#'

# correlate HAI vs cycle of IT
subsetData <- subset(mergedData,  cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
univScatter(data = subsetData, yData = "cTfh_ICOShiCD38hi_..FreqParent", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Cycle of IT vs +/+ cTfh at bL", yLabel= "ICOS+CD38+ cTfh freq", xLabel = "Cycle of immunotherapy") + coord_cartesian(ylim= c(0,8))
summary(fit <- lm(Cycle.of.Immunotherapy ~ cTfh_ICOShiCD38hi_..FreqParent, subsetData))

univScatter(data = subsetData, yData = "FCtfh_oW", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "cTfh responses", yLabel= "ICOS+CD38+ cTfh FoldChange", xLabel = "Cycle of immunotherapy")  + scale_fill_manual(values=c("#FFB18C"))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CycleIT_vs_cTfhFC.pdf", device='pdf')

univScatter(data = subsetData, yData = "FChai_late", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "HAI fold-change", yLabel= "H1N1pdm09 HAI fold-change", xLabel = "Cycle of immunotherapy")  + coord_cartesian(ylim = c(0,5)) + 
  scale_x_continuous (breaks=seq(0,40,5))+ scale_fill_manual(values=c("#FFB18C"))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CycleIT_vs_HAItiterFC.pdf", device='pdf', width=6, height=6)

summary(fit <- lm(FChai_late ~ FCtfh_oW, data = mergedData))

bivScatter(data1 = subset(mergedData, Cohort == "Healthy"), data2 = subset(mergedData, Cohort == "aPD1"), 
           name1 = "Healthy", name2 = "aPD1", xData = "FCtfh_oW", yData = "FChai_late", fillParam = "Cohort", title = "FC HAI vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "H1N1pdm09 titer fold-change")

# if we take away the sample with the highest Cook's distance mergedData  
outlier <- which(cooks.distance(fit) == max(cooks.distance(fit))) 
subsetData <- mergedData[ - as.numeric(names(outlier)), ]
bivScatter(data1 = subset(subsetData, Cohort == "Healthy"), data2 = subset(subsetData, Cohort == "aPD1") , 
           name1 = "Healthy", name2 = "aPD1", xData = "FCtfh_oW", yData = "FChai_late", fillParam = "Cohort", title = "FC HAI vs FC Tfh", xLabel = "ICOS+CD38+ cTfh fold-change", 
           yLabel = "H1N1pdm09 titer fold-change") #+ scale_y_continuous(trans = "log2", breaks = c(1,2^(1:8))) + coord_cartesian(ylim=c(0.1,128), xlim = c(0,7))
summary(fit <- lm(FChai_late ~ FCtfh_oW   + Cohort,   data = mergedData))


bivScatter(data1 = subset(mergedData, Cohort == "Healthy"), data2 = subset(mergedData, Cohort == "aPD1") , 
           name1 = "Healthy", name2 = "aPD1", xData = "FCPB_oW", yData = "FChai_late", fillParam = "Cohort", title = "FC HAI vs FC PB", xLabel = "PB fold-change", 
           yLabel = "H1N1pdm09 titer fold-change") #+ scale_y_continuous(trans = "log2", breaks = c(1,2^(1:8))) + coord_cartesian(ylim=c(0.1,128), xlim = c(0,7))
summary(fit <- lm(FChai_late ~ FCPB_oW   + Cohort,   data = mergedData))


subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("H1N1pdm09.HAI.titer")); melted <- melted[ which( !is.na(melted$value) ),]
overTime <- aggregate( melted$value, by= list(melted$variable, melted$TimeCategory, melted$Cohort), FUN=mean, na.rm = T); overTime
aPD1_seroprot <- length(which(melted$value > 40 & melted$Cohort == "aPD1" & melted$TimeCategory == "late") )
aPD1_total <- length(which(melted$Cohort == "aPD1" & melted$TimeCategory == "late") )
HC_seroprot <- length(which(melted$value > 40 & melted$Cohort == "Healthy" & melted$TimeCategory == "late") )
HC_total <- length(which(melted$Cohort == "Healthy" & melted$TimeCategory == "late") )
continTable <- matrix(c(aPD1_seroprot, aPD1_total - aPD1_seroprot, HC_seroprot, HC_total - HC_seroprot), ncol=2); continTable
fisher.test(continTable)

seroprot <- data.frame( Cohort = c("Healthy", "aPD1"), Seroprot = c(HC_seroprot / HC_total, aPD1_seroprot / aPD1_total))
seroprot$Cohort <- factor(seroprot$Cohort, levels = c("Healthy", "aPD1"))
ggplot(data=seroprot, aes(x=Cohort, y=Seroprot, fill=Cohort,width=0.6)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) +  
  geom_bar( position = position_dodge(), stat = "identity", color="black",size=0.1) + ggtitle("Seroprotection") + ylab("Proportion seroprotected") +  theme_bw() +
  theme(axis.text = element_text(size=28,hjust = 0.5,color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=32,hjust = 0.5)) + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45,hjust=1,vjust=1)) + 
  coord_cartesian(ylim = c(0,1)) + scale_y_continuous(breaks = seq(0,1,0.1))
# ggsave (filename = "D:/Pembro-Fluvac/Analysis/Images/Seroprotection_byCohort.pdf", device="pdf", width=4)

prePostTimeAveraged(melted, title = "HAI responses", xLabel = NULL, yLabel = "H1N1pdm09 HAI titer") + scale_y_continuous(trans = 'log2', breaks=c(2^(2:14)))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAIresponsesTimeAveraged.pdf", device = "pdf", width=7, height=7)
Anova(fit <- lm( value ~ Cohort * TimeCategory, data = melted)); tukey_hsd(fit)

subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(data = subsetData, xData = "TimeCategory", yData = "H1N1pdm09.HAI.titer", fillParam = "Cohort", title = "HAI titers over time", xLabel = "TimeCategory",
            yLabel = "HAI titer", groupby = "Subject") + scale_y_continuous(trans='log2', breaks=2^(2:13))

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="H1N1pdm09.HAI.titer", fillParam="Cohort", batch = "Year",title="All yrs: HAI titers at d0", yLabel="H1N1pdm09 titer") + 
  scale_y_continuous(trans='log2', breaks=c(2^(2:14) ) ) +  coord_cartesian(ylim = c(4,8192))
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late")
twoSampleBar(data=subsetData, xData="Cohort", yData="H1N1pdm09.HAI.titer", fillParam="Cohort", batch="Year", title="Late HAI titers", yLabel="H1N1pdm09 titer") + 
  scale_y_continuous(trans='log2' , breaks=c(2^(2:15))) + coord_cartesian(ylim = c(4,16000))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAItiter_late_byCohort.pdf", device="pdf", width=4.5)

subsetData <- subset(mergedData, Cohort != "nonPD1" )
seroconv <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("H1N1pdm09.HAI.titer")); 
seroconv$FC <- seroconv$`late`/seroconv$`baseline`
twoSampleBar(data=seroconv, xData="Cohort", yData="FC", fillParam="Cohort", title="Seroconversion", yLabel="Seroconversion factor") + 
  scale_y_continuous(trans='log2',breaks=c(2^(0:14)), limits = c(1,175))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/SeroconversionFactor.pdf", device = "pdf", width=4.5)


oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$CD27hiCD38hi_..FreqParent <-  oneWeek$CD27hiCD38hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent /100
subsetData1 <- subset(oneWeek, Cohort == "Healthy" & Year == "3" )    
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & Year == "3" & cTfh_ICOShiCD38hi_..FreqParent > 1)    #  *****   OUTLIER REMOVED
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD27hiCD38hi_..FreqParent", yData="FChai_late", 
           fillParam = "Cohort", title = "HAI vs PB response at oneWeek", xLabel= "Plasmablast d7 frequency (freq CD19)", yLabel = "Fold-change HAI (late)") #+ scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1))


univScatter(data = subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek"), yData = "FChai_late", xData = "CD19hi_CD23hi...FreqParent", fillParam = "dummy", 
            title = "HAI vs CD19+CD23+", yLabel= "HAI titer fold-change", xLabel = "CD23+ (% CD19+)")  

univScatter(data = subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), yData = "FChai_late", xData = "CD19hi_CD23hi...FreqParent", fillParam = "dummy", 
            title = "HAI vs CD19+CD23+", yLabel= "HAI titer fold-change", xLabel = "CD23+ (% CD19+)")  

subsetData1 <- subset(mergedData, Cohort == "Healthy" & TimeCategory == "oneWeek"); subsetData2 <- subset(mergedData, Cohort == "aPD1" & TimeCategory == "oneWeek")
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD19hi_CD23hi...FreqParent", yData="FChai_late", 
           fillParam = "Cohort", title = "HAI vs CD23 at oneWeek", xLabel= "CD23+ (% CD19)", yLabel = "Fold-change HAI (late)") #+ scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1))
subsetData1 <- subset(mergedData, Cohort == "Healthy" & TimeCategory == "baseline"); subsetData2 <- subset(mergedData, Cohort == "aPD1" & TimeCategory == "baseline")
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD19hi_CD23hi...FreqParent", yData="FChai_late", 
           fillParam = "Cohort", title = "HAI vs CD23 at baseline", xLabel= "CD23+ (% CD19)", yLabel = "Fold-change HAI (late)") #+ scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1))



	


#' ## ----------- CXCL13  --------------------
#'
#'



subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort','Year'), measure.vars = c("Plasma.CXCL13..pg.mL."))
prePostTimeAveraged(melted, title = "CXCL13", xLabel = NULL, yLabel = "Plasma CXCL13 (pg/mL)") + scale_y_continuous(breaks=seq(0,200,5), limits=c(50,110))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_TimeAveraged.pdf", device="pdf", width=4)
summary( fit <- aov(value ~ Cohort+TimeCategory + Cohort:TimeCategory, data=melted ) );  tukey_hsd(fit)

subsetData <- subset(melted, TimeCategory == "oneWeek")
rownames(subsetData) <- seq(1:nrow(subsetData))
twoSampleBar(data=subsetData, xData="Cohort", yData="value", fillParam="Cohort", title="CXCL13\nOne week", yLabel="plasma CXCL13 (pg/mL)") + 
  scale_y_continuous(limits = c(0,250))
aggregate( subsetData, by= list(subsetData$variable, subsetData$TimeCategory, subsetData$Cohort), FUN=mean, na.rm = T)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_oW.pdf", device="pdf",  width = 4)


FC_response <- dcast( subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" ), `Subject`+`Cohort`~`TimeCategory`, value.var = c("Plasma.CXCL13..pg.mL.")) 
FC_response$FCcxcl13 <- FC_response$oneWeek / FC_response$baseline
t.test(FC_response[which(FC_response$Cohort == "aPD1"), "baseline"], FC_response[which(FC_response$Cohort == "aPD1"), "oneWeek"], paired = T)
t.test(FC_response[which(FC_response$Cohort == "Healthy"), "baseline"], FC_response[which(FC_response$Cohort == "Healthy"), "oneWeek"], paired = T)


twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek"), xData="Cohort", yData="FCCXCL13_oW", fillParam="Cohort",FCplot =T, 
             title="CXCL13", yLabel="Fold-change at one week")+ scale_y_continuous(limits=c(0,3), breaks = seq(0,50,0.5))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_foldChange_oW.pdf", width = 4)


subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "Plasma.CXCL13..pg.mL.", fillParam = "Cohort", title = "CXCL13", xLabel = "TimeCategory",
            yLabel = "Plasma CXCL13 (pg/mL)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,20))     # limits = c(0,45)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_prePostTime.pdf", device="pdf")


univScatter(data = subset(mergedData, TimeCategory == "baseline" & Cohort == "aPD1"), yData = "Plasma.CXCL13..pg.mL.",  fillParam = "dummy",
            xData = "IgDlo_CD71hi_CD20loCD71hi.CD23hi...FreqParent", yLabel = "CXCL13 (pg/mL)", xLabel = "CD23hi (% ASC)", title = "CXCL13 vs CD23+ ASC at bL")

univScatter(data = subset(mergedData, TimeCategory == "late" & Cohort == "aPD1"), yData = "Plasma.CXCL13..pg.mL.",  fillParam = "dummy",
            xData = "IgDlo_CD71hi_CD20loCD71hi.CD23hi...FreqParent", yLabel = "CXCL13 (pg/mL)", xLabel = "CD23hi (% ASC)", title = "CXCL13 vs CD23+ ASC at bL")


univScatter(data = subset(mergedData, TimeCategory == "baseline" & Cohort == "aPD1"), xData = "Plasma.CXCL13..pg.mL.",  fillParam = "dummy", position = 'right',
            yData = "FCCXCL13_oW", xLabel = "baseline CXCL13 (pg/mL)", yLabel = "Fold-change CXCL13", title = "CXCL13 in aPD1") + scale_fill_manual(values=c("#FFB18C"))+
  scale_y_continuous(breaks = seq(0,5,1),limits=c(0,3))

# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_bL-vs-FC_aPD1.pdf")


univScatter(data = subset(mergedData, TimeCategory == "late" & Cohort == "aPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent!=0), xData = "Plasma.CXCL13..pg.mL.",  fillParam = "dummy",
            yData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", xLabel = "baseline CXCL13 (pg/mL)", yLabel = "ASC (% CD71+)", title = "ASC vs CXCL13 at bL")



univScatter(data = mergedData, yData = "FCCXCL13_oW", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "AllYrs: Cycle of IT vs CXCL13 FC", yLabel= "CXCL13 FoldChange", xLabel = "Cycle of immunotherapy")  + scale_y_continuous(limits = c(0,4))
univScatter(data = mergedData, yData = "FCCXCL13_oW", xData = "CD27hiCD38hi_..FreqParent", fillParam = "dummy", 
            title = "AllYrs: PB response vs CXCL13 FC", yLabel= "CXCL13 FoldChange", xLabel = "CD27+CD38+ frequency") # + coord_cartesian(ylim = c(0,2))
univScatter(data = mergedData, yData = "FCPB_oW", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "AllYrs: Cycle of IT vs PB FC", yLabel= "PB FoldChange", xLabel = "Cycle of immunotherapy")


subsetData <- subset(mergedData, Cohort != "nonPD1" )
crossTimeCategory <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c('Plasma.CXCL13..pg.mL.' , 'H1N1pdm09.HAI.titer') )
crossTimeCategory2 <- dcast(crossTimeCategory, Subject + Cohort ~ TimeCategory + variable)
univScatter(crossTimeCategory2, xData = "baseline_Plasma.CXCL13..pg.mL.", yData="late_H1N1pdm09.HAI.titer", fillParam = NULL, 
            title = "All yrs: Late HAI vs baseline CXCL13", xLabel= "baseline Plasma CXCL13 (pg/mL)", yLabel = "Late HAI titer - H1N1pdm09")   # nonstd appearance due to lack of fillParam
univScatter(crossTimeCategory2, xData = "oneWeek_Plasma.CXCL13..pg.mL.", yData="late_H1N1pdm09.HAI.titer", fillParam = NULL, 
            title = "All yrs: Late HAI vs d7 CXCL13", xLabel= "oneWeek Plasma CXCL13 (pg/mL)", yLabel = "Late HAI titer - H1N1pdm09")   # nonstd appearance due to lack of fillParam
summary(lm( late_H1N1pdm09.HAI.titer ~ oneWeek_Plasma.CXCL13..pg.mL. + baseline_H1N1pdm09.HAI.titer, data=crossTimeCategory2))


summary(lm( Plasma.CXCL13..pg.mL. ~ Cohort + Year + TimeCategory, data=mergedData))
summary(lm( FCCXCL13_oW ~ Cohort + FCPB_oW + FChai_late, data=mergedData))

bivScatter(data1 = subset(subsetData, Cohort == "Healthy" & TimeCategory == "baseline"), data2 = subset(subsetData, Cohort == "aPD1" & TimeCategory == "baseline"), 
           name1 = "HC", name2 = "aPD1", xData = "cTfh_ICOSloCD38lo_Foxp3hi...FreqParent", yData="FCCXCL13_oW", 
           fillParam = "Cohort", title = "CXCL13 FC vs Tfr at bL", xLabel= "ICOSloCD38lo Foxp3hi at bL", yLabel = "Fold-change CXCL13") + 
  scale_y_continuous(breaks=seq(0,99,1), limits = c(0,4))

bivScatter(data1 = subset(subsetData, Cohort == "Healthy" & TimeCategory == "baseline"), data2 = subset(subsetData, Cohort == "aPD1" & TimeCategory == "baseline"), 
           name1 = "HC", name2 = "aPD1", xData = "FChai_late", yData="FCCXCL13_oW", 
           fillParam = "Cohort", title = "CXCL13 FC vs HAI FC", xLabel= "FChai_late", yLabel = "Fold-change CXCL13") + 
  scale_y_continuous(breaks=seq(0,99,1), limits = c(0,4))




#' ## ----------- glycosylation data  --------------------
#'



# ***************************    Total afucosylation  *******************************
subsetData <- subset(mergedData, Cohort != "nonPD1" )

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.afucosylated", fillParam = "Cohort", title = "IgG1 afucosylation", xLabel = "TimeCategory",
            yLabel = "anti-H1 Afucosylation (% abund)", groupby = "Subject")  + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1afucosylation_overTime.pdf")

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.afucosylated", fillParam = "Cohort", title = "IgG2 afucosylation over time", xLabel = "TimeCategory",
            yLabel = "Afucosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30))
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.afucosylated", fillParam = "Cohort", title = "IgG3/4 afucosylation over time", xLabel = "TimeCategory",
            yLabel = "Afucosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30))

# ***************************    Total Bisection   *******************************
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.Bisection", fillParam = "Cohort", title = "IgG1 bisection", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG1", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1bisection_overTime.pdf")

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.Bisection", fillParam = "Cohort", title = "IgG2 bisection over time", xLabel = "TimeCategory",
            yLabel = "Bisection (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.Bisection", fillParam = "Cohort", title = "IgG3/4 bisection over time", xLabel = "TimeCategory",
            yLabel = "Bisection (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )

twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), xData="Cohort", yData="IgG1_Total.Bisection", fillParam="Cohort", 
             title="Baseline", yLabel="% anti-H1 IgG1")+ scale_y_continuous(limits = c(0,20), breaks = seq(0,50,5))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1bisection_bL.pdf", width = 4)
twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek"), xData="Cohort", yData="IgG1_Total.Bisection", fillParam="Cohort", 
             title="One Week", yLabel="% anti-H1 IgG1")+ scale_y_continuous(limits = c(0,20), breaks = seq(0,50,5))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1bisection_oW.pdf", width=4)


univScatter(data = subset(mergedData, TimeCategory == "baseline"), yData = "IgG1_Total.Bisection", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "IgG1 bisected at baseline", yLabel= "IgG1 bisected", xLabel = "Cycle of immunotherapy" )


# ***************************    Galactosylation *******************************

melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort','Year'), measure.vars = c("IgG1_Total.Galactosylation..G1.G2."))
prePostTimeAveraged(melted, title = "Total galactosylation", xLabel = NULL, yLabel = "% anti-H1 IgG1") + 
  scale_y_continuous(breaks=seq(0,99,1), limits = c(90,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_totGalactosylation_overTime.pdf")
melted %>%  anova_test(value ~ TimeCategory+Cohort )                        # two-way anova without interaction term 

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.Galactosylation..G1.G2.", fillParam = "Cohort", title = "IgG1 total galactosylation", 
            xLabel = "TimeCategory", yLabel = "% anti-H1 IgG1", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) ) + 
  coord_cartesian(ylim = c(90,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_Tot-galactosylation_overTime_perSubject.pdf", width=8)

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.Galactosylation..G1.G2.", fillParam = "Cohort", title = "IgG2 total galactosylation", 
            xLabel = "TimeCategory", yLabel = "% anti-H1 IgG2", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )+ 
  coord_cartesian(ylim = c(80,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG2_Tot-galactosylation_overTime_perSubject.pdf", width=8)

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.Galactosylation..G1.G2.", fillParam = "Cohort", title = "IgG3/4 total galactosylation", 
            xLabel = "TimeCategory", yLabel = "% anti-H1 IgG3/4", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )+ 
  coord_cartesian(ylim = c(50,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG34_Tot-galactosylation_overTime_perSubject.pdf", width=8)


twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), xData="Cohort", yData="IgG1_Total.Galactosylation..G1.G2.", fillParam="Cohort", 
             title="Total galactosylation\nbaseline", yLabel="% anti-H1 IgG1") + 
  scale_y_continuous(breaks = seq(0,100,2)) +   coord_cartesian(ylim = c(88,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_totGalactosylation_bL.pdf", width=4.5)

univScatter(data = subset(mergedData, TimeCategory == "baseline"), yData = "IgG1_Total.Galactosylation..G1.G2.", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Total galact at baseline", yLabel= "IgG1 total galac", xLabel = "Cycle of immunotherapy" )

univScatter(data = subset(mergedData, TimeCategory == "baseline" & Cohort == "aPD1"), xData = "IgG1_Total.Galactosylation..G1.G2.", yData = "Plasma.CXCL13..pg.mL.", fillParam = "dummy", 
            title = "Total galact vs CXCL13 at baseline", xLabel= "IgG1 total galac", yLabel = "CXCL13"  )
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_totGalactosylation_bL_vs_CXCL13.pdf", width=8)



prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.G0", fillParam = "Cohort", title = "IgG1 G0 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG1", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.G0", fillParam = "Cohort", title = "IgG2 G0 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG2", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.G0", fillParam = "Cohort", title = "IgG3/4 G0 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG3/4", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )

twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), xData="Cohort", yData="IgG1_Total.G0", fillParam="Cohort", 
             title="G0 Galactosylation\nbaseline", yLabel="% anti-H1 IgG1") + scale_y_continuous(breaks = seq(0,100,2)) + 
  coord_cartesian(ylim = c(0,10))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_G0_Galactosylation_bL.pdf", width=4.5)


univScatter(data = subset(mergedData, TimeCategory == "baseline"), yData = "IgG1_Total.G0", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "IgG1 G0 galact at baseline", yLabel= "IgG1 G0 galac", xLabel = "Cycle of immunotherapy" )



prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.G1", fillParam = "Cohort", title = "IgG1 G1 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G1 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.G1", fillParam = "Cohort", title = "IgG2 G1 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G1 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.G1", fillParam = "Cohort", title = "IgG3/4 G1 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G1 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )

twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), xData="Cohort", yData="IgG1_Total.G1", fillParam="Cohort", 
             title="G1 Galactosylation\nbaseline", yLabel="% anti-H1 IgG1") + scale_y_continuous(limits = c(0,70), breaks = seq(0,100,10)) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_G1_Galactosylation_bL.pdf", width=4.5)

univScatter(data = subset(mergedData, TimeCategory == "baseline"), yData = "IgG1_Total.G1", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "IgG1 G1 galact at baseline", yLabel= "IgG1 G1 galac", xLabel = "Cycle of immunotherapy" )



prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.G2", fillParam = "Cohort", title = "IgG1 G2 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G2 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,80) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.G2", fillParam = "Cohort", title = "IgG2 G2 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G2 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,80) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.G2", fillParam = "Cohort", title = "IgG3/4 G2 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G2 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,80) )

twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), xData="Cohort", yData="IgG1_Total.G2", fillParam="Cohort", 
             title="G2 Galactosylation\nbaseline", yLabel="% anti-H1 IgG1") + scale_y_continuous(limits = c(0,80), breaks = seq(0,100,10)) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_G2_Galactosylation_bL.pdf", width=4.5)




# ***************************    Sialylation *******************************

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort','Year'), measure.vars = c("IgG1_Total.sialylated"))
prePostTimeAveraged(melted, title = "Sialylation", xLabel = NULL, yLabel = "% anti-H1 IgG1") + scale_y_continuous(breaks=seq(0,99,2), limits = c(10,22))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_sialylation_overTime.pdf")

melted %>%  anova_test(value ~ TimeCategory+Cohort )                        # two-way anova without interaction term 


twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), xData="Cohort", yData="IgG1_Total.sialylated", fillParam="Cohort", 
             title="Baseline", yLabel="anti-H1 IgG1 Sialylation (% abund)") + scale_y_continuous(limits = c(0,40), breaks = seq(0,50,5))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1sialylation_bL.pdf", width=4)
twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek"), xData="Cohort", yData="IgG1_Total.sialylated", fillParam="Cohort", 
             title="One week", yLabel="anti-H1 IgG1 Sialylation (% abund)")+ scale_y_continuous(limits = c(0,40), breaks = seq(0,50,5))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1sialylation_oW.pdf", width=4)

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.sialylated", fillParam = "Cohort", title = "IgG1 sialylation", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG1", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,10), limits = c(0,40) )
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_sialylation_overTime_perSubject.pdf")
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.sialylated", fillParam = "Cohort", title = "IgG2 sialylation", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG2", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,10), limits = c(0,40) )
#  ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG2_sialylation_overTime_perSubject.pdf")
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.sialylated", fillParam = "Cohort", title = "IgG3/4 sialylation", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG3/4", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,10), limits = c(0,40) )
#  ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG34_sialylation_overTime_perSubject.pdf")


summary(lm( FChai_late ~ Cohort + Year + IgG1sial_oW, data=subsetData))
summary(lm( IgG1sial_oW ~ Cohort + FCtfh_oW, data=subsetData))

subsetData1 <- subset(oneWeek, Cohort == "Healthy" & Year == "3" )    
subsetData2 <- subset(oneWeek, Cohort == "aPD1" & Year == "3" & cTfh_ICOShiCD38hi_..FreqParent > 1)    #  *****   OUTLIER REMOVED
bivScatter(data1 = subset(subsetData, Cohort == "Healthy"), data2 = subset(subsetData, Cohort == "aPD1"), 
           name1 = "HC", name2 = "aPD1", xData = "IgG1sial_oW", yData="FChai_late", 
           fillParam = "Cohort", title = "HAI vs IgG1 sial FC", xLabel= "IgG1 sialylation fold-change", yLabel = "Fold-change HAI (late)") #+ scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1))
bivScatter(data1 = subset(subsetData, Cohort == "Healthy"), data2 = subset(subsetData, Cohort == "aPD1"), 
           name1 = "HC", name2 = "aPD1", xData = "IgG1sial_oW", yData="FCPB_oW", 
           fillParam = "Cohort", title = "PB FC vs IgG1 sial FC", xLabel= "IgG1 sialylation fold-change", yLabel = "Fold-change Plasmablasts") #+ scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1))
bivScatter(data1 = subset(subsetData, Cohort == "Healthy"), data2 = subset(subsetData, Cohort == "aPD1"), 
           name1 = "HC", name2 = "aPD1", xData = "IgG1_Total.sialylated", yData="FCtfh_oW", 
           fillParam = "Cohort", title = "Tfh FC vs IgG1 sial", xLabel= "IgG1 sialylation total", yLabel = "Fold-change Tfh") #+ scale_x_continuous(breaks=seq(0,99,5)) + scale_y_continuous(breaks=seq(0,99,1))
univScatter(data = subset(mergedData, TimeCategory == "oneWeek" & Cohort == "aPD1"), yData = "IgG1_Total.sialylated", xData = "FCtfh_oW", fillParam = "dummy", 
            title = "FCtfh_oW vs IgG1 sialylated", yLabel= "IgG1 sialylated", xLabel = "FCtfh oW" )


prePostTimeGene(singleGeneData, xData = "TimeCategory" , yData = "value", fillParam = "Cohort", groupby = "Subject", title = "FUT8 in ABC", xLabel = "TimeCategory", 
                yLabel = "log2 counts", paired = F)

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory != "late" )

t_test(data = subset(subsetData, Cohort == "Healthy"), IgG1_Total.sialylated~ TimeCategory, paired=T)
t_test(data = subset(subsetData, Cohort == "aPD1"), IgG1_Total.sialylated ~ TimeCategory, paired=T)

univScatter(data = subset(mergedData, TimeCategory == "baseline"), yData = "IgG1_Total.sialylated", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Cycle of IT vs IgG1 sialylated", yLabel= "1gG1 sialylated", xLabel = "Cycle of immunotherapy" )

univScatter(data = subset(mergedData, TimeCategory == "baseline" & Cohort == "aPD1"), yData = "IgG1_Total.sialylated", xData = "IgG1_Total.Galactosylation..G1.G2.", 
            fillParam = "dummy", title = "aPD1: Sialylation vs Glycosylation", yLabel= "Sialylation (% anti-H1 IgG1)", xLabel = "Total galactosylation (% anti-H1 IgG1)" )
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_sial-vs-TotGalactosyl_correlation_aPD1.pdf", width=8)


univScatter(data = subset(mergedData, TimeCategory == "baseline" & Cohort == "Healthy"), yData = "IgG1_Total.sialylated", xData = "IgG1_Total.Galactosylation..G1.G2.", 
            fillParam = "dummy", title = "Healthy: Sialylation vs Glycosylation", yLabel= "Sialylation (% anti-H1 IgG1)", xLabel = "Total galactosylation (% anti-H1 IgG1)" )
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_sial-vs-TotGalactosyl_correlation_HC.pdf", width=8)


# ***************************    Affinity   *******************************

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort','Year'), measure.vars = c("Lo.Hi.HA.affinity"))
prePostTimeAveraged(melted, title = "anti-HA Affinity", xLabel = NULL, yLabel = "Lo Hi HA affinity") + scale_y_continuous(breaks=seq(0,99,0.03))
# ggsave (filename = "D:/Pembro-Fluvac/Analysis/Images/antiHA-affinity_overTime_LoHi.pdf")
melted %>%  anova_test(value ~ TimeCategory+Cohort )                        # two-way anova without interaction term 

univScatter(data = subsetData, xData = "Lo.Hi.HA.affinity", yData = "H1N1pdm09.HAI.titer", fillParam = "dummy", 
            title = "HAI titer incr with affinity", xLabel= "Lo Hi HA affinity", yLabel = "HAI titer - H1N1pdm09") + scale_y_continuous(trans = "log", breaks = 2^(1:12)) +
  coord_cartesian(ylim = c(2,4096))
# ggsave (filename = "D:/Pembro-Fluvac/Analysis/Images/NeutAbtiter_vs_antiHA-affinity_LoHi.pdf")

prePostTime(subsetData, xData = "TimeCategory", yData = "Lo.Hi.HA.affinity", fillParam = "Cohort", title = "anti-HA affinity over time", xLabel = "TimeCategory",
            yLabel = "Lo Hi HA affinity (anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0,1) )
# ggsave (filename = "D:/Pembro-Fluvac/Analysis/Images/antiHA-affinity_overTime_LoHi_perSubject.pdf", width=8)

twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), xData="Cohort", yData="Lo.Hi.HA.affinity", fillParam="Cohort", 
             title="anti-HA affinity at baseline", yLabel="Lo Hi HA affinity (anti-H1)") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0,1.1) )
twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek"), xData="Cohort", yData="Lo.Hi.HA.affinity", fillParam="Cohort", 
             title="anti-HA affinity at oneWeek", yLabel="Lo Hi HA affinity (anti-H1)") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0,1.1) )
twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late"), xData="Cohort", yData="Lo.Hi.HA.affinity", fillParam="Cohort", 
             title="anti-HA affinity at late", yLabel="Lo Hi HA affinity (anti-H1)") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0,1.1) )

FC_Affinity <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("Lo.Hi.HA.affinity"))
FC_Affinity$FC <- FC_Affinity$`late`/FC_Affinity$`baseline`
data=FC_Affinity; xData="Cohort"; yData="FC"; fillParam="Cohort"; title="anti-HA Affinity"; yLabel="fold-change (late / baseline)";
# overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T)
# names(overTime)[which(names(overTime) == 'x')] <- yData
justforttest <- data[, c(xData,yData)]
fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p
annotationInfo <- paste0("P = ", round(pValue, 2)); my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28)))
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) + 
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5) + 
  geom_violin(draw_quantiles = 0.5) + 
  #geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity") + 
  geom_point(size=7, pch=21, fill="black", color="white", alpha=0.35, position = position_jitter(width=0.15, seed = 46)) + 
  ggtitle(title) + ylab(yLabel) +  theme_bw() +
  theme(axis.text = element_text(size=28,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=32,hjust = 0.5)) + 
  annotation_custom(my_grob) + theme(legend.position = "none") + scale_y_continuous(breaks=seq(0,99,0.25), limits = c(0.75,1.75))
# ggsave (filename = "D:/Pembro-Fluvac/Analysis/Images/antiHA-affinity_FC_LoHi.pdf", width=4)


subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("X1.Kd.HA")); melted$value <- 1/melted$value
prePostTimeAveraged(melted, title = "HA EC50/Kd by ELISA", xLabel = NULL, yLabel = "1 / Concentration (ug/mL)") + scale_y_continuous(breaks=seq(0,99,0.25))
summary(fit <- lm( data = mergedData, X1.Kd.HA ~ Cohort  + TimeCategory))
summary(fit <- aov(formula = X1.Kd.HA ~ Cohort  + TimeCategory + Cohort:TimeCategory, data = mergedData)); tukey_hsd(fit)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/AntibodyAffinity_overTime.pdf", device="pdf", height=8, width=8)
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"); subsetData$X1.Kd.HA <- 1 / subsetData$X1.Kd.HA
twoSampleBar(data=subsetData, xData="Cohort", yData="X1.Kd.HA", fillParam="Cohort", title="HA EC50/Kd at baseline", yLabel="1 / Concentration (ug/mL)",position = 'left')
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/AntibodyAffinity_baseLine.pdf", device="pdf", width=4)

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
univScatter(data = subsetData, xData = "Cycle.of.Immunotherapy", yData = "Lo.Hi.HA.affinity", fillParam = "TimeCategory", 
            title = "Affinity over time", xLabel= "Cycle of Immunotherapy", yLabel = "Lo Hi Affinity")  + scale_fill_manual(values=c("#FFB18C"))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/AntibodyAffinity_cycleofImmunotherapy.pdf", width=8)



#' ## ----------- Working with the medianFI for key transcription factors  --------------------
#'
#'

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="Ki67hiCD38hicTfh_medfi.Blimp1.", fillParam="Cohort", title="Blimp1 in +/+ at baseline", yLabel="Blimp1 median FI")
twoSampleBar(data=subsetData, xData="Cohort", yData="Ki67hiCD38hicTfh_medfi.Bcl6..batchEffect", fillParam="Cohort", title="Bcl6 in +/+ at baseline", yLabel="Bcl6 median FI")
subsetData$temp <- subsetData$Ki67hiCD38hicTfh_medfi.Bcl6..batchEffect / subsetData$Ki67hiCD38hicTfh_medfi.Blimp1.
twoSampleBar(data=subsetData, xData="Cohort", yData="temp", fillParam="Cohort", title="Bcl6/Blimp1 ratio in +/+ at baseline", yLabel="Bcl6/Blimp1 protein ratio")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBar(data=subsetData, xData="Cohort", yData="Ki67hiCD38hicTfh_medfi.Blimp1.", fillParam="Cohort", title="All yrs: Blimp1 in +/+ at d7", yLabel="Blimp1 median FI")
twoSampleBar(data=subsetData, xData="Cohort", yData="Ki67hiCD38hicTfh_medfi.Bcl6..batchEffect", fillParam="Cohort", title="All yrs: Bcl6 in +/+ at d7", yLabel="Bcl6 median FI")
subsetData$temp <- subsetData$Ki67hiCD38hicTfh_medfi.Bcl6..batchEffect / subsetData$Ki67hiCD38hicTfh_medfi.Blimp1.
twoSampleBar(data=subsetData, xData="Cohort", yData="temp", fillParam="Cohort", title="Bcl6/Blimp1 ratio in +/+ at d7", yLabel="Bcl6/Blimp1 protein ratio")

#' ## ----------- subset protein correlations  --------------------
#'
#'

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
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_PB-71hi-20lo-d7_CD32_univ.pdf")

oneWeek <- subset(mergedData, TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
subsetData1 <- subset(oneWeek, Cohort == "Healthy")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "ActBCells...medfi.CD32.", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "AllYrs: cTfh vs ABC CD32 at oneWeek", xLabel= "ABC - CD32 medFI", yLabel = "ICOS+CD38+ cTfh frequency")
univScatter(data = oneWeek, xData = "ActBCells...medfi.CD32.", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "AllYrs: cTfh vs ABC CD32 at oneWeek", xLabel= "ABC - CD32 medFI", yLabel = "ICOS+CD38+ cTfh frequency")


oneWeek <- subset(mergedData, TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
subsetData1 <- subset(oneWeek, Cohort == "Healthy")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDlo_CD71hi_CD20loCD71hi.CD86....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", 
                fillParam = "Cohort", title = "AllYrs: cTfh vs PB CD86 resp at oneWeek", xLabel= "CD20lo - CD86hi frequency", yLabel = "ICOS+CD38+ cTfh frequency")
univScatter(data = oneWeek, xData = "IgDlo_CD71hi_CD20loCD71hi.CD86....FreqParent", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "AllYrs: cTfh vs CD20lo CD86+ at oneWeek", xLabel= "CD20lo - CD86hi frequency", yLabel = "ICOS+CD38+ cTfh frequency")


oneWeek <- subset(mergedData, TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
subsetData1 <- subset(oneWeek, Cohort == "Healthy")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD19_NonnaiveB.CD27.CD38....medfi.CD32.", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "AllYrs: cTfh vs PB CD32 at oneWeek", xLabel= "CD27++CD38++ PB - CD32 medFI", yLabel = "ICOS+CD38+ cTfh frequency")
univScatter(data = oneWeek, xData = "CD19_NonnaiveB.CD27.CD38....medfi.CD32.", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "TimeCategory", 
            title = "AllYrs: cTfh vs CD27++CD38++ CD32 at oneWeek", xLabel= "CD27++CD38++ PB - CD32 medFI", yLabel = "ICOS+CD38+ cTfh frequency")
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs-PB-27hi38hi-d7_CD32medFI_univ.pdf")


# as predicted by the elastic net regression
oneWeek <- subset(mergedData, TimeCategory == "oneWeek" & cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
subsetData1 <- subset(oneWeek, Cohort == "Healthy")    ; subsetData2 <- subset(oneWeek, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "IgDloCD71hi..medfi.CD86.", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "AllYrs: cTfh vs CD86 resp at oneWeek", xLabel= "IgDloCD71hi - CD86 medFI", yLabel = "ICOS+CD38+ cTfh frequency")
univScatter(data = oneWeek, xData = "IgDloCD71hi..medfi.CD86.", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "TimeCategory", title = "AllYrs: cTfh vs CD86 resp at oneWeek", xLabel= "IgDloCD71hi - CD86 medFI", yLabel = "ICOS+CD38+ cTfh frequency")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek"); 
twoSampleBar(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_cTfh..._medfi.Ki67.", fillParam="Cohort", title="One Week", yLabel="medFI Ki67")
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek"); 
twoSampleBar(data=subsetData, xData="Cohort", yData="CD20loCD71hi...medfi.Ki67.", fillParam="Cohort", title="All yrs: at oneWeek", yLabel="medFI Ki67")
univScatter(data=subsetData, xData = "CD20loCD71hi...medfi.Ki67.", yData="cTfh_ICOShiCD38hi_..FreqParent", 
            fillParam = "TimeCategory", title = "AllYrs: cTfh vs CD86 resp at oneWeek", xLabel= "IgDloCD71hi - CD86 medFI", yLabel = "ICOS+CD38+ cTfh frequency")
subsetData1 <- subset(subsetData, Cohort == "Healthy")    ; subsetData2 <- subset(subsetData, Cohort == "aPD1")   
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "aPD1", xData = "CD20loCD71hi...medfi.Ki67.", yData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "AllYrs: cTfh vs CD86 resp at oneWeek", xLabel= "IgDloCD71hi - CD86 medFI", yLabel = "ICOS+CD38+ cTfh frequency")




#' ## ----------- irAE analysis  --------------------
#'
#'

index <- demog[which(demog$irAE == "Y" | demog$irAE == "N"),"Subject"]
mergedData.irAE <- mergedData[which(mergedData$Subject %in% index),]  
subsetData <- subset(mergedData.irAE, TimeCategory == "baseline" & Cohort == "aPD1")
subsetData$dateDiff <- round( as.numeric(difftime(subsetData$DateFirstIRAE,  subsetData$DateFirstFlowFile, units = "days")),digits = 0)
subsetData <- subsetData[-which(is.na(subsetData$dateDiff)), ]           # zero out the ones who did not have irAE
subsetData$Subject <- factor(subsetData$Subject, levels = subsetData$Subject[order(subsetData$dateDiff, decreasing = T)] )
ggplot(subsetData, aes(x=dateDiff, y=Subject)) + 
  geom_bar(stat="Identity", width=0.15, fill="#ff9a6a", size=0.01, color="black") + geom_point(aes(size=irAEgradeWorst)) + 
  xlab("Days until flu vaccine") + ylab("Subject") + labs(size = "irAE grade") + ggtitle("Time since first irAE") + scale_size(range=c(2,7),) + 
  theme_bw() + theme(axis.text.x = element_text(color = "black", size=24), axis.title = element_text(size=24), plot.title = element_text(size=32),
                     axis.text.y = element_blank(), legend.text = element_text(size=16), legend.title = element_text(size=16), 
                     ) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_dateDiff_irAEtoVaccine.pdf")


subsetData <- mergedData[which(!is.na(mergedData$irAE) & mergedData$irAE != "" ), ]
twoSampleBar(data = subsetData, yData="FCtfh_oW", yLabel = "Fold-change at one week", title="Post-vax\nICOS+CD38+ cTfh", xData = "irAE",fillParam = "irAE") + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_FCtfh.pdf", width=3.5)
twoSampleBar(data = subsetData, yData="H1N1pdm09.HAI.titer", yLabel = "Baseline H1N1pdm09 titer", title="Baseline HAI", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28)) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_baselineHAI.pdf", width=4)
twoSampleBar(data = subsetData, yData="FChai_late", yLabel = "FC HAI titer", title="Fold-change HAI", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("#FFE6da", "#ff9a6a"))
twoSampleBar(data = subsetData, yData="IgDlo_CD71hi_ActBCells...FreqParent", yLabel = "Baseline ABC (% CD71)", title="ABC", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("#FFE6da", "#ff9a6a"))
twoSampleBar(data = subsetData, yData="IgDlo_CD71hi_CD20loCD71hi...FreqParent", yLabel = "Baseline ASC (% CD71)", title="ASC", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("#FFE6da", "#ff9a6a"))
twoSampleBar(data = subsetData, yData="cTfh_ICOShiCD38hi_cTfh..._medfi.aIgG4..batchEffect", yLabel = "medianFI aIgG4", title="baseline\nICOS+CD38+ cTfh", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a"))  + theme(axis.title.x = element_text(size=28)) # + geom_label_repel(aes(label = Subject),size=3)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_aIgG4.pdf", width=4)
twoSampleBar(data = subsetData, yData="cTfh_ICOShiCD38hi_cTfh..._medfi.ICOS.", yLabel = "medianFI ICOS", title="baseline\nICOS+CD38+ cTfh", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_ICOS.pdf", width=4)
twoSampleBar(data = subsetData, yData="Ki67hiCD38hicTfh_medfi.CD71.", yLabel = "medianFI CD71 - baseline", title="ICOS+CD38+ cTfh", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_CD71.pdf", width=4)







