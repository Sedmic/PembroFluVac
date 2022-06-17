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

##' ## ----------- Save mergedData as RDS  --------------------

mergedData <- readRDS("D:/Pembro-Fluvac/Analysis/mergedData/mergedData.RDS")


##' ## ----------- Cohort 2 description --------------------
#'
#'
 

table(demog$Current.Immunotherapy[which(demog$Current.Immunotherapy != "Ipi/Nivo")])
x <- demog[which(demog$Cohort == "anti-PD-1"),] 
table(x$Sex)
x <- demog[which(demog$Cohort == "Healthy"),] 
table(x$Sex)


table(mergedData$Cohort, mergedData$Sex, mergedData$Year)  
table(serology$Cohort, serology$TimeCategory, serology$Year, !is.na(serology$H1N1pdm09.HAI.titer))
subsetData <- subset(demog, Current.Immunotherapy %in% c("Pembrolizumab","Nivolumab"))
table(subsetData$Current.Immunotherapy)
table(subsetData$Cohort, subsetData$Sex)
table(demog$Cohort, demog$Sex)
subsetData <- subset(mergedData, Cohort == "anti-PD-1" & TimeCategory == "baseline")
ggplot( data = subsetData, aes(x=Cycle.of.Immunotherapy)) + geom_histogram(color="white",fill="#FFB18C") + 
  ggtitle("Cohort 2 - anti-PD1 cohort") +  theme_bw() + 
  theme(axis.text=element_text(color="black",size=18), axis.title=element_text(color="black",size=24), plot.title = element_text(size=24)) + scale_x_continuous(breaks=seq(0,30,5)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CyclesImmunotherapy.pdf", device = 'pdf', height=3, width=8)
# write.csv(x = subsetData[,c("dbCode", "Cohort", "TimeCategory", "Cycle.of.Immunotherapy")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1b-2.csv")

median(mergedData[which(mergedData$Cohort == "anti-PD-1" & mergedData$TimeCategory == "baseline"), 'Cycle.of.Immunotherapy'], na.rm = T)  
sd(mergedData[which(mergedData$Cohort == "anti-PD-1" & mergedData$TimeCategory == "baseline"), 'Cycle.of.Immunotherapy'], na.rm = T)  

paste("Median Age anti-PD-1: ", median(mergedData[which(mergedData$Cohort == "anti-PD-1" & mergedData$TimeCategory == "baseline"), 'Age'], na.rm = T)  )
paste("Median Age Healthy: ", median(mergedData[which(mergedData$Cohort == "Healthy" & mergedData$TimeCategory == "baseline"), 'Age'], na.rm = T)   )
paste("Range Age anti-PD-1: ", range(mergedData[which(mergedData$Cohort == "anti-PD-1" & mergedData$TimeCategory == "baseline"), 'Age'], na.rm = T)  )
paste("Range Age Healthy: ", range(mergedData[which(mergedData$Cohort == "Healthy" & mergedData$TimeCategory == "baseline"), 'Age'], na.rm = T)   )



##'  ## ----------- Tfh dynamics   --------------------
#'
#'

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", title="cTfh +/+ at baseline", yLabel="ICOS+CD38+ cTfh frequency at d0")

## excluding PD-1 in the definition of cTfh
subsetData <- subset(mergedData, TimeCategory == "oneWeek" & Cohort != "nonPD1" & cTfh_ICOShiCD38hi_..FreqParent > 0)
twoSampleBar(data=subsetData, xData="Cohort", yData="FCtfh_oW", fillParam="Cohort", title="ICOS+CD38+\ncTfh", yLabel="Fold-change at one week", FCplot=T, ttest = F) + 
  scale_y_continuous(breaks=seq(0,99,1)) + theme(axis.title.y = element_text(vjust=1.9)) + 
  ggtitle(expression(paste(paste("ICO",S^'+', "CD3",8^'+' ), " cTfh")))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResponse_foldchange_AllYrs.pdf", device="pdf", width=3.8, height=7)
subsetData %>% group_by(Cohort) %>% shapiro_test(FCtfh_oW)
wilcox_test(subsetData, FCtfh_oW ~ Cohort)
# write.csv(x = subsetData[,c("dbCode", "Cohort", "TimeCategory", "FCtfh_oW")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/1d.csv")



## including PD-1 in the definition of cTfh
subsetData <- subset(mergedData, TimeCategory == "oneWeek" & Cohort != "nonPD1" & cTfh_ICOShiCD38hi_..FreqParent > 0)
twoSampleBar(data=subsetData, xData="Cohort", yData="X5PD1_ICOSCD38_oW", fillParam="Cohort", title="ICOS+CD38+\nCXCR5+PD1+", yLabel="Fold-change at one week", FCplot=T, ttest = F) + 
  scale_y_continuous(breaks=seq(0,99,1)) + theme(axis.title.y = element_text(vjust=1.9)) + 
  ggtitle(expression(paste(paste("PD-",1^'+',"ICO",S^'+', "CD3",8^'+' ), " cTfh")))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResponse_X5PD1ICOSCD38_foldchange_AllYrs.pdf", device="pdf", width=3.5)
subsetData %>% group_by(Cohort) %>% shapiro_test(X5PD1_ICOSCD38_oW)
wilcox_test(subsetData, X5PD1_ICOSCD38_oW ~ Cohort)
# write.csv(x = subsetData[,c("dbCode", "Cohort", "TimeCategory", "X5PD1_ICOSCD38_oW")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1k.csv")

subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(data = subsetData, xData = "TimeCategory", yData = "cTfh_ICOShiCD38hi_..FreqParent", fillParam = "Cohort", title = "cTfh response", xLabel = "TimeCategory",
            yLabel = "ICOS+CD38+ (% cTfh)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,150,2), limits=c(0,18))     + 
  ylab(expression(paste("ICO", S^'+', "CD3", 8^'+', " (% cTfh)")))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_PerSubject.pdf")
# write.csv(x = subsetData[,c("dbCode", "Cohort", "TimeCategory", "cTfh_ICOShiCD38hi_..FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1f.csv")

melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_..FreqParent"))
melted <- melted[which(melted$TimeCategory != 'late'),]
prePostTimeAveraged(data = melted, title = "cTfh", xLabel = NULL, yLabel = "ICOS+CD38+ (% cTfh)") + scale_y_continuous(breaks=seq(1,10,0.5), limits=c(3.0, 7))+ 
  ylab(expression(paste("ICO", S^'+', "CD3", 8^'+', " (% cTfh)")))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_freqcTfh.pdf", width=4)     # unpaired analysis
summary( fit <- lm(formula = value ~ Cohort * TimeCategory, data = melted) );  tukey_hsd(fit)       # testing using unpaired two-way ANOVA with Tukey's posthoc
melted %>% group_by(Cohort) %>% kruskal_test(formula = value ~ TimeCategory)                        # testing per cohort using Kruskall-Wallis without adjustment
index <- melted[melted$TimeCategory=="oneWeek", 'dbCode']
melted.paired <- melted[which(melted$dbCode %in% index), ]
prePostTimeAveraged(data = melted.paired, title = "cTfh", xLabel = NULL, yLabel = "ICOS+CD38+ (% cTfh)") + scale_y_continuous(breaks=seq(1,10,0.5), limits=c(3.0, 7))+ 
  ylab(expression(paste("ICO", S^'+', "CD3", 8^'+', " (% cTfh)")))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_freqcTfh_paired.pdf", width=4)
# write.csv(melted.paired, file = "cTfh_responses_freqParent.csv")                                  # paired two-way ANOVA results calculated in Prism
# write.csv(x = melted.paired, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1g.csv")


#************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$cTfh_ICOShiCD38hi_..FreqParent <- subsetData$cTfh_ICOShiCD38hi_..FreqParent * subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /10000
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_..FreqParent"))
melted <- melted[which(melted$TimeCategory != 'late'),]
prePostTimeAveraged(melted, title = "cTfh", xLabel = NULL, yLabel = "ICOS+CD38+ cTfh (% CD4)") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0.25,1))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_freqCD4.pdf", width=4)
summary( fit <- lm(formula = value ~ Cohort * TimeCategory, data = melted) ); tukey_hsd(fit)        # testing using unpaired two-way ANOVA with Tukey's posthoc
melted %>% group_by(Cohort) %>% kruskal_test(formula = value ~ TimeCategory)                        # testing per cohort using Kruskal-Wallis without adjustment
                                                                                                    
index <- melted[melted$TimeCategory=="oneWeek", 'dbCode']
melted.paired <- melted[which(melted$dbCode %in% index), ]
prePostTimeAveraged(data = melted.paired, title = "cTfh", xLabel = NULL, yLabel = "ICOS+CD38+ cTfh (% CD4+)") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0.25,1))+ 
  ylab(expression(paste("ICO", S^'+', "CD3", 8^'+', " (% CD", 4^'+', ")")))   # paired two-way ANOVA results calculated in Prism
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_freqCD4_paired.pdf", width=4)
# write.csv(melted.paired, file = "cTfh_responses_freqCD4.csv")
# write.csv(x = melted.paired, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1h.csv")


# ************** different denominator *********************
subsetData <- subset(mergedData, Cohort != "nonPD1")
subsetData$cTfh_ICOShiCD38hi_Ki67hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_Ki67hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent * 
  subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /1000000
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_Ki67hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ Ki67hi frequencies averaged", xLabel = NULL, yLabel = "+/+ Ki67hi (freq CD4)") 
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_Ki67hi...FreqParent ~ Cohort  + TimeCategory + Year))    # Year is a batch effect
tukey_hsd(aov(data = subsetData, formula = cTfh_ICOShiCD38hi_Ki67hi...FreqParent ~ Cohort:TimeCategory))  # proportion of CD4 not FreqParent

subsetData <- subset(subsetData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBar(data=subsetData, xData="Cohort", yData="cTfh_ICOShiCD38hi_Ki67hi...FreqParent", fillParam="Cohort", title="Ki67 in ICOS+CD38+ cTfh at oW", 
             yLabel="+/+ Ki67hi (% CD4) ") 
data1 = subset(subsetData, Cohort == "Healthy");  data2 = subset(subsetData, Cohort == "anti-PD-1")
bivScatter(data1 = data1, data2 = data2, name1 = "Healthy", name2 = "anti-PD-1",
           xData = "FCtfh_oW", yData = "cTfh_ICOShiCD38hi_Ki67hi...FreqParent", fillParam = "Cohort",
           title = expression(paste("Ki6",7^"+", " cTfh at one week")),
           xLabel = expression(paste("Fold-change in ICO", S^'+', "CD3", 8^'+'," cTfh")), 
           yLabel = expression(paste("Ki6",7^'+',"ICO", S^'+', "CD3", 8^'+', " cTfh (% CD4)")),
           nonparam = T) + scale_y_continuous(limits=c(0,3.5))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Ki67-vs-cTfhFoldChange_biv.pdf")
temp <- rbind(data1, data2)
# write.csv(temp[,c("dbCode","TimeCategory","FCtfh_oW","cTfh_ICOShiCD38hi_Ki67hi...FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e5g.csv")


# ************** relationship to cycle of immunotherapy *********************
subsetData <- subset(mergedData,  cTfh_ICOShiCD38hi_..FreqParent > 1 & Cohort != "nonPD1")    #  *****   OUTLIER REMOVED
univScatter(data = subsetData, yData = "cTfh_ICOShiCD38hi_..FreqParent", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "Cycle of IT vs +/+ cTfh at bL", yLabel= "ICOS+CD38+ cTfh freq", xLabel = "Cycle of immunotherapy") + coord_cartesian(ylim= c(0,8)) 
summary(fit <- lm(Cycle.of.Immunotherapy ~ cTfh_ICOShiCD38hi_..FreqParent, subsetData))

univScatter(data = subsetData, yData = "FCtfh_oW", xData = "Cycle.of.Immunotherapy", fillParam = "dummy", 
            title = "cTfh responses", yLabel= "Fold-change at one week", xLabel = "Cycle of immunotherapy", nonparam=T)  + scale_fill_manual(values=c("#FFB18C"))+ 
  ggtitle(expression(paste("ICO",S^'+',"CD3",8^'+'," cTfh")))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CycleIT_vs_cTfhFC.pdf", device='pdf')
# write.csv(x = subsetData[,c("dbCode","Cohort","Cycle.of.Immunotherapy","FCtfh_oW")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1l.csv")



##' ## ----------- BAA Year 3 analysis - healthy young and elderly adults --------------------
#'

BAAyear3long <- read.csv(file="D:/Pembro-Fluvac/Analysis/mergedData/BAAyear3long.csv"); names(BAAyear3long)[6] <- 'TimePoint'
BAAyear3long$Cohort <- factor(BAAyear3long$Cohort, levels=c("Young","Elderly"))
temp <- BAAyear3long[!duplicated(BAAyear3long$Subject),c("Subject","Cohort")]
temp$dbCode <- paste("Participant",1:nrow(temp))
BAAyear3long <- merge(x=BAAyear3long, y=temp, all.x=T, by=c("Subject", "Cohort"))
prePostTime(data = BAAyear3long, xData = 'TimePoint', yData = 'ICOShiCD38hi_FreqParent', fillParam = 'Cohort', groupby = 'dbCode', title = "cTfh with Aging",
            xLabel = 'Age Group', yLabel = "ICOS+CD38+ (% cTfh)") + scale_y_continuous(breaks=seq(0,60,5),limits = c(0,55)) + 
  scale_fill_manual(values = c("#EABD81", "#7B4D9A"))  + ylab(expression(paste(paste("ICO",S^'+', "CD3",8^'+' ), " (% cTfh)")))
tukey_hsd( aov(data = BAAyear3long, formula = `ICOShiCD38hi_FreqParent` ~ `Cohort`:`TimePoint`))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/BAAyear3_FCtfh.pdf")
# write.csv(BAAyear3long[,c("dbCode","TimePoint","ICOShiCD38hi_FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1m.csv")

BAAyear3wide <- read.csv(file="D:/Pembro-Fluvac/Analysis/mergedData/BAAyear3wide.csv") 
BAAyear3wide$Cohort <- factor(BAAyear3wide$Cohort, levels=c("Young","Elderly"))

twoSampleBar(data = BAAyear3wide, xData = "Cohort", yData = "FCtfh", fillParam = "Cohort", title = "ICOS+CD38+ cTfh with aging", yLabel = "Fold-change at one week", FCplot = T) + 
  scale_y_continuous(breaks=seq(0,6,1), limits=c(0,6))+ theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())



## ******************** Subsets based on CXCR3 and CCR6 ***********************************

subsetData <- subset(mergedData, Cohort != "nonPD1" )
subsetData$cTfh_ICOShiCD38hi_CXCR3hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_CXCR3hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent * 
  subsetData$cTfh_..FreqParent * subsetData$Live_CD3..CD4..Nonnaive...FreqParent /1000000
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CXCR3hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ CXCR3hi frequencies averaged", xLabel = NULL, yLabel = "+/+ CXCR3hi (freq CD4)") 
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_CXCR3hi...FreqParent ~ Cohort + Year + TimeCategory))


subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CXCR3lo_CCR6lo...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ X3lo R6lo", xLabel = NULL, yLabel = "X3loR6lo (%ICOS+CD38+ cTfh)") 
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory != "late")
prePostTime(subsetData, xData = "TimeCategory", yData = "cTfh_ICOShiCD38hi_CXCR3lo_CCR6lo...FreqParent", fillParam = "Cohort", title = "X3loR6lo", xLabel = "TimeCategory",
            yLabel = "CXCR3lo CCR6lo (% ICOS+CD38+ cTfh)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,10), limits=c(0,60))    +
  ylab(expression(paste( "CXCR", 3^"lo", "CCR", 6^"lo", " (% ICO", S^'+', "CD3", 8^"+", " cTfh)")))
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ICOShiCD38hi_X3loR6lo_prepostTime.pdf", width=5, height=7)
# write.csv(x = subsetData[,c("dbCode", "Cohort","TimeCategory","cTfh_ICOShiCD38hi_CXCR3lo_CCR6lo...FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1o-1.csv")


temp <- reshape2::dcast(subset(subsetData, Cohort == "Healthy"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_CXCR3lo_CCR6lo...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)
temp <- reshape2::dcast(subset(subsetData, Cohort == "anti-PD-1"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_CXCR3lo_CCR6lo...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)


melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CXCR3lo_CCR6hi...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ X3lo R6hi", xLabel = NULL, yLabel = "X3loR6hi (%ICOS+CD38+ cTfh)") 
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory != "late")
prePostTime(subsetData, xData = "TimeCategory", yData = "cTfh_ICOShiCD38hi_CXCR3lo_CCR6hi...FreqParent", fillParam = "Cohort", title = "X3loR6hi", xLabel = "TimeCategory",
            yLabel = "CXCR3lo CCR6hi (% ICOS+CD38+ cTfh)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,10), limits=c(0,50))     +
  ylab(expression(paste( "CXCR", 3^"lo", "CCR", 6^"hi", " (% ICO", S^'+', "CD3", 8^"+", " cTfh)")))
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ICOShiCD38hi_X3loR6hi_prepostTime.pdf", width=5, height=7)
# write.csv(x = subsetData[,c("dbCode", "Cohort","TimeCategory","cTfh_ICOShiCD38hi_CXCR3lo_CCR6hi...FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1o-2.csv")

temp <- reshape2::dcast(subset(subsetData, Cohort == "Healthy"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_CXCR3lo_CCR6hi...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)
temp <- reshape2::dcast(subset(subsetData, Cohort == "anti-PD-1"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_CXCR3lo_CCR6hi...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)

melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_CXCR3hi_CCR6lo...FreqParent"))
prePostTimeAveraged(melted, title = "+/+ X3hi R6lo", xLabel = NULL, yLabel = "X3hiR6lo (%ICOS+CD38+ cTfh)") 
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory != "late")
prePostTime(subsetData, xData = "TimeCategory", yData = "cTfh_ICOShiCD38hi_CXCR3hi_CCR6lo...FreqParent", fillParam = "Cohort", title = "X3hiR6lo", xLabel = "TimeCategory",
            yLabel = "CXCR3hi CCR6lo (% ICOS+CD38+ cTfh)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,10), limits = c(0,70))     +
  ylab(expression(paste( "CXCR", 3^"hi", "CCR", 6^"lo", " (% ICO", S^'+', "CD3", 8^"+", " cTfh)")))
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ICOShiCD38hi_X3hiR6lo_prepostTime.pdf", width=5, height=7)
# write.csv(x = subsetData[,c("dbCode", "Cohort","TimeCategory","cTfh_ICOShiCD38hi_CXCR3hi_CCR6lo...FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1o-3.csv")

temp <- reshape2::dcast(subset(subsetData, Cohort == "Healthy"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_CXCR3hi_CCR6lo...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)
temp <- reshape2::dcast(subset(subsetData, Cohort == "anti-PD-1"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_CXCR3hi_CCR6lo...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)



## ******************** Tfr analyses ***********************************
subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_cTfh..._medfi.CD25."))
twoSampleBar(data=subset(subsetData, TimeCategory == "baseline"), xData="Cohort", yData="cTfh_ICOShiCD38hi_cTfh..._medfi.CD25.", fillParam="Cohort", 
             title="medianFI CD25 in ICOS+CD38+ cTfh at baseline", yLabel="medianFI CD25 in ICOS+CD38+ cTfh", position="left")  
summary(lm( data = subsetData, cTfh_ICOShiCD38hi_cTfh..._medfi.CD25. ~ Cohort + Year + TimeCategory))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_Foxp3hi...FreqParent"))
prePostTimeAveraged(melted, title = "Foxp3+ ICOS+CD38+ cTfh", xLabel = NULL, yLabel = "Foxp3+ (% ICOS+CD38+ cTfh)") # + scale_y_continuous(breaks=seq(0,99,0.25))
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory != "late")
prePostTime(subsetData, xData = "TimeCategory", yData = "cTfh_ICOShiCD38hi_Foxp3hi...FreqParent", fillParam = "Cohort", title = "cTfr", xLabel = "TimeCategory",
            yLabel = "Foxp3+ (% ICOS+CD38+ cTfh)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,10), limits = c(0,40))    + 
  ylab(expression(paste("Foxp",3^'+', " (%ICO",S^"+","CD3",8^"+", " cTfh)")))
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ICOShiCD38hi_Foxp3_prepostTime.pdf", width=5, height=7)
# write.csv(x = subsetData[,c("dbCode", "Cohort","TimeCategory","cTfh_ICOShiCD38hi_Foxp3hi...FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1q.csv")
temp <- reshape2::dcast(subset(subsetData, Cohort == "Healthy"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_Foxp3hi...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)
temp <- reshape2::dcast(subset(subsetData, Cohort == "anti-PD-1"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_Foxp3hi...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)

##  different denominator
subsetData$cTfh_ICOShiCD38hi_Foxp3hi...FreqParent <- subsetData$cTfh_ICOShiCD38hi_Foxp3hi...FreqParent * subsetData$cTfh_ICOShiCD38hi_..FreqParent /100
prePostTime(subsetData, xData = "TimeCategory", yData = "cTfh_ICOShiCD38hi_Foxp3hi...FreqParent", fillParam = "Cohort", title = "cTfr", xLabel = "TimeCategory",
            yLabel = "Foxp3+ICOS+CD38+ (% CXCR5+)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,1), limits=c(0,3))    + 
  ylab(expression(paste("Foxp",3^'+', "ICO",S^"+","CD3",8^"+", " (% CXCR",5^"+",")")))
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ICOShiCD38hi_Foxp3_prepostTime_freqX5.pdf", width=5, height=7)
temp <- reshape2::dcast(subset(subsetData, Cohort == "Healthy"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_Foxp3hi...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)
temp <- reshape2::dcast(subset(subsetData, Cohort == "anti-PD-1"), formula =  dbCode  ~  TimeCategory, value.var="cTfh_ICOShiCD38hi_Foxp3hi...FreqParent"  )
wilcox.test(temp$oneWeek, temp$baseline, paired=T)

# this suggests the cTfr differences are not due to global differences in Foxp3 expression but rather shifts within the CXCR5+ compartment

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("cTfh_ICOShiCD38hi_Foxp3hi...FreqParent"))
temp <- reshape2::dcast(melted, dbCode + Cohort  ~ TimeCategory, value.var = 'value'); #  temp <- temp[, -which(names(temp) == "late")]
temp$FC <- temp$oneWeek / temp$baseline
temp <- temp[-which(is.na(temp$FC)), ]
twoSampleBar(data = temp, xData = "Cohort", yData = "FC", fillParam="Cohort", title="cTfr", yLabel = "Fold-change in \nFoxp3+ ICOS+CD38+ cTfh", 
             nonparam = T) + ylab(expression(paste("Fold-change in cTfr")))  
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/ICOShiCD38hi_Foxp3_fold-change.pdf", width=5, height=8)
wilcox_test(temp, FC ~ Cohort)
# write.csv(x = temp[,c("dbCode", "Cohort","FC")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1r.csv")



##' ## ----------- PB response and correlations with cTfh response --------------------
#'
#'
## ***************    CD27++CD38++     ********************

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBar(data=subsetData, xData="Cohort", yData="CD27hiCD38hi_..FreqParent", fillParam="Cohort", title="All yrs: Plasmablast at d7", yLabel="CD19+CD27++CD38++ frequency at d7")

subsetData <- subset(mergedData, TimeCategory != "late"  & Cohort != "nonPD1")
FC_PBresponse <- dcast(subsetData, `dbCode`+`Cohort`~`TimeCategory`, value.var = c("CD27hiCD38hi_..FreqParent")); 
FC_PBresponse$FC <- FC_PBresponse$`oneWeek`/FC_PBresponse$`baseline`
twoSampleBar(data=FC_PBresponse, xData="Cohort", yData="FC", fillParam="Cohort", title="Plasmablasts", yLabel="Fold-change at one week",
             FCplot=T) +   scale_y_continuous(breaks=seq(0,99,1))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/PBresponse_foldchange_AllYrs.pdf", device="pdf", width=4.1, height=7)
# write.csv(x = FC_PBresponse, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e2b.csv")



oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
univScatter(data = oneWeek, yData = "FCPB_oW", xData = "FCtfh_oW", fillParam = "TimeCategory", 
            title = "+/+ cTfh vs CD27++CD38++ at oneWeek", yLabel= "CD27++CD38++ (freq IgD-)", xLabel = "ICOS+CD38+ cTfh frequency")
subsetData1 <- subset(oneWeek, Cohort == "Healthy" )    
subsetData2 <- subset(oneWeek, Cohort == "anti-PD-1")    
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "anti-PD-1", yData = "FCPB_oW", xData="FCtfh_oW", 
           fillParam = "Cohort", title = "Cohort 2 - One Week", yLabel= "Fold-change CD27++CD38++ PB", xLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T, 
           statsOff = F) + 
  scale_x_continuous(breaks=seq(0,99,1),limits=c(0,7)) + scale_y_continuous(breaks=seq(0,99,2),limits=c(0,15)) + 
  xlab(expression(paste("Fold-change in ICO",S^'+',"CD3",8^'+'," cTfh"))) + ylab(expression(paste("Fold-change in Plasmablasts"))) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "anti-PD-1", yData = "FCPB_oW", xData="FCtfh_oW", 
           fillParam = "Cohort", title = "Cohort 2 - One Week", yLabel= "Fold-change CD27++CD38++ PB", xLabel = "Fold-change ICOS+CD38+ cTfh", nonparam = T, 
           statsOff = T) + 
  scale_x_continuous(breaks=seq(0,99,1),limits=c(0,7)) + scale_y_continuous(breaks=seq(0,99,2),limits=c(0,15)) + 
  xlab(expression(paste("Fold-change in ICO",S^'+',"CD3",8^'+'," cTfh"))) + ylab(expression(paste("Fold-change in Plasmablasts"))) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-vs-PB_foldchange.pdf", device="pdf", width=8, height=6)
temp <- rbind(subsetData1,subsetData2); names(temp)[grep(pattern="FCPB_oW",names(temp))] <- "Fold-change in plasmablasts"; names(temp)[grep(pattern = "FCtfh_oW",names(temp))] <- "Fold-change in ICOS+CD38+ cTfh"
# write.csv(x = temp[,c("dbCode", "Cohort","Fold-change in plasmablasts","Fold-change in ICOS+CD38+ cTfh")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/1f-2.csv")



## ***************    CD71+ CD20lo ASC    ********************
baseline <- subset(mergedData, TimeCategory == "baseline")
baseline$IgDlo_CD71hi_ActBCells...FreqParent <-  baseline$IgDlo_CD71hi_ActBCells...FreqParent * baseline$IgDlo_CD71hi_..FreqParent * baseline$CD19hi_NotNaiveB_freqParent / 10000
baseline$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  baseline$IgDlo_CD71hi_CD20loCD71hi...FreqParent * baseline$IgDlo_CD71hi_..FreqParent * baseline$CD19hi_NotNaiveB_freqParent / 10000

subsetData1 <- subset(baseline, Cohort == "Healthy"   )   
subsetData2 <- subset(baseline, Cohort == "anti-PD-1"  )      
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "anti-PD-1", yData = "IgDlo_CD71hi_ActBCells...FreqParent", xData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "ABC at baseline", yLabel= "ABC (% CD19)", xLabel = "ICOS+CD38+ (% cTfh)", nonparam=T) + 
  coord_cartesian(xlim=c(0,7.5),ylim=c(0,3)) + xlab(expression(paste("ICO",S^"+","CD3",8^"+"," (% cTfh)")))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d0_vs_ABC-71hi-20hi-d0_freqCD19_biv.pdf", device="pdf")
temp <- rbind(subsetData1,subsetData2); 
# write.csv(x = temp[,c("dbCode", "Cohort","IgDlo_CD71hi_ActBCells...FreqParent","cTfh_ICOShiCD38hi_..FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4q-1.csv")


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline" & IgDlo_CD71hi_ActBCells...FreqParent > 0 & IgDlo_CD71hi_CD20loCD71hi...FreqParent >0)
twoSampleBar(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_CD20loCD71hi...FreqParent", fillParam="Cohort", title="baseline ASC", 
             yLabel="CD20- ASC (% CD71+)")+ scale_y_continuous(limits=c(0,100), breaks=seq(0,100,10))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/ASC_atBaseline_barPlot.pdf", width=4)
# sourceData exported as part of e4p.csv
# ************** different denominator *********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  oneWeek$IgDlo_CD71hi_CD20loCD71hi...FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000
subsetData1 <- subset(oneWeek, Cohort == "Healthy" )   
subsetData2 <- subset(oneWeek, Cohort == "anti-PD-1"  )    
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "anti-PD-1", yData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", xData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "cTfh vs CD20lo response at oneWeek", yLabel= "ASC (% CD19)", xLabel = "ICOS+CD38+ (% cTfh)", nonparam=T) + 
  coord_cartesian(xlim = c(0,11.2), ylim=c(0,3)) + xlab(expression(paste("ICO",S^"+","CD3",8^"+"," (% cTfh)")))+ theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())
 # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ASC-71hi-20lo-d7_freqCD19_biv.pdf", device="pdf")
temp <- rbind(subsetData1,subsetData2); 
# write.csv(x = temp[,c("dbCode", "Cohort","IgDlo_CD71hi_CD20loCD71hi...FreqParent","cTfh_ICOShiCD38hi_..FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4r-1.csv")

subsetData <- subset(mergedData, Cohort != "nonPD1" & IgDlo_CD71hi_CD20loCD71hi...FreqParent != 0 )   # exclude zero values due to presumptive technical artifact
prePostTime(subsetData, xData = "TimeCategory", yData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", fillParam = "Cohort", title = "ASC response by subject", 
            xLabel = "TimeCategory", yLabel = "CD20- ASC (% IgD-CD71+)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,150,5))     


## ***************    CD71+ CD20hi ABC    ********************
baseline <- subset(mergedData, TimeCategory == "baseline")
baseline$IgDlo_CD71hi_ActBCells...FreqParent <-  baseline$IgDlo_CD71hi_ActBCells...FreqParent * baseline$IgDlo_CD71hi_..FreqParent * baseline$CD19hi_NotNaiveB_freqParent / 10000
baseline$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  baseline$IgDlo_CD71hi_CD20loCD71hi...FreqParent * baseline$IgDlo_CD71hi_..FreqParent * baseline$CD19hi_NotNaiveB_freqParent / 10000
subsetData1 <- subset(baseline, Cohort == "Healthy"   )   
subsetData2 <- subset(baseline, Cohort == "anti-PD-1"  )      
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "anti-PD-1", yData = "IgDlo_CD71hi_CD20loCD71hi...FreqParent", xData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "ASC at baseline", yLabel= "ASC (% CD19)", xLabel = "ICOS+CD38+ (% cTfh)", nonparam=T) + 
  coord_cartesian(xlim=c(0,7.5),ylim=c(0,3)) + xlab(expression(paste("ICO",S^"+","CD3",8^"+"," (% cTfh)")))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d0_vs_ASC-71hi-20lo-d0_freqCD19_biv.pdf", device="pdf")
temp <- rbind(subsetData1,subsetData2); 
# write.csv(x = temp[,c("dbCode", "Cohort","IgDlo_CD71hi_CD20loCD71hi...FreqParent","cTfh_ICOShiCD38hi_..FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4q-2.csv")


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline" & IgDlo_CD71hi_ActBCells...FreqParent > 0 & IgDlo_CD71hi_CD20loCD71hi...FreqParent >0)
twoSampleBar(data=subsetData, xData="Cohort", yData="IgDlo_CD71hi_ActBCells...FreqParent", fillParam="Cohort", title="baseline ABC", 
             yLabel="CD20+ ABC (% CD71+)") + scale_y_continuous(limits=c(0,100), breaks=seq(0,100,10))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/ABC_atBaseline_barPlot.pdf", width=4)
# write.csv(x = subsetData[,c("Cohort", "IgDlo_CD71hi_ActBCells...FreqParent", "IgDlo_CD71hi_CD20loCD71hi...FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4p.csv")
# ************** different denominator *********************
oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$IgDlo_CD71hi_ActBCells...FreqParent <-  oneWeek$IgDlo_CD71hi_ActBCells...FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000

subsetData1 <- subset(oneWeek, Cohort == "Healthy"   )   
subsetData2 <- subset(oneWeek, Cohort == "anti-PD-1"  )      
bivScatter(data1 = subsetData1, data2 = subsetData2, name1 = "HC", name2 = "anti-PD-1", yData = "IgDlo_CD71hi_ActBCells...FreqParent", xData="cTfh_ICOShiCD38hi_..FreqParent", 
           fillParam = "Cohort", title = "cTfh vs ABC at oneWeek", yLabel= "ABC (% CD19)", xLabel = "ICOS+CD38+ (% cTfh)", nonparam=T) + 
  coord_cartesian(xlim=c(0,11.2),ylim=c(0,3)) + xlab(expression(paste("ICO",S^"+","CD3",8^"+"," (% cTfh)")))
 # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_ABC-71hi-20hi-d7_freqCD19_biv.pdf", device="pdf")
temp <- rbind(subsetData1,subsetData2); 
# write.csv(x = temp[,c("dbCode", "Cohort","IgDlo_CD71hi_CD20loCD71hi...FreqParent","cTfh_ICOShiCD38hi_..FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4r-2.csv")


subsetData <- subset(mergedData, Cohort != "nonPD1" & IgDlo_CD71hi_ActBCells...FreqParent != 0 )   # exclude zero values due to presumptive technical artifact
prePostTime(subsetData, xData = "TimeCategory", yData = "IgDlo_CD71hi_ActBCells...FreqParent", fillParam = "Cohort", title = "ABC response by subject", 
            xLabel = "TimeCategory", yLabel = "CD20+ ABC (% IgD-CD71+)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,150,5)) 


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline" & IgDlo_CD71hi_ActBCells...FreqParent > 0 & IgDlo_CD71hi_CD20loCD71hi...FreqParent >0)
twoSampleBar(data=subsetData, xData="Cohort", yData="ASC_ABCratio", fillParam="Cohort", title="baseline ASC/ABC ratio", 
             yLabel="ASC/ABC ratio") + scale_y_continuous(limits=c(0,8), breaks=seq(0,100,1))



# ************** correlations to Tfh responses *********************
data1 = mergedData[which(mergedData$Cohort == "Healthy"),]; data2 = mergedData[ which(mergedData$Cohort == "anti-PD-1") ,]
bivScatter(data1 = data1, data2 = data2, name1 = "Healthy", name2 = "anti-PD-1", xData = "Age", yData = "FCtfh_oW", fillParam = "Cohort", title = "FC Tfh vs Age", 
           xLabel = "Age", yLabel = "ICOS+CD38+ cTfh fold-change", nonparam=T) + scale_y_continuous(limits=c(0,12))

bivScatter(data1 = data1, data2 = data2, name1 = "Healthy", name2 = "anti-PD-1", xData = "Age", yData = "FCPB_oW", fillParam = "Cohort", title = "FC PB vs Age", 
           xLabel = "Age", yLabel = "CD27++CD38++ pB fold-change", nonparam=T) + scale_y_continuous(limits=c(0,12))

bivScatter(data1 = data1, data2 = data2, name1 = "Healthy", name2 = "anti-PD-1", xData = "Age", yData = "FCCXCL13_oW", fillParam = "Cohort", title = "FC CXCL13 vs Age", 
           xLabel = "Age", yLabel = "Plasma CXCL13 fold-change", nonparam=T) + scale_y_continuous(limits=c(0,5))


subsetData <- subset(mergedData, Cohort == 'Healthy' & TimeCategory == "baseline")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "FCCXCL13_oW","Age"), collapse = "|"), names(subsetData))]  
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Cohort <- "Healthy"
cor.matrix <- merge( x = cor.matrix, y = cor.matrix.pmat[,"Age"], by = "row.names"); names(cor.matrix)[grep("y",names(cor.matrix))] <- "Pvalue"

subsetData <- subset(mergedData, Cohort == 'anti-PD-1' & TimeCategory == "baseline")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "FCCXCL13_oW","Age"), collapse = "|"), names(subsetData))]  
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Cohort <- "anti-PD-1"
cor.matrix2 <- merge( x = cor.matrix2, y = cor.matrix2.pmat[,"Age"], by = "row.names"); names(cor.matrix2)[grep("^y",names(cor.matrix2))] <- "Pvalue"

temp <- as.data.frame(rbind(cor.matrix, cor.matrix2)); temp <- temp[ -grep( paste( c("Age"), collapse = "|"), temp$Row.names),]

temp$Labels <- str_replace(temp$Labels, pattern = "FCCXCL13_oW", replacement = "Fold-change \n plasma CXCL13\nat 1 week" )
temp$Labels <- str_replace(temp$Labels, pattern = "FCtfh_oW", replacement = "Fold-change \n ICOS+CD38+ cTfh\nat 1 week" )
temp$Labels <- str_replace(temp$Labels, pattern = "FCPB_oW", replacement = "Fold-change \n plasmablasts\nat 1 week" )
temp$Labels <- factor(temp$Labels, levels = c("Fold-change \n plasma CXCL13\nat 1 week" , "Fold-change \n plasmablasts\nat 1 week","Fold-change \n ICOS+CD38+ cTfh\nat 1 week"))
ggplot( data = temp, aes(x = Labels,y = Age, fill = Cohort)) + geom_bar(stat='identity',position = position_dodge(width=0.5), width=0.05, color="black", size=0.1) + 
  geom_point(aes(fill=Cohort, size=Pvalue), pch=21, color="black", stroke=0.2, position = position_dodge(width=0.5)) + 
  theme_bw() + scale_size(range = c(18,1), breaks = c(0,0.05,0.1,0.3,0.7), limits = c(0,0.8), trans = 'sqrt') + guides(size = guide_legend(reverse=TRUE)) + 
  scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + ylab("Kendall's tau vs\nAge") + xlab(" ") + ggtitle(" ") + geom_hline(yintercept=0, linetype = "dashed") +
  theme(axis.text.y = element_text(size = 28, color="black"), plot.title = element_text(size=28), axis.text.x = element_text(size=22, color="black", angle=45, hjust=1,vjust=1), 
        axis.title.x = element_text(size=28, color="black"), legend.text = element_text(size=22), legend.title = element_text(size=22)) + 
  coord_flip() + scale_y_continuous(limits = c(-1,0.5), breaks = seq(-1,1,0.25))+ theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Fold-changes_Age_correlations.pdf", device = "pdf", width=8)
# write.csv(x = temp, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e1n.csv")



cTfhCorrel <- function( data, subset )
{
  z <- vector()
  z[1] <- subset
  z[2] <- round(cor(data[ which(data$Cohort == "Healthy"), subset], data[ which(data$Cohort == "Healthy") ,"cTfh_ICOShiCD38hi_..FreqParent"], 
                    method = "kendall", use = "complete.obs"), 2)
  z[3] <- round(cor(data[ which(data$Cohort == "anti-PD-1"), subset], data[ which(data$Cohort == "anti-PD-1"),"cTfh_ICOShiCD38hi_..FreqParent"], 
                    method = "kendall", use = "complete.obs"), 2)
  return(z)
}
result <- data.frame(CellType = 0, Healthy = 0, aPD1 = 0)
result[1,] <- cTfhCorrel(baseline, subset = "IgDlo_CD71hi_ActBCells...FreqParent")
result[2,] <- cTfhCorrel(baseline, subset = "IgDlo_CD71hi_CD20loCD71hi...FreqParent")
result$CellType <- str_replace(result$CellType, "IgDlo_CD71hi_ActBCells...FreqParent", "ABC \n(%CD19)")
result$CellType <- str_replace(result$CellType, "IgDlo_CD71hi_CD20loCD71hi...FreqParent", "ASC \n(%CD19)")
result <- melt(result, id.vars = "CellType", measure.vars = c("Healthy","aPD1")); result$value <- as.numeric(result$value)
result$CellType <- factor(result$CellType, levels = c("ASC \n(%CD19)",  "ABC \n(%CD19)"))
ggplot(result, aes(x=CellType, fill=variable, color="bl")) + geom_bar(aes(y = value), stat='identity', position="dodge" ) + theme_bw()  + 
  scale_fill_manual(values=c("#7FAEDB", "#FFB18C")) + scale_color_manual(values="black") + ggtitle(expression(paste("cTfh correlation"))) + 
  ylab("Kendall t") + 
  theme(axis.text.x = element_text(angle = 0), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
        axis.title.x = element_blank(), axis.text = element_text(size=22, color="black"), legend.position = "none") + 
  scale_y_continuous(breaks = seq(-1,1,0.1), limits=c(-0.2,0.6)) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d0_vs_various-d0_pearsons_freqCD19.pdf", device = "pdf", width=4, height=6)
# write.csv(x = result, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4q-3.csv")



oneWeek <- subset(mergedData, TimeCategory == "oneWeek")
oneWeek$IgDlo_CD71hi_ActBCells...FreqParent <-  oneWeek$IgDlo_CD71hi_ActBCells...FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000
oneWeek$IgDlo_CD71hi_CD20loCD71hi...FreqParent <-  oneWeek$IgDlo_CD71hi_CD20loCD71hi...FreqParent * oneWeek$IgDlo_CD71hi_..FreqParent * oneWeek$CD19hi_NotNaiveB_freqParent / 10000

result <- data.frame(CellType = 0, Healthy = 0, aPD1 = 0)
result[1,] <- cTfhCorrel(oneWeek, subset = "IgDlo_CD71hi_ActBCells...FreqParent")
result[2,] <- cTfhCorrel(oneWeek, subset = "IgDlo_CD71hi_CD20loCD71hi...FreqParent")
result$CellType <- str_replace(result$CellType, "IgDlo_CD71hi_ActBCells...FreqParent", "ABC \n(% CD19)")
result$CellType <- str_replace(result$CellType, "IgDlo_CD71hi_CD20loCD71hi...FreqParent", "ASC \n(% CD19)")
result <- melt(result, id.vars = "CellType", measure.vars = c("Healthy","aPD1")); result$value <- as.numeric(result$value)
result$CellType <- factor(result$CellType, levels = c("ASC \n(% CD19)",  "ABC \n(% CD19)"))
ggplot(result, aes(x=CellType, fill=variable, color="bl")) + geom_bar(aes(y = value), stat='identity', position="dodge" ) + theme_bw()  + 
  scale_fill_manual(values=c("#7FAEDB", "#FFB18C")) + scale_color_manual(values="black") + ggtitle("cTfh correlation") + 
  ylab("Kendall t") + 
  theme(axis.text.x = element_text(angle = 0), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
        axis.title.x = element_blank(), axis.text = element_text(size=22, color="black"), legend.position = "none") + 
  scale_y_continuous(breaks = seq(-1,1,0.1), limits=c(-0.2,0.6))+ theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfh-d7_vs_various-d7_pearsons_freqCD19.pdf", device = "pdf", width=4, height=6)
# write.csv(x = result, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4r-3.csv")



##' ## ----------- CXCL13  --------------------
#'

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort','Year'), measure.vars = c("Plasma.CXCL13..pg.mL."))
prePostTimeAveraged(melted, title = "CXCL13", xLabel = NULL, yLabel = "Plasma CXCL13 (pg/mL)") + scale_y_continuous(breaks=seq(0,200,5), limits=c(50,110))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_TimeAveraged.pdf", device="pdf", width=4)
summary( fit <- aov(value ~ Cohort+TimeCategory + Cohort:TimeCategory, data=melted ) );  tukey_hsd(fit)
# write.csv(melted, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/1g.csv")

subsetData <- subset(melted, TimeCategory == "oneWeek")
rownames(subsetData) <- seq(1:nrow(subsetData))
twoSampleBar(data=subsetData, xData="Cohort", yData="value", fillParam="Cohort", title="CXCL13\nOne week", yLabel="plasma CXCL13 (pg/mL)") + 
  scale_y_continuous(limits = c(0,250))
aggregate( subsetData, by= list(subsetData$variable, subsetData$TimeCategory, subsetData$Cohort), FUN=mean, na.rm = T)


FC_response <- dcast( subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1" ), `Subject`+`Cohort`~`TimeCategory`, value.var = c("Plasma.CXCL13..pg.mL.")) 
FC_response$FCcxcl13 <- FC_response$oneWeek / FC_response$baseline
t.test(FC_response[which(FC_response$Cohort == "anti-PD-1"), "baseline"], FC_response[which(FC_response$Cohort == "anti-PD-1"), "oneWeek"], paired = T)
t.test(FC_response[which(FC_response$Cohort == "Healthy"), "baseline"], FC_response[which(FC_response$Cohort == "Healthy"), "oneWeek"], paired = T)
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBar(data=subsetData, xData="Cohort", yData="FCCXCL13_oW", fillParam="Cohort",FCplot =T, 
             title="CXCL13", yLabel="Fold-change at one week")+ scale_y_continuous(limits=c(0,3), breaks = seq(0,50,0.5))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_foldChange_oW.pdf", width = 4)
# write.csv(subsetData[,c("dbCode","Cohort","FCCXCL13_oW")], file="D:/Pembro-FluVac/Analysis/sourceDataExports/1h.csv")

subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "Plasma.CXCL13..pg.mL.", fillParam = "Cohort", title = "CXCL13", xLabel = "TimeCategory",
            yLabel = "Plasma CXCL13 (pg/mL)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,20))     # limits = c(0,45)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/CXCL13_prePostTime.pdf", device="pdf")
# write.csv(subsetData[,c("dbCode","Cohort", "TimeCategory","Plasma.CXCL13..pg.mL.")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e2e.csv")
univScatter(data = subset(mergedData, TimeCategory == "baseline" & Cohort == "anti-PD-1"), xData = "Plasma.CXCL13..pg.mL.",  fillParam = "dummy", position = 'right',
            yData = "FCCXCL13_oW", xLabel = "baseline CXCL13 (pg/mL)", yLabel = "Fold-change CXCL13", title = "CXCL13 in anti-PD-1") + scale_fill_manual(values=c("#FFB18C"))+
  scale_y_continuous(breaks = seq(0,5,1),limits=c(0,3))

summary(lm( Plasma.CXCL13..pg.mL. ~ Cohort + Year + TimeCategory, data=mergedData))
summary(lm( FCCXCL13_oW ~ Cohort + FCPB_oW + FChai_late, data=mergedData))


subsetData <- subset(mergedData, Cohort == 'Healthy' & TimeCategory == "baseline")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "FCCXCL13_oW","FCH1N1_hai_late","FCH3N2_hai_late","FCFluB_hai_late"), collapse = "|"), names(subsetData))]  
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Cohort <- "Healthy"
cor.matrix <- merge( x = cor.matrix, y = cor.matrix.pmat[,"FCCXCL13_oW"], by = "row.names"); names(cor.matrix)[grep("y",names(cor.matrix))] <- "Pvalue"

subsetData <- subset(mergedData, Cohort == 'anti-PD-1' & TimeCategory == "baseline")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "FCCXCL13_oW","FCH1N1_hai_late","FCH3N2_hai_late","FCFluB_hai_late"), collapse = "|"), names(subsetData))]  
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Cohort <- "anti-PD-1"
cor.matrix2 <- merge( x = cor.matrix2, y = cor.matrix2.pmat[,"FCCXCL13_oW"], by = "row.names"); names(cor.matrix2)[grep("^y",names(cor.matrix2))] <- "Pvalue"

temp <- as.data.frame(rbind(cor.matrix, cor.matrix2)); temp <- temp[ -grep( paste( c("FCCXCL13"), collapse = "|"), temp$Row.names),]
temp$Labels <- str_replace(temp$Labels, pattern = "FCtfh_oW", replacement = "Fold-change\n ICOS+CD38+ cTfh\nat 1 week" )
temp$Labels <- str_replace(temp$Labels, pattern = "FCPB_oW", replacement = "Fold-change\n plasmablasts\nat 1 week" )
temp$Labels <- str_replace(temp$Labels, pattern = "FCH1N1_hai_late", replacement = "Fold-change\n H1N1 HAI\nat late" )
temp$Labels <- str_replace(temp$Labels, pattern = "FCH3N2_hai_late", replacement = "Fold-change\n H3N2 HAI\nat late" )
temp$Labels <- str_replace(temp$Labels, pattern = "FCFluB_hai_late", replacement = "Fold-change\n FluB HAI\nat late" )
temp$Labels <- factor(temp$Labels, levels = c("Fold-change\n FluB HAI\nat late","Fold-change\n H3N2 HAI\nat late","Fold-change\n H1N1 HAI\nat late" ,
                                              "Fold-change\n plasmablasts\nat 1 week", "Fold-change\n ICOS+CD38+ cTfh\nat 1 week"))
ggplot( data = temp, aes(x = Labels,y = FCCXCL13_oW, fill = Cohort)) + geom_bar(stat='identity',position = position_dodge(width=0.5), width=0.05, color="black", size=0.1) + 
  geom_point(aes(fill=Cohort, size=Pvalue), pch=21, color="black", stroke=0.2, position = position_dodge(width=0.5)) + 
  theme_bw() + scale_size(range = c(18,1), breaks = c(0,0.05,0.1,0.3,0.7), limits = c(0,0.8), trans = 'sqrt') + guides(size = guide_legend(reverse=TRUE)) + 
  scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + ylab("Kendall's tau vs\nFold-change in CXCL13") + xlab(" ") + ggtitle(" ") + geom_hline(yintercept=0, linetype = "dashed") +
  theme(axis.text.y = element_text(size = 28, color="black"), plot.title = element_text(size=28), axis.text.x = element_text(size=22, color="black", angle=45, hjust=1,vjust=1), 
        axis.title.x = element_text(size=28, color="black"), legend.text = element_text(size=22), legend.title = element_text(size=22)) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())+ 
  coord_flip() + scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.25))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Fold-changes_CXCL13_correlations.pdf", device = "pdf", width=9, height=10)
# write.csv(temp, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e2h.csv")


##' ## ----------- Vaccine response determinants and biomarkers --------------------
#'

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
  a<-twoSampleBar(data=subsetData, xData="Cohort", yData="H1N1pdm09.HAI.titer", fillParam="Cohort", batch="none", title="Baseline - H1N1", yLabel="H1N1 titer", nonparam=T) +
  scale_y_continuous(trans='log2', breaks=c(2^(2:14) ) ) +  coord_cartesian(ylim = c(4,8192))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAItiter_late_byCohort.pdf", device="pdf", width=4.5)
  b<-twoSampleBar(data=subsetData, xData="Cohort", yData="H3N2.HAI.titer", fillParam="Cohort", batch="none", title="Baseline - H3N2", yLabel="H3N2 titer", nonparam=T) +
  scale_y_continuous(trans='log2', breaks=c(2^(2:14) ) ) +  coord_cartesian(ylim = c(4,8192))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAItiter_late_byCohort.pdf", device="pdf", width=4.5)
  c<-twoSampleBar(data=subsetData, xData="Cohort", yData="FluB.HAI.titer", fillParam="Cohort", batch="none", title="Baseline - FluB", yLabel="FluB titer", nonparam=T) +
  scale_y_continuous(trans='log2', breaks=c(2^(2:14) ) ) +  coord_cartesian(ylim = c(4,8192))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAItiter_late_byCohort.pdf", device="pdf", width=4.5)
grid.arrange(a,b,c,nrow=1)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAI_baselinetiters_3strains.pdf", arrangeGrob(a,b,c,nrow=1), width=14, height=8)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","H1N1pdm09.HAI.titer", "H3N2.HAI.titer", "FluB.HAI.titer")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e3b.csv")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
subsetData %>% group_by(Cohort) %>% get_summary_stats(FCH1N1_hai_late, FCH3N2_hai_late, FCFluB_hai_late)

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
  a<-twoSampleBar(data=subsetData, xData="Cohort", yData="FCH1N1_hai_late", fillParam="Cohort", title="H1N1", yLabel="Fold-change in HAI", nonparam = T) + 
  scale_y_continuous(trans='log2',breaks=c(2^(0:14)), limits = c(1,175))
  b<-twoSampleBar(data=subsetData, xData="Cohort", yData="FCH3N2_hai_late", fillParam="Cohort", title="H3N2", yLabel="Fold-change in HAI", nonparam = T) + 
  scale_y_continuous(trans='log2',breaks=c(2^(0:14)), limits = c(1,175))
  c<-twoSampleBar(data=subsetData, xData="Cohort", yData="FCFluB_hai_late", fillParam="Cohort", title="Influenza B", yLabel="Fold-change in HAI", nonparam = T) + 
  scale_y_continuous(trans='log2',breaks=c(2^(0:14)), limits = c(1,175))
grid.arrange(a,b,c,nrow=1)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Seroconversion_3strains.pdf", arrangeGrob(a,b,c,nrow=1), width=11)
# write.csv(subsetData[,c("dbCode","Cohort","FCH1N1_hai_late", "FCH3N2_hai_late", "FCFluB_hai_late")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2a.csv")


plotSeroProtection <- function(STRAIN, PLOTTITLE)
{
  subsetData <- subset(mergedData, Cohort != "nonPD1" )
  melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = STRAIN); melted <- melted[ which( !is.na(melted$value) ),]
  overTime <- aggregate( melted$value, by= list(melted$variable, melted$TimeCategory, melted$Cohort), FUN=mean, na.rm = T); overTime
  aPD1_seroprot <- length(which(melted$value > 40 & melted$Cohort == "anti-PD-1" & melted$TimeCategory == "late") )
  aPD1_total <- length(which(melted$Cohort == "anti-PD-1" & melted$TimeCategory == "late") )
  HC_seroprot <- length(which(melted$value > 40 & melted$Cohort == "Healthy" & melted$TimeCategory == "late") )
  HC_total <- length(which(melted$Cohort == "Healthy" & melted$TimeCategory == "late") )
  continTable <- matrix(c(aPD1_seroprot, aPD1_total - aPD1_seroprot, HC_seroprot, HC_total - HC_seroprot), ncol=2); continTable
  print( fisher.test(continTable))
  
  seroprot <- data.frame( Cohort = c("Healthy", "anti-PD-1"), Seroprot = c(HC_seroprot / HC_total, aPD1_seroprot / aPD1_total))
  seroprot$Cohort <- factor(seroprot$Cohort, levels = c("Healthy", "anti-PD-1"))
  return( ggplot(data=seroprot, aes(x=Cohort, y=Seroprot, fill=Cohort,width=0.6)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) +  
    geom_bar( position = position_dodge(), stat = "identity", color="black",size=0.1) + ggtitle(PLOTTITLE) + ylab("Proportion seroprotected") +  theme_bw() +
    theme(axis.text = element_text(size=28,hjust = 0.5,color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=32,hjust = 0.5)) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=45,hjust=1,vjust=1)) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) +
    coord_cartesian(ylim = c(0,1)) + scale_y_continuous(breaks = seq(0,1,0.1) )    ) 
}

  a <- plotSeroProtection(STRAIN = "H1N1pdm09.HAI.titer", PLOTTITLE="H1N1"); a
  b <- plotSeroProtection(STRAIN = "H3N2.HAI.titer", PLOTTITLE="H3N2"); b
  c <- plotSeroProtection(STRAIN = "FluB.HAI.titer", PLOTTITLE="Influenza B"); c
grid.arrange(a,b,c,nrow=1)
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Seroprotection_3strains.pdf", arrangeGrob(a,b,c,nrow=1), width=11)
  # write.csv(a$data, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2b-1.csv")
  # write.csv(b$data, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2b-2.csv")
  # write.csv(c$data, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2b-3.csv")

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = "H1N1pdm09.HAI.titer"); melted <- melted[ which( !is.na(melted$value) ),]
  a<-prePostTimeAveraged(melted, title = "H1N1", xLabel = NULL, yLabel = "H1N1 HAI titer") + scale_y_continuous(trans = 'log2', breaks=c(2^(1:14)), limits = c(2,2500)) +
  geom_hline(yintercept = 40, linetype = "dashed", color="grey50") + annotate("text",x = 1.5, y=30, label="Seroprotection", size=7,color="grey50")
Anova(fit <- lm( value ~ Cohort * TimeCategory, data = melted)); tukey_hsd(fit)
write.csv(melted, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2c-1.csv")
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = "H3N2.HAI.titer"); melted <- melted[ which( !is.na(melted$value) ),]
  b<-prePostTimeAveraged(melted, title = "H3N2", xLabel = NULL, yLabel = "H3N2 HAI titer") + scale_y_continuous(trans = 'log2', breaks=c(2^(1:14)), limits = c(2,2500)) +
  geom_hline(yintercept = 40, linetype = "dashed", color="grey50") + annotate("text",x = 1.5, y=30, label="Seroprotection", size=7,color="grey50")
Anova(fit <- lm( value ~ Cohort * TimeCategory, data = melted)); tukey_hsd(fit)
# write.csv(melted, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2c-2.csv")
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = "FluB.HAI.titer"); melted <- melted[ which( !is.na(melted$value) ),]
  c<-prePostTimeAveraged(melted, title = "Influenza B", xLabel = NULL, yLabel = "Influenza B HAI titer") + scale_y_continuous(trans = 'log2', breaks=c(2^(1:14)), limits = c(2,2500)) +
  geom_hline(yintercept = 40, linetype = "dashed", color="grey50") + annotate("text",x = 1.5, y=30, label="Seroprotection", size=7,color="grey50")
Anova(fit <- lm( value ~ Cohort * TimeCategory, data = melted)); tukey_hsd(fit)
# write.csv(melted, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2c-3.csv")
grid.arrange(a,b,c,nrow=1)
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAIresponsesTimeAveraged_3strains.pdf", arrangeGrob(a,b,c,nrow=1), width=14)



subsetData <- subset(mergedData, Cohort != "nonPD1" )
  a<-prePostTime(data = subsetData, xData = "TimeCategory", yData = "H1N1pdm09.HAI.titer", fillParam = "Cohort", title = "H1N1", xLabel = "TimeCategory",
            yLabel = "H1N1pdm09 HAI titer", groupby = "dbCode") + scale_y_continuous(trans='log2', breaks=2^(1:13))
  b<-prePostTime(data = subsetData, xData = "TimeCategory", yData = "H3N2.HAI.titer", fillParam = "Cohort", title = "H3N2", xLabel = "TimeCategory",
               yLabel = "H3N2 HAI titer", groupby = "dbCode") + scale_y_continuous(trans='log2', breaks=2^(1:13))
  c<-prePostTime(data = subsetData, xData = "TimeCategory", yData = "FluB.HAI.titer", fillParam = "Cohort", title = "FluB", xLabel = "TimeCategory",
               yLabel = "FluB HAI titer", groupby = "dbCode") + scale_y_continuous(trans='log2', breaks=2^(1:13))
grid.arrange(a,b,c,nrow=1)
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/HAItiter_overTime_3strains.pdf", arrangeGrob(a,b,c,nrow=1), width=14)
# write.csv(a$data, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e3a-1.csv")
# write.csv(b$data, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e3a-2.csv")
# write.csv(c$data, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e3a-3.csv")


###' ## ------------ Correlations with HAI titer
subsetData <- subset(mergedData, Cohort == 'Healthy' & TimeCategory == "oneWeek")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "cTfh_ICOShiCD38hi_..FreqParent", "CD27hiCD38hi_..FreqParent", "FCH1N1_hai_late"), collapse = "|"), names(subsetData))]  
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Cohort <- "Healthy"
cor.matrix <- merge( x = cor.matrix, y = cor.matrix.pmat[,"FCH1N1_hai_late"], by = "row.names"); names(cor.matrix)[grep("y",names(cor.matrix))] <- "Pvalue"

subsetData <- subset(mergedData, Cohort == 'anti-PD-1'& TimeCategory == "oneWeek")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "cTfh_ICOShiCD38hi_..FreqParent", "CD27hiCD38hi_..FreqParent","FCH1N1_hai_late"), collapse = "|"), names(subsetData))]  
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Cohort <- "anti-PD-1"
cor.matrix2 <- merge( x = cor.matrix2, y = cor.matrix2.pmat[,"FCH1N1_hai_late"], by = "row.names"); names(cor.matrix2)[grep("^y",names(cor.matrix2))] <- "Pvalue"

temp <- as.data.frame(rbind(cor.matrix, cor.matrix2)); temp <- temp[ -grep( "FCH1N1", temp$Row.names),]

temp$Labels <- str_replace(temp$Labels, pattern = "FCtfh_oW", replacement = "Fold-change ICOS+CD38+ cTfh\nat one week" )
temp$Labels <- str_replace(temp$Labels, pattern = "FCPB_oW", replacement = "Fold-change plasmablasts\nat one week" )
temp$Labels <- str_replace(temp$Labels, pattern = "cTfh_ICOShiCD38hi_..FreqParent", replacement = "ICOS+CD38+ cTfh (%CXCR5+CD4+)\nat one week" )
temp$Labels <- str_replace(temp$Labels, pattern = "CD27hiCD38hi_..FreqParent", replacement = "Plasmablasts (% CD19) \nat one week" )
temp$Labels <- factor(temp$Labels, levels = c( "Plasmablasts (% CD19) \nat one week","Fold-change plasmablasts\nat one week", 
                                               "ICOS+CD38+ cTfh (%CXCR5+CD4+)\nat one week","Fold-change ICOS+CD38+ cTfh\nat one week"))
ggplot( data = temp, aes(x = Labels,y = FCH1N1_hai_late, fill = Cohort)) + geom_bar(stat='identity',position = position_dodge(width=0.5), width=0.05, color="black", size=0.1) + 
  geom_point(aes(fill=Cohort, size=Pvalue), pch=21, color="black", stroke=0.2, position = position_dodge(width=0.5)) + 
  theme_bw() + scale_size(range = c(18,1), breaks = c(0,0.05,0.1,0.3,0.7), limits = c(0,0.8), trans = 'sqrt') + guides(size = guide_legend(reverse=TRUE)) + 
  scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + ylab("Kendall's tau vs\nFold-change in H1N1 HAI") + xlab(" ") + ggtitle("H1N1") + geom_hline(yintercept=0, linetype = "dashed") +
  theme(axis.text.y = element_text(size = 28, color="black"), plot.title = element_text(size=28), axis.text.x = element_text(size=22, color="black", angle=45, hjust=1,vjust=1), 
        axis.title.x = element_text(size=28, color="black"), legend.text = element_text(size=22), legend.title = element_text(size=22)) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) + 
  coord_flip() + scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.25)) 
 # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Fold-changes_HAI_correlations.pdf", device = "pdf", width=13, height=9)
# write.csv(temp, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e3e-1.csv")

subsetData <- subset(mergedData, Cohort == 'Healthy' & TimeCategory == "oneWeek")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "cTfh_ICOShiCD38hi_..FreqParent", "CD27hiCD38hi_..FreqParent", "FCH3N2_hai_late"), collapse = "|"), names(subsetData))]  
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Cohort <- "Healthy"
cor.matrix <- merge( x = cor.matrix, y = cor.matrix.pmat[,"FCH3N2_hai_late"], by = "row.names"); names(cor.matrix)[grep("y",names(cor.matrix))] <- "Pvalue"

subsetData <- subset(mergedData, Cohort == 'anti-PD-1'& TimeCategory == "oneWeek")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "cTfh_ICOShiCD38hi_..FreqParent", "CD27hiCD38hi_..FreqParent","FCH3N2_hai_late"), collapse = "|"), names(subsetData))]  
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Cohort <- "anti-PD-1"
cor.matrix2 <- merge( x = cor.matrix2, y = cor.matrix2.pmat[,"FCH3N2_hai_late"], by = "row.names"); names(cor.matrix2)[grep("^y",names(cor.matrix2))] <- "Pvalue"

temp <- as.data.frame(rbind(cor.matrix, cor.matrix2)); temp <- temp[ -grep( "FCH3N2", temp$Row.names),]

temp$Labels <- str_replace(temp$Labels, pattern = "FCtfh_oW", replacement = "Fold-change ICOS+CD38+ cTfh\nat one week" )
temp$Labels <- str_replace(temp$Labels, pattern = "FCPB_oW", replacement = "Fold-change plasmablasts\nat one week" )
temp$Labels <- str_replace(temp$Labels, pattern = "cTfh_ICOShiCD38hi_..FreqParent", replacement = "ICOS+CD38+ cTfh (%CXCR5+CD4+)\nat one week" )
temp$Labels <- str_replace(temp$Labels, pattern = "CD27hiCD38hi_..FreqParent", replacement = "Plasmablasts (% CD19) \nat one week" )
temp$Labels <- factor(temp$Labels, levels = c( "Plasmablasts (% CD19) \nat one week","Fold-change plasmablasts\nat one week", 
                                              "ICOS+CD38+ cTfh (%CXCR5+CD4+)\nat one week","Fold-change ICOS+CD38+ cTfh\nat one week"))
ggplot( data = temp, aes(x = Labels,y = FCH3N2_hai_late, fill = Cohort)) + geom_bar(stat='identity',position = position_dodge(width=0.5), width=0.05, color="black", size=0.1) + 
 geom_point(aes(fill=Cohort, size=Pvalue), pch=21, color="black", stroke=0.2, position = position_dodge(width=0.5)) + 
 theme_bw() + scale_size(range = c(18,1), breaks = c(0,0.05,0.1,0.3,0.7), limits = c(0,0.8), trans = 'sqrt') + guides(size = guide_legend(reverse=TRUE)) + 
 scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + ylab("Kendall's tau vs\nFold-change in H3N2 HAI") + xlab(" ") + ggtitle("H3N2") + geom_hline(yintercept=0, linetype = "dashed") +
 theme(axis.text.y = element_text(size = 28, color="black"), plot.title = element_text(size=28), axis.text.x = element_text(size=22, color="black", angle=45, hjust=1,vjust=1), 
       axis.title.x = element_text(size=28, color="black"), legend.text = element_text(size=22), legend.title = element_text(size=22))+ theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) + 
 coord_flip() + scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.25))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Fold-changes_H3N2_correlations.pdf", device = "pdf", width=13, height=9)
# write.csv(temp, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e3e-2.csv")

 
subsetData <- subset(mergedData, Cohort == 'Healthy' & TimeCategory == "oneWeek")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "cTfh_ICOShiCD38hi_..FreqParent", "CD27hiCD38hi_..FreqParent", "FCFluB_hai_late"), collapse = "|"), names(subsetData))]  
cor.matrix <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix <- as.data.frame(cor.matrix); cor.matrix$Labels <- row.names(cor.matrix); cor.matrix$Cohort <- "Healthy"
cor.matrix <- merge( x = cor.matrix, y = cor.matrix.pmat[,"FCFluB_hai_late"], by = "row.names"); names(cor.matrix)[grep("y",names(cor.matrix))] <- "Pvalue"

subsetData <- subset(mergedData, Cohort == 'anti-PD-1'& TimeCategory == "oneWeek")
subsetData <- subsetData[,grep( paste(c("FCtfh_oW","FCPB_oW", "cTfh_ICOShiCD38hi_..FreqParent", "CD27hiCD38hi_..FreqParent","FCFluB_hai_late"), collapse = "|"), names(subsetData))]  
cor.matrix2 <- cor(subsetData, method="kendall" , use="pairwise.complete.obs" )
cor.matrix2.pmat <- ggcorrplot::cor_pmat(subsetData, method="kendall", use="pairwise.complete.obs"  )
cor.matrix2 <- as.data.frame(cor.matrix2); cor.matrix2$Labels <- row.names(cor.matrix2); cor.matrix2$Cohort <- "anti-PD-1"
cor.matrix2 <- merge( x = cor.matrix2, y = cor.matrix2.pmat[,"FCFluB_hai_late"], by = "row.names"); names(cor.matrix2)[grep("^y",names(cor.matrix2))] <- "Pvalue"

temp <- as.data.frame(rbind(cor.matrix, cor.matrix2)); temp <- temp[ -grep( "FCFluB", temp$Row.names),]

temp$Labels <- str_replace(temp$Labels, pattern = "FCtfh_oW", replacement = "Fold-change ICOS+CD38+ cTfh\nat one week" )
temp$Labels <- str_replace(temp$Labels, pattern = "FCPB_oW", replacement = "Fold-change plasmablasts\nat one week" )
temp$Labels <- str_replace(temp$Labels, pattern = "cTfh_ICOShiCD38hi_..FreqParent", replacement = "ICOS+CD38+ cTfh (%CXCR5+CD4+)\nat one week" )
temp$Labels <- str_replace(temp$Labels, pattern = "CD27hiCD38hi_..FreqParent", replacement = "Plasmablasts (% CD19) \nat one week" )
temp$Labels <- factor(temp$Labels, levels = c( "Plasmablasts (% CD19) \nat one week","Fold-change plasmablasts\nat one week", 
                                               "ICOS+CD38+ cTfh (%CXCR5+CD4+)\nat one week","Fold-change ICOS+CD38+ cTfh\nat one week"))
ggplot( data = temp, aes(x = Labels,y = FCFluB_hai_late, fill = Cohort)) + geom_bar(stat='identity',position = position_dodge(width=0.5), width=0.05, color="black", size=0.1) + 
  geom_point(aes(fill=Cohort, size=Pvalue), pch=21, color="black", stroke=0.2, position = position_dodge(width=0.5)) + 
  theme_bw() + scale_size(range = c(18,1), breaks = c(0,0.05,0.1,0.3,0.7), limits = c(0,0.8), trans = 'sqrt') + guides(size = guide_legend(reverse=TRUE)) + 
  scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + ylab("Kendall's tau vs\nFold-change in FluB HAI") + xlab(" ") + ggtitle("Influenza B") + geom_hline(yintercept=0, linetype = "dashed") +
  theme(axis.text.y = element_text(size = 28, color="black"), plot.title = element_text(size=28), axis.text.x = element_text(size=22, color="black", angle=45, hjust=1,vjust=1), 
        axis.title.x = element_text(size=28, color="black"), legend.text = element_text(size=22), legend.title = element_text(size=22))+ theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) + 
  coord_flip() + scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.25))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Fold-changes_FluB_correlations.pdf", device = "pdf", width=13, height=9)
# write.csv(temp, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e3e-3.csv")

 
##' ## ----------- glycosylation --------------------




# ***************************    Total afucosylation  *******************************


# nonspecific afucosylation
subsetData <- subset(mergedData, Cohort == "anti-PD-1" )
univScatter(data = subsetData, xData = "IgG1_Total.afucosylated", yData="Nonsp_afucosylated",fillParam = 'dummy', title = "Nonsp vs spec afucosylation", 
            xLabel = "anti-H1 afucosylation", yLabel = "Nonspecific afucosylation") + scale_x_continuous(limits=c(0,16))+ scale_fill_manual(values=c("#FFB18C"))

# H1 specific afucosylation
subsetData <- subset(mergedData, Cohort != "nonPD1" )

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.afucosylated", fillParam = "Cohort", title = "IgG1 afucosylation", xLabel = "TimeCategory",
            yLabel = "anti-H1 Afucosylation (% abund)", groupby = "dbCode")  + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1afucosylation_overTime.pdf", width=6)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.afucosylated")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4b.csv")

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.afucosylated", fillParam = "Cohort", title = "IgG2 afucosylation over time", xLabel = "TimeCategory",
            yLabel = "Afucosylation (% abundance, anti-H1)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30))
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.afucosylated", fillParam = "Cohort", title = "IgG3/4 afucosylation over time", xLabel = "TimeCategory",
            yLabel = "Afucosylation (% abundance, anti-H1)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30))

# ***************************    Total Bisection   *******************************
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.Bisection", fillParam = "Cohort", title = "IgG1 bisection", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG1", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.Bisection", fillParam = "Cohort", title = "IgG2 bisection over time", xLabel = "TimeCategory",
            yLabel = "Bisection (% abundance, anti-H1)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.Bisection", fillParam = "Cohort", title = "IgG3/4 bisection over time", xLabel = "TimeCategory",
            yLabel = "Bisection (% abundance, anti-H1)", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )

twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), xData="Cohort", yData="IgG1_Total.Bisection", fillParam="Cohort", 
             title="Baseline bisection", yLabel="% anti-H1 IgG1")+ scale_y_continuous(limits = c(0,20), breaks = seq(0,50,5))
twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek"), xData="Cohort", yData="IgG1_Total.Bisection", fillParam="Cohort", 
             title="One Week bisection", yLabel="% anti-H1 IgG1")+ scale_y_continuous(limits = c(0,20), breaks = seq(0,50,5))


# ***************************    Galactosylation *******************************

# nonspecific afucosylation
subsetData <- subset(mergedData, Cohort == "anti-PD-1" )
univScatter(data = subsetData, xData = "IgG1_Total.Galactosylation..G1.G2.", yData="Nonsp_GalactosylationG1.G2",fillParam = 'dummy', title = "Nonsp vs spec galactosylation", 
            xLabel = "anti-H1 galactosylation", yLabel = "Nonspecific galactosylation") + scale_x_continuous(limits=c(90,100))+ scale_fill_manual(values=c("#FFB18C"))

# specific galactosylation
subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort','Year'), measure.vars = c("IgG1_Total.Galactosylation..G1.G2."))
prePostTimeAveraged(melted, title = "Total galactosylation", xLabel = NULL, yLabel = "% anti-H1 IgG1") + 
  scale_y_continuous(breaks=seq(0,100,1), limits = c(90,100))
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_totGalactosylation_overTime.pdf", width=4.5, height=7)
# melted %>%  anova_test(value ~ TimeCategory+Cohort )                        # two-way anova without interaction term 
Anova(fit <- lm( value ~ Cohort * TimeCategory, data = melted)); tukey_hsd(fit)
# write.csv(melted, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2d.csv")

subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.Galactosylation..G1.G2.", fillParam = "Cohort", title = "IgG1 total galactosylation", 
            xLabel = "TimeCategory", yLabel = "% anti-H1 IgG1", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) ) + 
  coord_cartesian(ylim = c(90,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_Tot-galactosylation_overTime_perSubject.pdf", width=6.5)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.afucosylated")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4f-1.csv")

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.Galactosylation..G1.G2.", fillParam = "Cohort", title = "IgG2 total galactosylation", 
            xLabel = "TimeCategory", yLabel = "% anti-H1 IgG2", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )+ 
  coord_cartesian(ylim = c(80,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG2_Tot-galactosylation_overTime_perSubject.pdf", width=6.5)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.afucosylated")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4f-2.csv")

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.Galactosylation..G1.G2.", fillParam = "Cohort", title = "IgG3/4 total galactosylation", 
            xLabel = "TimeCategory", yLabel = "% anti-H1 IgG3/4", groupby = "dbCode") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )+ 
  coord_cartesian(ylim = c(50,100))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG34_Tot-galactosylation_overTime_perSubject.pdf", width=6.5)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.afucosylated")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4f-3.csv")

subsetData <- subset(mergedData, Cohort != "nonPD1"  & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="IgG1_Total.Galactosylation..G1.G2.", fillParam="Cohort", title="Total galactosylation\nbaseline", 
             yLabel="% anti-H1 IgG1") + 
  scale_y_continuous(breaks = seq(0,100,2)) +   coord_cartesian(ylim = c(88,102)) 
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_totGalactosylation_bL.pdf", width=5, height=8)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.Galactosylation..G1.G2.")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2e.csv")

prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.G0", fillParam = "Cohort", title = "IgG1 G0 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG1", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.G0", fillParam = "Cohort", title = "IgG2 G0 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG2", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.G0", fillParam = "Cohort", title = "IgG3/4 G0 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG3/4", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,30) )

subsetData <- subset(mergedData, Cohort != "nonPD1"  & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="IgG1_Total.G0", fillParam="Cohort", 
             title="G0 - baseline", yLabel="% anti-H1 IgG1", nonparam=T) + scale_y_continuous(breaks = seq(0,100,2)) + 
  coord_cartesian(ylim = c(0,10))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_G0_Galactosylation_bL.pdf", width=4.5)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.G0")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4e-1.csv")



prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.G1", fillParam = "Cohort", title = "IgG1 G1 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G1 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.G1", fillParam = "Cohort", title = "IgG2 G1 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G1 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.G1", fillParam = "Cohort", title = "IgG3/4 G1 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G1 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,100) )

subsetData <- subset(mergedData, Cohort != "nonPD1"  & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="IgG1_Total.G1", fillParam="Cohort", 
             title="G1 - baseline", yLabel="% anti-H1 IgG1", nonparam=T) + scale_y_continuous(limits = c(0,70), breaks = seq(0,100,10)) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_G1_Galactosylation_bL.pdf", width=4.5)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.G1")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4e-2.csv")


prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.G2", fillParam = "Cohort", title = "IgG1 G2 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G2 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,80) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.G2", fillParam = "Cohort", title = "IgG2 G2 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G2 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,80) )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.G2", fillParam = "Cohort", title = "IgG3/4 G2 galactosylation over time", xLabel = "TimeCategory",
            yLabel = "G2 Galactosylation (% abundance, anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,5), limits = c(0,80) )

subsetData <- subset(mergedData, Cohort != "nonPD1"  & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="IgG1_Total.G2", fillParam="Cohort", 
             title="G2 - baseline", yLabel="% anti-H1 IgG1", nonparam=T) + scale_y_continuous(limits = c(0,80), breaks = seq(0,100,10)) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_G2_Galactosylation_bL.pdf", width=4.5)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.G2")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4e-3.csv")




# ***************************    Sialylation *******************************

subsetData1 <- subset(mergedData, TimeCategory == "baseline" & Cohort == "Healthy")
  subsetData2 <- subset(mergedData, TimeCategory == "baseline" & Cohort == "anti-PD-1")
bivScatter(data1 = subsetData1, data2 = subsetData2,
           name1 = "Healthy", name2 = "anti-PD-1", xData = "IgG1_Total.Galactosylation..G1.G2.", yData = "IgG1_Total.sialylated", fillParam = "Cohort", 
           title = "Sialylation vs Galactosylation", xLabel = "Total galactosylation (% anti-H1 IgG)", yLabel = "Sialylation (% anti-H1 IgG)", nonparam = T) + 
  coord_cartesian(ylim = c(0,30), xlim = c(90,99))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_sial-vs-TotGalactosyl_correlation_twoCohorts.pdf", width=8, height=8)
temp <- rbind(subsetData1,subsetData2); 
# write.csv(x = temp[,c("dbCode", "Cohort","IgG1_Total.Galactosylation..G1.G2.","IgG1_Total.sialylated")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4c.csv")


# nonspecific sialylation
subsetData <- subset(mergedData, Cohort == "anti-PD-1" )
univScatter(data = subsetData, xData = "IgG1_Total.sialylated", yData="Nonsp_sialylated",fillParam = 'dummy', title = "Nonsp vs spec sialylation", 
            xLabel = "anti-H1 sialylation", yLabel = "Nonspecific sialylation") + scale_x_continuous(limits=c(5,30))+ scale_fill_manual(values=c("#FFB18C"))
  

# specific sialylation
subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort','Year'), measure.vars = c("IgG1_Total.sialylated"))
prePostTimeAveraged(melted, title = "Sialylation", xLabel = NULL, yLabel = "% anti-H1 IgG1") + scale_y_continuous(breaks=seq(0,99,2), limits = c(10,22))
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_sialylation_overTime.pdf", width=4.5, height=7)
# write.csv(melted, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2f.csv")

melted %>%  anova_test(value ~ TimeCategory+Cohort )                        # two-way anova without interaction term 
Anova(fit <- lm( value ~ Cohort * TimeCategory, data = melted)); tukey_hsd(fit)
# write.csv(melted, file = "./sialylation_forTwoWayAnova.csv")

subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
twoSampleBar(data=subsetData, xData="Cohort", yData="IgG1_Total.sialylated", fillParam="Cohort", 
             title="Sialylation\nbaseline", yLabel="% anti-H1 IgG1") + scale_y_continuous(limits = c(0,40), breaks = seq(0,50,5))
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1sialylation_bL.pdf", width=5, height=8)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.sialylated")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2g.csv")


subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek")
twoSampleBar(data=subsetData, xData="Cohort", yData="IgG1_Total.sialylated", fillParam="Cohort", 
             title="One week", yLabel="Sialylation (% anti-H1 IgG1)")+ scale_y_continuous(limits = c(0,40), breaks = seq(0,50,5))

subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG1_Total.sialylated", fillParam = "Cohort", title = "IgG1 sialylation", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG1", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,10), limits = c(0,40) )
 # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG1_sialylation_overTime_perSubject.pdf", width=6)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG1_Total.sialylated")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4f-1.csv")
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG2_Total.sialylated", fillParam = "Cohort", title = "IgG2 sialylation", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG2", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,10), limits = c(0,40) )
 # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG2_sialylation_overTime_perSubject.pdf", width=6)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG2_Total.sialylated")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4f-2.csv")
prePostTime(subsetData, xData = "TimeCategory", yData = "IgG34_Total.sialylated", fillParam = "Cohort", title = "IgG3/4 sialylation", xLabel = "TimeCategory",
            yLabel = "% anti-H1 IgG3/4", groupby = "Subject") + scale_y_continuous(breaks=seq(0,240,10), limits = c(0,40) )
 # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/IgG34_sialylation_overTime_perSubject.pdf", width=6)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","IgG34_Total.sialylated")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4f-3.csv")


summary(lm( FChai_late ~ Cohort + Year + IgG1sial_oW, data=subsetData))
summary(lm( IgG1sial_oW ~ Cohort + FCtfh_oW, data=subsetData))



 
##' ## -----------    Affinity   -------------------------
subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"); subsetData$X1.Kd.HA <- 1 / subsetData$X1.Kd.HA
twoSampleBar(data=subsetData, xData="Cohort", yData="X1.Kd.HA", fillParam="Cohort", title="HA EC50/Kd\nBaseline", yLabel="1 / Concentration (ug/mL)",position = 'left') + 
  scale_y_continuous(breaks = seq(0,100,1)) + theme(plot.title = element_text(size=32)) +
  ylab(expression(paste("1 / H1 Concentration (",mu,"g/mL)"))) + ggtitle(expression(paste("E",C[50],"/",K[d], " baseline")))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/AntibodyAffinity_baseLine.pdf", device="pdf", width=4.1, height=8)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","X1.Kd.HA")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4k.csv")


subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('Subject', 'TimeCategory', 'Cohort','Year'), measure.vars = c("Lo.Hi.HA.affinity"))
prePostTimeAveraged(melted, title = "anti-HA Affinity", xLabel = NULL, yLabel = "Lo-Hi HA affinity") + scale_y_continuous(breaks=seq(0,99,0.03), limits = c(0.6,0.9))
 # ggsave (filename = "D:/Pembro-Fluvac/Analysis/Images/antiHA-affinity_overTime_LoHi.pdf", width=4.5, height=7)
melted %>%  anova_test(value ~ TimeCategory+Cohort )                        # two-way anova without interaction term 
Anova(fit <- lm( value ~ Cohort * TimeCategory, data = melted)); tukey_hsd(fit)
# write.csv(melted, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4l.csv")

subsetData <- subset(mergedData, Cohort != "nonPD1" )
prePostTime(subsetData, xData = "TimeCategory", yData = "Lo.Hi.HA.affinity", fillParam = "Cohort", title = "anti-HA Affinity", xLabel = "TimeCategory",
            yLabel = "Lo-Hi HA affinity (anti-H1)", groupby = "Subject") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0,1) )
 # ggsave (filename = "D:/Pembro-Fluvac/Analysis/Images/antiHA-affinity_overTime_LoHi_perSubject.pdf", width=8)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","Lo.Hi.HA.affinity")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e4m.csv")

twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline"), xData="Cohort", yData="Lo.Hi.HA.affinity", fillParam="Cohort", 
             title="anti-HA affinity at baseline", yLabel="Lo Hi HA affinity (anti-H1)") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0,1.1) )
twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "oneWeek"), xData="Cohort", yData="Lo.Hi.HA.affinity", fillParam="Cohort", 
             title="anti-HA affinity at oneWeek", yLabel="Lo Hi HA affinity (anti-H1)") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0,1.1) )
twoSampleBar(data=subset(mergedData, Cohort != "nonPD1" & TimeCategory == "late"), xData="Cohort", yData="Lo.Hi.HA.affinity", fillParam="Cohort", 
             title="anti-HA affinity at late", yLabel="Lo Hi HA affinity (anti-H1)") + scale_y_continuous(breaks=seq(0,10,0.25), limits = c(0,1.1) )

FC_Affinity <- dcast(subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("Lo.Hi.HA.affinity"))
FC_Affinity$FC <- FC_Affinity$`late`/FC_Affinity$`baseline`
data=FC_Affinity; xData="Cohort"; yData="FC"; fillParam="Cohort"; title="anti-HA Affinity"; yLabel="fold-change (late / baseline)";
justforttest <- data[, c(xData,yData)]
fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
pValue <- fit$p
annotationInfo <- paste0("P = ", round(pValue, 2)); my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28)))
ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) + 
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5) + 
  geom_violin(draw_quantiles = 0.5) + ggbeeswarm::geom_quasirandom(size=7, pch=21, fill="black", color="white", alpha=0.35) + 
  ggtitle(title) + ylab(yLabel) +  theme_bw() +
  theme(axis.text = element_text(size=28,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=32,hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1,vjust=1)) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank())+ 
  annotation_custom(my_grob) + theme(legend.position = "none") + scale_y_continuous(breaks=seq(0,99,0.25), limits = c(0.75,1.75))
 

subsetData <- subset(mergedData, Cohort != "nonPD1" )
melted <- melt(subsetData, id.vars = c('dbCode', 'TimeCategory', 'Cohort'), measure.vars = c("X1.Kd.HA")); melted$value <- 1/melted$value
prePostTimeAveraged(melted, title = "HA EC50/Kd", xLabel = NULL, yLabel = "1 / Concentration (ug/mL)") + scale_y_continuous(breaks=seq(0,99,0.5), limits = c(0.5,3.5)) + 
  ylab(expression(paste("1 / H1 Concentration (",mu,"g/mL)"))) + ggtitle(expression(paste("E",C[50],"/",K[d])))
 # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/AntibodyAffinity_overTime.pdf", device="pdf", width=4.5, height=7)
summary(fit <- lm( data = mergedData, X1.Kd.HA ~ Cohort  + TimeCategory))
summary(fit <- aov(formula = X1.Kd.HA ~ Cohort  + TimeCategory + Cohort:TimeCategory, data = mergedData)); tukey_hsd(fit) %>% print(n=42)
# write.csv(melted, file="D:/Pembro-Fluvac/Analysis/sourceDataExports/2h.csv")



subsetData <- subset(mergedData, Cohort != "nonPD1" & TimeCategory == "baseline")
univScatter(data = subsetData, xData = "Cycle.of.Immunotherapy", yData = "Lo.Hi.HA.affinity", fillParam = "TimeCategory", 
            title = "Affinity over time", xLabel= "Cycle of Immunotherapy", yLabel = "Lo Hi Affinity")  + scale_fill_manual(values=c("#FFB18C"))

##' ## ----------- irAE analysis  --------------------
#'
#'

index <- demog[which(demog$irAE == "Y" | demog$irAE == "N"),"Subject"]
mergedData.irAE <- mergedData[which(mergedData$Subject %in% index),]  
subsetData <- subset(mergedData.irAE, TimeCategory == "baseline" & Cohort == "anti-PD-1")
subsetData$dateDiff <- round( as.numeric(difftime(subsetData$DateFirstIRAE,  subsetData$DateFirstFlowFile, units = "days")),digits = 0)
subsetData <- subsetData[-which(is.na(subsetData$dateDiff)), ]           # zero out the ones who did not have irAE
subsetData$Subject <- factor(subsetData$Subject, levels = subsetData$Subject[order(subsetData$dateDiff, decreasing = T)] )
ggplot(subsetData, aes(x=dateDiff, y=Subject)) + 
  geom_bar(stat="Identity", width=0.15, fill="#ff9a6a", size=0.01, color="black") + geom_point(aes(size=irAEgradeWorst)) + 
  xlab("Days until flu vaccine") + ylab("Participant") + labs(size = "irAE grade") + ggtitle("Time since first irAE") + scale_size(range=c(2,7),) + 
  theme_bw() + theme(axis.text.x = element_text(color = "black", size=24), axis.title = element_text(size=24), plot.title = element_text(size=32),
                     axis.text.y = element_blank(), legend.text = element_text(size=16), legend.title = element_text(size=16), 
                     ) + theme(panel.grid.minor = element_blank(), panel.grid.major=element_blank()) + geom_vline(xintercept = 0, size=0.5)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_dateDiff_irAEtoVaccine.pdf")
# write.csv(subsetData[,c("dbCode","dateDiff","irAEgradeWorst")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7a.csv")



subsetData <- mergedData[which(!is.na(mergedData$irAE) & mergedData$irAE != "" ), ]
twoSampleBar(data = subsetData, yData="FCtfh_oW", yLabel = "Fold-change at one week", title="Post-vax\nICOS+CD38+ cTfh", xData = "irAE",fillParam = "irAE", nonparam = F) + 
  theme(axis.title.x = element_text(size=28),  axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0)) +
  scale_fill_manual(values = c("grey90", "#ff9a6a")) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_FCtfh.pdf", width=3.5)
# write.csv(subsetData[,c("dbCode","TimeCategory","irAE","FCtfh_oW")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/4c.csv")

subsetData <- mergedData[which(!is.na(mergedData$irAE) & mergedData$irAE != "" ), ]
subsetData %>% group_by(irAE) %>% get_summary_stats(FCtfh_oW)
twoSampleBar(data = subsetData, yData="H1N1pdm09.HAI.titer", yLabel = "Baseline H1N1pdm09 titer", title="Baseline HAI", xData = "irAE",fillParam = "irAE")+ 
  theme(axis.title.x = element_text(size=28),  axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0), plot.title = element_text(size=28)) + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) 
twoSampleBar(data = subsetData, yData="FChai_late", yLabel = "FC HAI titer", title="Fold-change HAI", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a"))
twoSampleBar(data = subsetData, yData="IgDlo_CD71hi_ActBCells...FreqParent", yLabel = "Baseline ABC (% CD71)", title="ABC", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a"))
twoSampleBar(data = subsetData, yData="IgDlo_CD71hi_CD20loCD71hi...FreqParent", yLabel = "Baseline ASC (% CD71)", title="ASC", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a"))

twoSampleBar(data = subsetData, yData="cTfh_ICOShiCD38hi_cTfh..._medfi.ICOS.", yLabel = "medianFI ICOS", title="baseline\nICOS+CD38+ cTfh", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28),  axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_ICOS.pdf", width=4)
# write.csv(subsetData[,c("dbCode","TimeCategory","irAE","cTfh_ICOShiCD38hi_cTfh..._medfi.ICOS.")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7d.csv")

twoSampleBar(data = subsetData, yData="cTfh_ICOShiCD38hi_cTfh..._medfi.aIgG4..batchEffect", yLabel = "medianFI aIgG4", title="baseline\nICOS+CD38+ cTfh", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a"))  + theme(axis.title.x = element_text(size=28),  axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0)) # + geom_label_repel(aes(label = Subject),size=3)
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_aIgG4.pdf", width=4)
# write.csv(subsetData[,c("dbCode","TimeCategory","irAE","cTfh_ICOShiCD38hi_cTfh..._medfi.aIgG4..batchEffect")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7e.csv")




subsetData <- mergedData[which(!is.na(mergedData$irAE) & mergedData$irAE != "" ), ]
twoSampleBar(data = subsetData, yData="IgG1_Total.sialylated", yLabel = "Sialylation (% anti-H1 IgG1)", title="Sialylation", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))
twoSampleBar(data = subsetData, yData="IgG1sial_oW", yLabel = "Fold-change sialylation", title="Fold-change sialylation", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))
twoSampleBar(data = subsetData, yData="IgG1_Total.Galactosylation..G1.G2.", yLabel = "Galactosylation (% anti-H1 IgG1)", title="Galactosylation", xData = "irAE",fillParam = "irAE")+ 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))
twoSampleBar(data = subsetData, yData="Lo.Hi.HA.affinity", yLabel = "Lo Hi Affinity", title="Affinity", xData = "irAE",fillParam = "irAE") + 
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))


subsetData <- mergedData[which(!is.na(mergedData$irAE) & mergedData$irAE != "" & mergedData$TimeCategory == "baseline" & mergedData$Cohort == "anti-PD-1"), ]
twoSampleBar(data = subsetData, yData="cTfh_ICOShiCD38hi_..FreqParent", yLabel = expression(paste("ICO",S^"+", "CD3",8^"+"," (% cTfh)")), 
             title="Baseline\nICOS+CD38+ cTfh", xData = "irAE", fillParam = "irAE", nonparam = T) + 
  theme(axis.title.x = element_text(size=28),  axis.text.x = element_text(hjust = 0.5, vjust=0.5, angle=0)) +
  scale_fill_manual(values = c("grey90", "#ff9a6a")) 
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_cTfh_Baseline.pdf", width=3.5)
# write.csv(subsetData[,c("dbCode","TimeCategory","irAE","cTfh_ICOShiCD38hi_..FreqParent")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7j.csv")



subsetData <- subset(mergedData, Cohort != "nonPD1" )
twoSampleBar(data = subsetData, yData="ANA.titer", yLabel = "ANA titer", title="Anti-nuclear\nantibodies", xData = "Cohort",fillParam = "Cohort", nonparam = T)
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/ANA_baseline_fullGroup.pdf", width=4)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","ANA.titer")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7k-1.csv")

twoSampleBar(data = subsetData, yData="Anti.dsDNA.titer", yLabel = "Anti-dsDNA titer", title="Anti-dsDNA\nantibodies", xData = "Cohort",fillParam = "Cohort", nonparam = T)
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/Anti-dsDNA_baseline_fullGroup.pdf", width=4)
# write.csv(subsetData[,c("dbCode","Cohort","TimeCategory","Anti.dsDNA.titer")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7k-2.csv")

subsetData <- mergedData[which(!is.na(mergedData$irAE) & mergedData$irAE != "" ), ]
twoSampleBar(data = subsetData, yData="ANA.titer", yLabel = "ANA titer", title="Anti-nuclear\nantibodies", xData = "irAE",fillParam = "irAE", nonparam = T, position = "right") +
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_ANA_baseline.pdf", width=4)
# write.csv(subsetData[,c("dbCode","TimeCategory","irAE","ANA.titer")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7l-1.csv")

twoSampleBar(data = subsetData, yData="Anti.dsDNA.titer", yLabel = "Anti-dsDNA titer", title="Anti-dsDNA\nantibodies", xData = "irAE",fillParam = "irAE", nonparam=T, position = "right") +
  scale_fill_manual(values = c("grey90", "#ff9a6a")) + theme(axis.title.x = element_text(size=28))
  # ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/irAE_anti-dsDNA_baseline.pdf", width=4)
# write.csv(subsetData[,c("dbCode","TimeCategory","irAE","Anti.dsDNA.titer")], file="D:/Pembro-Fluvac/Analysis/sourceDataExports/e7l-2.csv")
