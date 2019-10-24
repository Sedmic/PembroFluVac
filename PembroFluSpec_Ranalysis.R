library("grid")
library("ggplot2")
library("gplots")
library("viridis")
library("reshape2")
library("ggpubr")
library("ggcorrplot")
library("gridExtra")
library("RColorBrewer")
library("scales")
library("MASS")



## ------------------------------------------- PLOTTING FUNCTIONS ------------------------------------------------------------------------------


twoSampleBox <- function (data, xData, yData, fillParam, title, yLabel)
{
  fit <- lm(data[,yData] ~ data[,xData])
  pValue <- summary(fit)$coefficients[2,"Pr(>|t|)"];  CI <- confint(fit)[2,]; CI <- round(CI,2)
  if (pValue < 0.01)
  {
    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1), "\n", "Mean 95%CI: [",CI[1], ", ",CI[2], "]")
  }
  if (pValue >= 0.01)
  {
    annotationInfo <- paste0("P = ", round(pValue, 2),"\n", "Mean 95%CI: [",CI[1], ", ",CI[2], "]")
  }
  my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
  return (
    ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) +  
      geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
      ggtitle(title) + ylab(yLabel) +  
      theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")    
  )
}


univScatter <- function(data, xData, yData, fillParam, title, xLabel, yLabel)
{
  pearson <- round(cor(data[,xData], data[,yData], method = "pearson", use = "complete.obs"), 2)
  pValue <- cor.test(data[,xData], data[,yData], method="pearson")
  if (pValue$p.value < 0.01)
  {
    annotationInfo <- paste0("Pearson r = ", pearson,"\n","P = ", formatC(pValue$p.value, format="e", digits=1))
  }
  if (pValue$p.value >= 0.01)
  {
    annotationInfo <- paste0("Pearson r = ", pearson,"\n","P = ", round(pValue$p.value,2))
  }
  
  my_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
  return (
    ggplot(data ) + 
      geom_smooth(aes_string(x=xData, y=yData), method='lm') + aes(color="one",  fill = "two") +
      geom_point(aes_string(x=xData, y=yData), size=5, fill=quote('black')) + theme_bw() + 
      scale_color_manual(name="Val", values=c(one="black")) + 
      scale_fill_manual(name="fills", values=c(two="grey90")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")   
  )
}





setwd("D:/Pembro-Fluvac/Analysis")
mergedData <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_freqParent_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)

serology <- read.csv(file= "D:/Pembro-Fluvac/Analysis/mergedData/SerologyMeasurements.csv", stringsAsFactors = F, header = T)
temp1 <- merge(x = mergedData, y=serology, all = T, suffixes = c(".Tflow",".Serology"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year'))

Bflow <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_freqParent_allyrs.csv", stringsAsFactors = F, header = T, row.names = 1)
temp2 <- merge(x = temp1, y=Bflow, all = T, suffixes = c(".Tflow",".Bflow"), by= c('Label', 'TimeCategory','Cohort','Subject','TimePoint','Year'))

mergedData <- temp2

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
a <- univScatter(subsubsetData, xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam = NULL, 
                      title = "Yr3: aPD1 subjects at one week", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")
a + scale_x_continuous(breaks=seq(0,99,0.25)) + scale_y_continuous(breaks=seq(0,99,1))

subsetData <- subset(mergedData, Cohort == "Healthy" & Year == "3" )
subsubsetData <- subset(subsetData, TimeCategory == "oneWeek")
a <- univScatter(subsubsetData, xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam = NULL, 
                      title = "Yr3: Healthy subjects at one week", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")
a + scale_x_continuous(breaks=seq(0,99,0.25)) + scale_y_continuous(breaks=seq(0,99,1))



subsetData <- subset(mergedData, Cohort == "aPD1" & TimeCategory == "oneWeek" )
univScatter(subsetData, xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam = NULL, 
            title = "All yrs: aPD1 subjects at one week", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")


subsetData <- subset(mergedData, Cohort == "Healthy" & TimeCategory == "oneWeek" )
a <- univScatter(subsetData, xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam = NULL, 
            title = "All yrs: Healthy subjects at one week", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")
a + scale_x_continuous(breaks=seq(0,99,0.25)) + scale_y_continuous(breaks=seq(0,99,1))



## ------------------------------------------------- Scatterplots of PB vs cTfh  ----------------------------------------------------------


subsetData <- subset(mergedData, Cohort == "aPD1" & TimeCategory == "oneWeek" )
univScatter(subsetData, xData = "HAI.titer..H1N1pdm09.", yData="Plasma.CXCL13..pg.mL.", fillParam = NULL, 
                 title = "All yrs: aPD1 subjects at one week", xLabel= "HAI titer - H1N1pdm09", yLabel = "Plasma CXCL13 (pg/mL)")
 



























## ------------------------------------------- Deprecated  ------------------------------------------------------------------------------




