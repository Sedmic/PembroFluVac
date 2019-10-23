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
library("Rlabkey")
labkey.setDefaults(apiKey="apikey|2280ef3ab44cc0f3e090af59e71c2df2")

labkey.data <- labkey.selectRows(
  baseUrl="http://localhost:8080/labkey", 
  folderPath="/PembroFluVac", 
  schemaName="study", 
  queryName="allStudyVisits", 
  viewName="allDataMerged", 
  colFilter=NULL, 
  containerFilter=NULL
)


# data to exclude: 
# Cohort:  nonPD1
# samples:  Tpanel 1718:  PBMC_19616-007_d7_002.fcs

finalDataset <- labkey.data
# finalDataset <- subset(labkey.data, !(`Cohort`=="nonPD1"))
# finalDataset <- finalDataset[-which(finalDataset$`Participant ID` == "19616-007" & finalDataset$Year==2 & finalDataset$`Time Category`=="oneWeek"),]
finalDataset$`Time Category` <- factor(finalDataset$`Time Category`, levels=c("baseline","oneWeek","late"))
allD7 <- subset(finalDataset, `Time Category` == "oneWeek")
Yr2data <- subset(finalDataset, `Year`==2)
Yr2D0 <- subset(Yr2data, `Time Category`=="baseline")
Yr2D7 <- subset(Yr2data,  `Time Category` == "oneWeek")
Yr3data <- subset(finalDataset, `Year`==3)
Yr3D7 <- subset(Yr3data,  `Time Category` == "oneWeek")

Yr2D0D7 <- subset(Yr2data, `Time Category` == "baseline" | `Time Category` == "oneWeek"); Yr2D0D7 <- Yr2D0D7[order(Yr2D0D7$`Participant ID`),]


###################################################################################################### --------------- Yr3D7 ICOS+CD38+ by Cohort ----------------
fit <- lm(Yr3D7$`ICOS+CD38+ c Tfh Freq`~Yr3D7$Cohort)
summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))

ggplot(data=Yr3D7, aes(x=`Cohort`, y=`ICOS+CD38+ c Tfh Freq`, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) +  
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
  ggtitle("Year 3 day 7: cTfh responses") + ylab("ICOS+CD38+ cTfh freq at d7") + scale_y_continuous(breaks=seq(1:12), limits=c(0.5,11)) +  
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  annotation_custom(my_grob1) + theme(legend.position = "none")
# ggsave(filename = "Images/cTfhresponses_Yr3D7.pdf")

###################################################################################################### --------------- Yr3D7 plasmablast Cohort ----------------
fit <- lm(Yr3D7$`CD19+/CD27+CD38+ Freq`~Yr3D7$Cohort)
summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=Yr3D7, aes(x=`Cohort`, y=`CD19+/CD27+CD38+ Freq`, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
  ggtitle("Year 3 day 7: Plasmablast responses") + ylab("Plasmablast freq at d7") + scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  annotation_custom(my_grob1) + theme(legend.position = "none")
# ggsave(filename = "Images/Plasmablastresponses_Yr3D7.pdf")

###################################################################################################### --------------- Alld7 ICOS+CD38+ by Cohort ----------------
temp <- allD7[-which(allD7$`Participant ID` == "19616-007" & allD7$Year==2),]
fit <- lm(temp$`ICOS+CD38+ c Tfh Freq`~temp$Cohort)
summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=temp, aes(x=`Cohort`, y=`ICOS+CD38+ c Tfh Freq`, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
  ggtitle("All years day 7: cTfh responses") + ylab("ICOS+CD38+ cTfh freq at d7") + scale_y_continuous(breaks=seq(1:16), limits=c(1,16)) +  
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  annotation_custom(my_grob1) + theme(legend.position = "none")
# ggsave(filename = "Images/cTfhresponses_allD7.pdf")

###################################################################################################### --------------- Alld7 plasmablast by Cohort ----------------
temp <- allD7[-which(allD7$`Participant ID` == "19616-007" & temp$Year==2),]
fit <- lm(temp$`CD19+/CD27+CD38+ Freq`~temp$Cohort)
summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=temp, aes(x=`Cohort`, y=`CD19+/CD27+CD38+ Freq`, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
  ggtitle("All years day 7: Plasmablast responses") + ylab("Plasmablast freq at d7") + scale_y_continuous(breaks=seq(1:5), limits=c(0,5)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  annotation_custom(my_grob1) + theme(legend.position = "none")
# ggsave(filename = "Images/Plasmablastresponses_allD7.pdf")

###################################################################################################### --------------- Yr2D0 CXCL13 by Cohort ----------------
mean(Yr2D0$`Plasma CXCL13 (pg/m L)`, na.rm = T)+3*sd(Yr2D0$`Plasma CXCL13 (pg/m L)`, na.rm = T)               # one outlier, FS1718-120 d7 is >3sd away so will exclude
temp <- Yr2D0[-which(Yr2D0$`Plasma CXCL13 (pg/m L)`== max(Yr2D0$`Plasma CXCL13 (pg/m L)`, na.rm = T)),]
fit <- lm(temp$`Plasma CXCL13 (pg/m L)`~temp$Cohort)
summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=temp, aes(x=`Cohort`, y=`Plasma CXCL13 (pg/m L)`, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
  ggtitle("Plasma CXCL13 at day 0") + ylab("CXCL13 (pg/mL)") + # scale_y_continuous(breaks=seq(1:5), limits=c(0,5)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  annotation_custom(my_grob1) + theme(legend.position = "none")
# ggsave(filename = "Images/CXCL13plasma_Yr2D0.pdf")

###################################################################################################### --------------- Alld7 CXCL13 by Cohort ----------------
mean(allD7$`Plasma CXCL13 (pg/m L)`, na.rm = T)+3*sd(allD7$`Plasma CXCL13 (pg/m L)`, na.rm = T)               # one outlier, FS1718-120 d7 is >3sd away so will exclude
temp <- allD7[-which(allD7$`Plasma CXCL13 (pg/m L)`== max(allD7$`Plasma CXCL13 (pg/m L)`, na.rm = T)),]
fit <- lm(temp$`Plasma CXCL13 (pg/m L)`~temp$Cohort)
summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=temp, aes(x=`Cohort`, y=`Plasma CXCL13 (pg/m L)`, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
  ggtitle("Plasma CXCL13 at day 7") + ylab("CXCL13 (pg/mL)") + # scale_y_continuous(breaks=seq(1:5), limits=c(0,5)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  annotation_custom(my_grob1) + theme(legend.position = "none")
# ggsave(filename = "Images/CXCL13plasma_allD7.pdf")


###################################################################################################### --------------- FC ICOS+CD38+ by Cohort ----------------

FC_Tfhresponse <- dcast(finalDataset, `Participant ID`+`Cohort`~`Time Category`, value.var = c("ICOS+CD38+ c Tfh Freq")); FC_Tfhresponse$FC <- FC_Tfhresponse$oneWeek/FC_Tfhresponse$baseline
fit <- lm(FC_Tfhresponse$FC~FC_Tfhresponse$Cohort)
summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=FC_Tfhresponse, aes(x=`Cohort`, y=FC, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
  ggtitle("All years: ICOS+CD38+ cTfh responses") + ylab("ICOS+CD38+ cTfh fold-change") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  annotation_custom(my_grob1) + theme(legend.position = "none")
# ggsave(filename = "Images/cTfhResponses_FCallyears.pdf")

###################################################################################################### --------------- FC Plasmablast by Cohort ----------------


FC_PBresponse <- dcast(finalDataset, `Participant ID`+`Cohort`~`Time Category`, value.var = c("CD19+/CD27+CD38+ Freq")); FC_PBresponse$FC <- FC_PBresponse$oneWeek/FC_PBresponse$baseline
fit <- lm(FC_PBresponse$FC~FC_PBresponse$Cohort)
summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=FC_PBresponse, aes(x=`Cohort`, y=FC, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
  ggtitle("All years: ICOS+CD38+ cTfh responses") + ylab("ICOS+CD38+ cTfh fold-change") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  annotation_custom(my_grob1) + theme(legend.position = "none")
# ggsave(filename = "Images/PlasmablastResponses_FCallyears.pdf")



###################################################################################################### --------------- D0-D7 plots by Cohort ----------------

temp <- Yr3data[-which(Yr3data$`Time Category`=="late"),]
ggplot(temp, aes(x=`Time Category`, y=`ICOS+CD38+ c Tfh Freq`, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) +  theme_bw() +
  geom_boxplot(size=1) +  # geom_line(aes(group=`Participant ID`),size=1) + 
  ggtitle("All years: ICOS+CD38+ cTfh responses") + ylab("ICOS+CD38+ cTfh freq") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) 
# ggsave(filename = "Images/cTfhResponses_BeforeAfterallyears.pdf",width=7, height=7)

ggplot(temp, aes(x=`Time Category`, y=`CD19+/CD27+CD38+ Freq`, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) +  theme_bw() +
  geom_boxplot(size=1) +  # geom_line(aes(group=`Participant ID`),size=1) + 
  ggtitle("All years: Plasmablast responses") + ylab("Plasmablast freq") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) 
# ggsave(filename = "Images/PlasmablastResponses_BeforeAfterallyears.pdf",width=7, height=7)



###################################################################################################### --------------- HAI titers by Cohort ----------------

temp <- finalDataset[-which(finalDataset$`Time Category`=="oneWeek"),]
ggplot(temp, aes(x=`Time Category`, y=`HAI titer (H1N1pdm09)`, fill=`Cohort`)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) +  theme_bw() +
  geom_boxplot() + #geom_line(aes(group=`Participant ID`),size=1) + 
  ggtitle("All years: HAI titers") + ylab("H1N1pdm09 titer") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) 








setwd("D:/Pembro-Fluvac/Analysis")
UPennYr3 <- read.csv(file = "D:/Pembro-Fluvac/18-19season/Flow cytometry/Tcell/Analysis/mergeData_freqParent.csv", stringsAsFactors = F, header = T, row.names = 1)
UPennYr3 <- UPennYr3[which(UPennYr3$Cohort != "nonPD1"),]

fit <- lm(UPennYr3$`cTfh_ICOShiCD38hi_..FreqParent` ~ UPennYr3$Cohort);  summary(fit)
pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("P = ", pValue,"\n","95%CI: [",CI[1],", ",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
a <- twoSampleBox(data=UPennYr3[which(UPennYr3$Cohort != "nonPD1" & UPennYr3$TimeCategory == "oneWeek"),], xData="Cohort", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam="Cohort", 
             title="ICOS+CD38+ cTfh at d7", yLabel="ICOS+CD38+ cTfh frequency at d7", myGrob = my_grob1)
a + scale_y_continuous(breaks=seq(1:12), limits=c(0.5,11)) 


subsetData <- UPennYr3[which(UPennYr3$Cohort != "nonPD1" & UPennYr3$TimeCategory == "oneWeek"), ]
fit <- lm(data[, xData] ~ data[, yData] ) 
pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);  CI <- confint(fit)[2,]; CI <- round(CI,2)
annotationInfo <- paste0("P = ", pValue,"\n","95%CI: [",CI[1], ", ",CI[2], "]")
my_grob = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
twoSampleBox(data= subsetData, xData="Cohort", yData="CD19_CD27.CD38....FreqParent", fillParam="Cohort", title="Plasmablasts at d7", yLabel="Plasmablast freq at d7")

subsetData <- subset(UPennYr3, Cohort == "aPD1" )
subsubsetData <- subset(subsetData, TimeCategory == "oneWeek")
fit <- lm(subsubsetData$`CD19_CD27.CD38....FreqParent` ~ subsubsetData$`cTfh_ICOShiCD38hi_..FreqParent` );  summary(fit)
a <- twoSampleScatter(subsubsetData, xData = "CD19_CD27.CD38....FreqParent", yData="cTfh_ICOShiCD38hi_..FreqParent", fillParam = NULL, 
                      title = "aPD1 subjects at one week", xLabel= "Plasmablast frequency", yLabel = "ICOS+CD38+ cTfh frequency")
a + scale_x_continuous(breaks=seq(0,99,0.25)) + scale_y_continuous(breaks=seq(0,99,1))


## ------------------------------------------- PLOTTING FUNCTIONS ------------------------------------------------------------------------------


twoSampleBox <- function (data, xData, yData, fillParam, title, yLabel, my_grob)
{
  return (
    ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) +  
      geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
      ggtitle(title) + ylab(yLabel) +  
      theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")    
  )
}


twoSampleScatter <- function(data, xData, yData, fillParam, title, xLabel, yLabel)
{
  pearson <- round(cor(data[,xData], data[,yData], method = "pearson"), 2)
  pValue <- cor.test(data[,xData], data[,yData], method="pearson")
  annotationInfo <- paste0("Pearson r = ", pearson,"\n","P = ", round(pValue$p.value,2)); my_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
  return (
    ggplot(data ) + 
      geom_point(aes_string(x=xData, y=yData), size=5, fill=quote('black')) + theme_bw() + 
      geom_smooth(aes_string(x=xData, y=yData), method='lm') + aes(color="one",  fill = "two") +
      scale_color_manual(name="Val", values=c(one="black")) + 
      scale_fill_manual(name="fills", values=c(two="grey90")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")   
  )
}













