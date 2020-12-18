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
# library("Rlabkey")
# labkey.setDefaults(apiKey="apikey|2280ef3ab44cc0f3e090af59e71c2df2")
# labkey.data <- labkey.selectRows(  baseUrl="http://localhost:8080/labkey", folderPath="/PembroFluVac", schemaName="study", queryName="Rockefeller_freqParent", viewName="demo_freqParent", colFilter=NULL, containerFilter=NULL)
# dataset <- labkey.data
# dataset$`Visit` <- factor(dataset$`Visit`, levels=c("1","2","3"))

# FC_Tfhresponse <- dcast(dataset, `Participant ID`+`Trx group`~`Visit`, value.var = c("ICOS+CD38+ c Tfh freq")); 
# FC_Tfhresponse$FC <- FC_Tfhresponse$`2`/FC_Tfhresponse$`1`
#  exclude subject 87 because day 7 timepoint had unreasonably few CD4 recorded to be confident of the point estimate
# FC_Tfhresponse <- FC_Tfhresponse[-which(FC_Tfhresponse$`Participant ID` == 86),]
# FC_Tfhresponse$`Trx group`<- factor(FC_Tfhresponse$`Trx group`, levels=c("Chemotherapy","PD-L1 only", "TKI only", "PD1 only", "PD1 + TKI"))
# fit <- lm(FC_Tfhresponse$FC~FC_Tfhresponse$`Trx group`); summary(fit)
# summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
# annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
# ggplot(data=FC_Tfhresponse, aes(x=`Trx group`, y=FC, fill=`Trx group`)) + theme_bw() + 
#   geom_boxplot(outlier.shape=NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black") + 
#   ggtitle("Rockefeller: ICOS+CD38+ cTfh responses") + ylab("ICOS+CD38+ cTfh fold-change") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
#   theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
#   # annotation_custom(my_grob1) + 
#   theme(legend.position = "none")
  

# FC_Tfhresponse <- dcast(dataset, `Participant ID`+`Any A PD1`~`Visit`, value.var = c("ICOS+CD38+ c Tfh freq")); 
# FC_Tfhresponse$FC <- FC_Tfhresponse$`2`/FC_Tfhresponse$`1`
# #  exclude subject 87 because day 7 timepoint had unreasonably few CD4 recorded to be confident of the point estimate
# FC_Tfhresponse <- FC_Tfhresponse[-which(FC_Tfhresponse$`Participant ID` == 86),]
# FC_Tfhresponse$`Any A PD1`<- factor(FC_Tfhresponse$`Any A PD1`, levels=c("F","T"))
# fit <- lm(FC_Tfhresponse$FC~FC_Tfhresponse$`Any A PD1`); summary(fit)
# # summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
# # annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
# ggplot(data=FC_Tfhresponse, aes(x=`Any A PD1`, y=FC, fill=`Any A PD1`)) + theme_bw() + 
#   geom_boxplot(outlier.shape=NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black") + 
#   ggtitle("Rockefeller: ICOS+CD38+ cTfh responses") + ylab("ICOS+CD38+ cTfh fold-change") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
#   theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
#   # annotation_custom(my_grob1) + 
#   theme(legend.position = "none")


# FC_Tfhresponse <- dcast(dataset, `Participant ID`+`Treatment`~`Visit`, value.var = c("ICOS+CD38+ c Tfh freq")); 
# FC_Tfhresponse$FC <- FC_Tfhresponse$`2`/FC_Tfhresponse$`1`
# #  exclude subject 87 because day 7 timepoint had unreasonably few CD4 recorded to be confident of the point estimate
# FC_Tfhresponse <- FC_Tfhresponse[-which(FC_Tfhresponse$`Participant ID` == 86),]
# FC_Tfhresponse$`Treatment`<- factor(FC_Tfhresponse$`Treatment`, levels=c("Carbo Gem","axitinib","cabozantinib","pazopanib",
#                                                                         "Atezolizumab","Pembrolizumab","Nivolumab",
#                                                                         "Pembrolizumab/lenvatinib","Nivolumab/cabozantinib"))
# fit <- lm(FC_Tfhresponse$FC~FC_Tfhresponse$`Treatment`); summary(fit)
# # summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
# # annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
# ggplot(data=FC_Tfhresponse, aes(x=`Treatment`, y=FC, fill=`Treatment`)) + theme_bw() + 
#   geom_boxplot(outlier.shape=NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black") + 
#   ggtitle("Rockefeller: ICOS+CD38+ cTfh responses") + ylab("ICOS+CD38+ cTfh fold-change") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
#   theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), 
#         plot.title = element_text(size=28,hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1)) + 
#



#' ## ----------- Dataset assembly  --------------------
#'
#'
#'
#'
setwd("D:/Pembro-Fluvac/Analysis")

mergedDataRock <- read.csv(file="D:/Pembro-Fluvac/Analysis/mergedData/Rockefeller_SubjectTable.csv")

serology <- read.csv(file= "D:/Pembro-Fluvac/Analysis/mergedData/Rockefeller_Serology.csv", header = T)
print("Data by year and cohort for serology")
temp1 <- merge(x = mergedDataRock, y=serology, all = T, suffixes = c(".Tflow",".Serology"), by= c('Subject'))
table(serology$Cohort, serology$TimeCategory)

TflowRock <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/Rockefeller_Tflow_freqParent.csv", header = T)
print("Data by year and cohort for T cell flow cytometry for Rockefeller Cohort")
# temp <- sapply(TflowRock, sd)
# fit <- lm(data=TflowRock[ , - which( colnames(TflowRock) == "Label" | colnames(TflowRock) == "fcsfile"  | colnames(TflowRock) == 'CD19..CD27.CD38..TCF1hi.freq.1')], formula = Batch  ~ .)
# a <- as.data.frame(summary(fit)$coefficients)
# a <- rownames(a[which(a$`Pr(>|t|)` < 0.05), ] )
# if (length (a) > 0) { colnames(TflowRock) [ grep( paste(a, collapse="|"), colnames(TflowRock)) ] <- paste0( colnames(TflowRock)[ grep( paste(a, collapse="|"), colnames(TflowRock)) ], ".batchEffect" ) }

temp2 <- merge(x = temp1, y=TflowRock, all=T, suffixes = c(".temp1",".TflowRock"), by=c('Label'))
mergedDataRock <- temp2   # annotation_custom(my_grob1) + 
#   theme(legend.position = "none") 
mergedDataRock <- mergedDataRock[-which(mergedDataRock$Subject == "13-226-91"),]    # ambiguity about what time points are being studied
mergedDataRock <- mergedDataRock[-which(mergedDataRock$Subject == "13-226-67"),]    # ambiguity about what time points are being studied
mergedDataRock <- mergedDataRock[-which(mergedDataRock$Label == "13-226-51E"),]     # two draws at oneWeek, this one is at d10

table(mergedDataRock$Cohort, mergedDataRock$TimeCategory)

mergedDataRock$logHAI <- log(mergedDataRock$H1N1pdm09.HAI.titer)

#' ## ----------- calculate fold-changes for HAI, Tfh, and PB --------------------
subsetData <- subset(mergedDataRock, TimeCategory != "oneWeek" &  Cohort!="")
FC_response <- dcast(subsetData, `Subject` + `Cohort` ~`TimeCategory`, value.var = c("logHAI"))
FC_response$FChai_late <- FC_response$`late`/FC_response$`baseline`; FC_response$baseline <- FC_response$late <- NULL

subsetData <- subset(mergedDataRock, TimeCategory != "late" &  Cohort!="")
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("ICOS.CD38..cTfh.freq")) 
FC_response2$FCtfh_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject",); FC_response$baseline <- FC_response$oneWeek <- NULL
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("CD19..CD27.CD38..freq")) 
FC_response2$FCPB_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject"); FC_response$baseline <- FC_response$oneWeek <- NULL
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("Plasma.CXCL13..pg.mL.")) 
FC_response2$FCCXCL13_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject"); FC_response$baseline <- FC_response$oneWeek <- NULL
FC_response2 <- dcast( subsetData, `Subject`+`Cohort`~`TimeCategory`, value.var = c("IgG1_Total.sialylated")) 
FC_response2$IgG1sial_oW <- FC_response2$`oneWeek`/FC_response2$`baseline`; FC_response2$Cohort <- NULL
FC_response <- merge(x=FC_response, y=FC_response2, by = "Subject"); FC_response$baseline <- FC_response$oneWeek <- NULL

mergedDataRock <- merge(x = mergedDataRock, y= FC_response, all=T, by = c('Subject', 'Cohort'))

# write.csv(mergedDataRock, file = "D:/Pembro-Fluvac/Analysis/mergedDataRock/allmergedDataRock.csv")


# mergedDataRock <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedDataRock/allmergedDataRock.csv", stringsAsFactors = F, header = T, row.names = 1)

mergedDataRock$dummy <- "dummy"
mergedDataRockComplete <- temp <- mergedDataRock

temp$Cohort <- factor(temp$Cohort, levels = c("Chemo","TKI","aPDL1","aPD1","aPD1_TKI","IpiNivo"))

data=subset(temp, Cohort != "nonPD1" & TimeCategory == "oneWeek"); xData = "Cohort"; yData = "FCtfh_oW"; fillParam = 'Cohort'; title = "ICOS+CD38+ cTfh responses"; yLabel = "Fold-change at one week";
overTime <- aggregate(x = data[,yData], by= list(Cohort = data$Cohort), FUN=mean, na.rm = T)  # make a table of all the summary statistics for the bars
names(overTime)[which(names(overTime) == 'x')] <- yData

ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + #scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) + 
  geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity") + 
  geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.15)) + 
  ggtitle(title) + ylab(yLabel) +  theme_bw() +
  theme(axis.text = element_text(size=28,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=32,hjust = 0.5)) 


mergedDataRock <- mergedDataRock[ grep(paste0(c("aPD1$","Chemo","aPDL1"), collapse = '|'), mergedDataRock$Cohort), ]   # simplify down to aPD1 vs non-aPD1
mergedDataRock$Cohort <- ifelse(mergedDataRock$Cohort == 'aPD1', 'aPD1', 'non-aPD1')
mergedDataRock$Cohort <- factor(mergedDataRock$Cohort, levels = c("non-aPD1","aPD1"))
mergedDataRock.nolate <- mergedDataRock[-which(mergedDataRock$TimeCategory == 'late'),]
mergedDataRock$TimeCategory <- factor(mergedDataRock$TimeCategory, levels = c("baseline", "oneWeek","late"))
mergedDataRock.nolate$TimeCategory <- factor(mergedDataRock.nolate$TimeCategory, levels = c("baseline","oneWeek"))


#' ## ----------- Tfh analysis --------------------
#'
#'
subsetData <- subset(mergedDataRock, TimeCategory == 'oneWeek')
twoSampleBar(data=subsetData, xData="Cohort", yData="FCtfh_oW", fillParam="Cohort", title="ICOS+CD38+ cTfh responses", yLabel="Fold-change at one week", position="left")
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_PerSubject_ROCK.pdf", width = 8)

prePostTime(data = mergedDataRock.nolate, xData = "TimeCategory", yData = "ICOS.CD38..cTfh.freq", fillParam = "Cohort", title = "cTfh response by subject", xLabel = "TimeCategory",
            yLabel = "ICOS+CD38+ (% cTfh)", groupby = "Subject")  + scale_y_continuous(breaks=seq(0,150,1), limits=c(0,18))     

melted <- melt(mergedDataRock.nolate, id.vars = c('Subject', 'TimeCategory', 'Cohort'), measure.vars = c("ICOS.CD38..cTfh.freq"))
melted <- melted[which(melted$TimeCategory != 'late'),]
prePostTimeAveraged(data = melted, title = "cTfh", xLabel = NULL, yLabel = "ICOS+CD38+ (% cTfh)")# + scale_y_continuous(breaks=seq(1,10,0.5), limits=c(3.5, 7))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/cTfhResp_overTime_freqcTfh.pdf", width=4)
summary( fit <- lm(formula = value ~ Cohort * TimeCategory, data = melted) );  tukey_hsd(fit)       # testing using unpaired two-way ANOVA with Tukey's posthoc
melted %>% group_by(Cohort) %>% kruskal_test(formula = value ~ TimeCategory)                        # testing per cohort using Kruskall-Wallis without adjustment




















