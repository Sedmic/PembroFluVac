library(ggplot2)
library(reshape2)
library(stringr)


#  -------------------------------------  MERGE AND CONCATENATE T CELL DATA --------------------------------------------------
a <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Tcell/Analysis/FreqParent")
files <-   list.files("./", pattern=').csv$', full=F)   # this list is all of the duplicates that Windows Explorer appends with (x) postfix to filename
mergeDataYr3 <- read.csv(file = "FreqParent.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr3$Batch = "0"

for (i in 1:length(files))
{
  temp <- read.csv(files[i], stringsAsFactors = F, header=T, row.names=1)
  temp$Batch = i
  mergeDataYr3 <- rbind(mergeDataYr3, temp)
}

mergeDataYr2 <- read.csv(file = "D:/Pembro-Fluvac/17-18season/TcellPanel/Analysis/FreqParent.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr2$Batch = i+1

setwd("D:/Pembro-Fluvac/16-17season/Analysis/freqParent/")
files <-   list.files("./", pattern=').csv$', full=F) # this list is all of the duplicates that Windows Explorer appends with (x) postfix to filename
mergeDataYr1 <- read.csv(file = "byParent-2019edition.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr1$Batch = i+2
for (j in 1:length(files))
{
  temp <- read.csv(files[j], stringsAsFactors = F, header=T, row.names=1)
  temp$Batch = i + 2 + j
  mergeDataYr1 <- rbind(mergeDataYr1, temp)
}

# mergeData <- mergeDataYr1
cleanColumnNames <- function (mergeData)
{
  mergeData <- mergeData[-grep(pattern="Mean", rownames(mergeData)),]
  mergeData <- mergeData[-grep(pattern="SD", rownames(mergeData)),]
  mergeData <- mergeData[,-which(colnames(mergeData) == "X.1")]
  
  temp <- rownames(mergeData); temp <- str_replace(temp, "PBMCs_", ""); temp <- str_replace(temp, "PBMC_", ""); 
  temp <- str_replace(temp, "day", "d"); temp <- substr(temp, 1, nchar(temp)-8)   # clean up rownames a bit, lop off the last 8 characters _0xx.fcs
  rownames(mergeData) <- temp
  # add columns for Cohort, Subject, TimePoint (numeric), TimeCategory (categorical), and Year
  mergeData$Label <- rownames(mergeData)
  mergeData$Subject <- word(mergeData$Label,1,sep = "\\_")  # take first chunk prior to _ character
  mergeData$Cohort <- 0;  
  mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "19")] <- "aPD1";  
  mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "FS")] <- "Healthy"; 
  mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "20")] <- "nonPD1"
  mergeData$TimePoint <- as.numeric(word(mergeData$Label, 2, sep= "d"))  # parse out the number of days since vaccination
  mergeData$TimeCategory <- 0;  
  mergeData$TimeCategory[which(mergeData$TimePoint == 0)] <- "baseline"  
  mergeData$TimeCategory[which(mergeData$TimePoint > 0 & mergeData$TimePoint <= 14)] <- "oneWeek"; 
  mergeData$TimeCategory[which(mergeData$TimePoint > 14)] <- "late"

  temp <- colnames(mergeData)
  temp <- str_replace(temp,"Freq..of.Parent....", "FreqParent")
  temp <- str_replace(temp,"Time.Lymphocytes.Singlets.Live.CD3..CD4..Nonnaive.CXCR5..", "cTfh_")
  temp <- str_replace(temp,"Time.Lymphocytes.Singlets.Live.CD19..", "CD19_")
  temp <- str_replace(temp,"Time.Lymphocytes.Singlets.Live.", "Live_")
  temp[grep(pattern = "FreqParent.1", temp)] <- str_replace( temp[grep(pattern = "FreqParent.1", temp)],"cTfh_ICOS.CD38..","cTfh_ICOSloCD38lo_") 
  temp <- str_replace( temp,"cTfh_ICOSloCD38lo.","cTfh_ICOSloCD38lo_") 
  temp <- str_replace(temp,"cTfh_ICOS.CD38..","cTfh_ICOShiCD38hi_")
  temp[grep(pattern = "Live_CD3..CD4....FreqParent.1", temp)] <- str_replace( temp[grep(pattern = "Live_CD3..CD4....FreqParent.1", temp)], "Live_CD3..CD4.", "Live_CD3hi_CD4lo_" )
  temp <- str_replace(temp,"FreqParent.1","FreqParent")
  temp <- str_replace(temp,"cTfh\\.","")
  colnames(mergeData) <- temp
  return(mergeData)  
}

Yr3data <- cleanColumnNames(mergeDataYr3);  Yr3data$Year <- 3
Yr2data <- cleanColumnNames(mergeDataYr2);  Yr2data$Year <- 2
Yr1data <- cleanColumnNames(mergeDataYr1);  Yr1data$Year <- 1

combinedYrs <- merge(x=Yr3data, y=Yr2data, all=T)
# colnames(combinedYrs)[which(! (colnames(combinedYrs) %in% colnames(Yr1data)))]

combinedYrs <-  merge(x=combinedYrs, y=Yr1data, all=T)
rownames(combinedYrs) <- combinedYrs$Label

# write.csv(mergeDataYr3, file = "../mergeData_freqParent.csv")   # year 3 only
# write.csv(mergeData, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_freqParent_allyrs.csv")  # all years combined


setwd(a) 


# -------------------------------------------------------------------------------------------------------------------------------


#  -------------------------------------  MERGE AND CONCATENATE B CELL DATA --------------------------------------------------
a <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/Analysis/FreqParent")
files <-   list.files("./", pattern=').csv$', full=F)
mergeData <- read.csv(file = "FreqParent.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeData$Batch = "0"

for (i in 1:length(files))
{
  temp <- read.csv(files[i], stringsAsFactors = F, header=T, row.names=1)
  temp$Batch = i
  mergeData <- rbind(mergeData, temp)
}

mergeData <- mergeData[-grep(pattern="Mean", rownames(mergeData)),]
mergeData <- mergeData[-grep(pattern="SD", rownames(mergeData)),]
mergeData <- mergeData[,-which(colnames(mergeData) == "X.1")]

temp <- rownames(mergeData); temp <- str_replace(temp, "PBMCs_", ""); temp <- substr(temp, 1, nchar(temp)-8)   # clean up rownames a bit
rownames(mergeData) <- temp
# add columns for Cohort, Subject, TimePoint (numeric), TimeCategory (categorical), and Year
mergeData$Label <- rownames(mergeData)
mergeData$Subject <- word(mergeData$Label,1,sep = "\\_")  # take first chunk prior to _ character
mergeData$Cohort <- 0;  
mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "19")] <- "aPD1";  
mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "FS")] <- "Healthy"; 
mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "20")] <- "nonPD1"
mergeData$TimePoint <- as.numeric(word(mergeData$Label,2,sep= "d"))  # parse out the number of days since vaccination
mergeData$TimeCategory[which(mergeData$TimePoint == 0)] <- "baseline";  
mergeData$TimeCategory[which(mergeData$TimePoint > 0 & mergeData$TimePoint <= 14)] <- "oneWeek"; 
mergeData$TimeCategory[which(mergeData$TimePoint > 14)] <- "late"
mergeData$Year <- 3

temp <- colnames(mergeData)
temp <- str_replace(temp,"Freq..of.Parent....", "FreqParent")
temp <- str_replace(temp,"Lymphocytes.Singlets.Live.CD14.CD4..CD19..", "CD19hi_")
temp <- str_replace(temp,"Lymphocytes.Singlets.Live.CD4..CD3..", "CD4hi_")
colnames(mergeData) <- temp

# write.csv(mergeData, file = "../mergeData_freqParent.csv")
# write.csv(mergeData, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_freqParent.csv")

setwd(a) 

# -------------------------------------------------------------------------------------------------------------------------------




melted <- melt(mergeData[,-c(119:123)], id.vars = "Label");  melted$value <- as.numeric(melted$value)
ggplot(melted, aes(x=value)) +   geom_histogram(bins = 20) +   facet_wrap(~variable, scales = 'free') + theme_bw()
ggplot(melted, aes(x=variable, y=log(value))) +   geom_boxplot() + theme_bw()
