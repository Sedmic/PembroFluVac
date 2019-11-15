library(ggplot2)
library(reshape2)
library(stringr)


#  -------------------------------------  MERGE AND CONCATENATE T CELL DATA --------------------------------------------------
a <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Tcell/Analysis/FreqParent")
files <-   list.files("./", pattern=').csv$', full=F)   # this list is all of the duplicates that Windows Explorer appends with (x) postfix to filename
mergeDataYr3 <- read.csv(file = "FreqParent.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr3$TflowBatch = "0"

for (i in 1:length(files))
{
  temp <- read.csv(files[i], stringsAsFactors = F, header=T, row.names=1)
  temp$TflowBatch = i
  mergeDataYr3 <- rbind(mergeDataYr3, temp)
}

mergeDataYr2 <- read.csv(file = "D:/Pembro-Fluvac/17-18season/TcellPanel/Analysis/FreqParent_aPD1.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr2$TflowBatch = i
temp <- read.csv(file = "D:/Pembro-Fluvac/17-18season/TcellPanel/Analysis/FreqParent_healthy.csv", stringsAsFactors = F, header=T, row.names=1)
temp$TflowBatch = i+2
mergeDataYr2 <- rbind(mergeDataYr2, temp)


setwd("D:/Pembro-Fluvac/16-17season/Analysis/freqParent/")
files <-   list.files("./", pattern=').csv$', full=F) # this list is all of the duplicates that Windows Explorer appends with (x) postfix to filename
mergeDataYr1 <- read.csv(file = "byParent-2019edition.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr1$TflowBatch = i+3
for (j in 1:length(files))
{
  temp <- read.csv(files[j], stringsAsFactors = F, header=T, row.names=1)
  temp$TflowBatch =  i + j + 3
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

# write.csv(Yr3data, file = "../mergeData_freqParent.csv")   # year 3 only
# write.csv(combinedYrs, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_freqParent_allyrs.csv")  # all years combined
rm(mergeDataYr3); rm(mergeDataYr2); rm(mergeDataYr1); rm(combinedYrs)


setwd(a) 


# -------------------------------------------------------------------------------------------------------------------------------


#  -------------------------------------  MERGE AND CONCATENATE B CELL DATA --------------------------------------------------
a <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/Analysis/FreqParent")
files <-   list.files("./", pattern=').csv$', full=F)
mergeDataYr3 <- read.csv(file = "FreqParent.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr3$BflowBatch = "0"

for (i in 1:length(files))
{
  temp <- read.csv(files[i], stringsAsFactors = F, header=T, row.names=1)
  temp$BflowBatch = i
  mergeDataYr3 <- rbind(mergeDataYr3, temp)
}
rownames(mergeDataYr3)[grep("-32_d21", rownames(mergeDataYr3))] <- "PBMCs_19616-032_d21_030.fcs"
rownames(mergeDataYr3)[grep("-31_d41", rownames(mergeDataYr3))] <- "PBMCs_19616-031_d41_003.fcs"



# *************************************************************************
# incomplete, need to rationally extract B cell data from year 2 and year 1
# *************************************************************************


cleanColumnNames <- function (mergeData)
{
  mergeData <- mergeData[-grep(pattern="Mean", rownames(mergeData)),]
  mergeData <- mergeData[-grep(pattern="SD", rownames(mergeData)),]
  mergeData <- mergeData[,-which(colnames(mergeData) == "X.1")]
  
  temp <- rownames(mergeData); temp <- str_replace(temp, "PBMCs_", ""); 
  temp <- substr(temp, 1, nchar(temp)-8)   # clean up rownames a bit
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
  
  
  temp <- colnames(mergeData)
  temp <- str_replace(temp,"Freq..of.Parent....", "FreqParent")
  temp <- str_replace(temp,"Lymphocytes.Singlets.Live.CD14.CD4..CD19..", "CD19hi_")
  temp <- str_replace(temp,"Lymphocytes.Singlets.Live.CD4..CD3..", "CD4hi_")
  temp <- str_replace(temp,"CD19hi_NonnaiveB...FreqParent", "CD19hi_NotNaiveB_freqParent")
  temp <- str_replace(temp,"CD19hi_NonnaiveB.", "")
  temp <- str_replace(temp,"IgD.CD71..", "IgDlo_CD71hi_")
  temp <- str_replace(temp,"CD27.CD38..", "CD27hiCD38hi_")
  colnames(mergeData) <- temp
  return (mergeData)
}
Yr3data <- cleanColumnNames(mergeDataYr3);  Yr3data$Year <- 3

# combinedYrs <- merge(x=Yr3data, y=Yr2data, all=T)


# write.csv(Yr3data, file = "../mergeData_freqParent.csv")   # year 3 only
# write.csv(Yr3data, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_freqParent_allyrs.csv")    # all years combined
rm(mergeDataYr3)
setwd(a) 

# -------------------------------------------------------------------------------------------------------------------------------


#  -------------------------------------  MERGE AND CONCATENATE T MFI DATA --------------------------------------------------
a <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Tcell/Analysis/medFI")
files <-   list.files("./", pattern=').csv$', full=F)   # this list is all of the duplicates that Windows Explorer appends with (x) postfix to filename
mergeDataYr3 <- read.csv(file = "medianFI.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr3$TmfiBatch = "0"

for (i in 1:length(files))
{
  temp <- read.csv(files[i], stringsAsFactors = F, header=T, row.names=1)
  temp$TmfiBatch = i
  mergeDataYr3 <- rbind(mergeDataYr3, temp)
}

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
  temp <- str_replace(temp,"Median..Comp.820_60.UV.A", "_medfi-CD20")
  temp <- str_replace(temp,"Median..Comp.800_30.Violet.A", "_medfi-ICOS")
  temp <- str_replace(temp,"Median..Comp.780_60.YG.A", "_medfi-Ki67")
  temp <- str_replace(temp,"Median..Comp.780_60.Red.A", "_medfi-CD36")
  temp <- str_replace(temp,"Median..Comp.750_30.Violet.A", "_medfi-CD4")
  temp <- str_replace(temp,"Median..Comp.740_35.UV.A", "_medfi-CCR6")
  temp <- str_replace(temp,"Median..Comp.730_45.Red.A", "_medfi-CXCR5")
  temp <- str_replace(temp,"Median..Comp.710_50.YG.A", "_medfi-Foxp3")
  temp <- str_replace(temp,"Median..Comp.710_50.Violet.A", "_medfi-CXCR3")
  temp <- str_replace(temp,"Median..Comp.710_50.Blue.A", "_medfi-CD127")
  temp <- str_replace(temp,"Median..Comp.660_20.Red.A", "_medfi-Tcf1")
  temp <- str_replace(temp,"Median..Comp.660_40.YG.A", "_medfi-CXCR4")
  temp <- str_replace(temp,"Median..Comp.660_20.Violet.A", "_medfi-CD25")
  temp <- str_replace(temp,"Median..Comp.660_20.UV.A", "_medfi-CD38")
  temp <- str_replace(temp,"Median..Comp.610_20.YG.A", "_medfi-CD28")
  temp <- str_replace(temp,"Median..Comp.610_20.Violet.A", "_medfi-CD45RA")
  temp <- str_replace(temp,"Median..Comp.586_15.YG.A", "_medfi-aIgG4")
  temp <- str_replace(temp,"Median..Comp.586_15.Violet.A", "_medfi-Dead")
  temp <- str_replace(temp,"Median..Comp.586_15.UV.A", "_medfi-CD19")
  temp <- str_replace(temp,"Median..Comp.515_30.UV.A", "_medfi-CD3")
  temp <- str_replace(temp,"Median..Comp.515_20.Violet.A", "_medfi-CD16")
  temp <- str_replace(temp,"Median..Comp.515_20.Blue.A", "_medfi-SLAM")
  temp <- str_replace(temp,"Median..Comp.470_15.Violet.A", "_medfi-PD1")  
  temp <- str_replace(temp,"Median..Comp.450_50.Violet.A", "_medfi-CTLA4")  
  temp <- str_replace(temp,"Median..Comp.379_28.UV.A", "_medfi-CD27")  
  temp[ grep(x = temp, pattern = "..1$") ] <- str_replace( temp[ grep(x = temp, pattern = "..1$") ], "Time.Lymphocytes.Singlets.Live.CD3..CD4..Nonnaive.CXCR5..ICOS.CD38..cTfh.", 
                       replacement = "cTfh_ICOSloCD38lo_") 
  temp <- str_replace(temp,"Time.Lymphocytes.Singlets.Live.CD3..CD4..Nonnaive.CXCR5..", "cTfh_")
  temp <- str_replace(temp,"Time.Lymphocytes.Singlets.Live.CD19..", "CD19_")
  temp <- str_replace(temp,"Time.Lymphocytes.Singlets.Live.", "Live_")
  temp <- str_replace(temp,"cTfh_ICOS.CD38..","cTfh_ICOShiCD38hi_")
  temp <- str_replace(temp,"CD19_CD27.CD38..","CD19hi_CD27hiCD38hi_")
  temp <- str_replace(temp, "_.._", "_")
  colnames(mergeData) <- temp
  return(mergeData)  
}

Yr3data <- cleanColumnNames(mergeDataYr3);  Yr3data$Year <- 3
combinedYrs <- Yr3data
# write.csv(Yr3data, file = "../mergeData_medianFI.csv")   # year 3 only
# write.csv(combinedYrs, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_medianFI_allyrs.csv")  # all years combined

rm(mergeDataYr3); rm(combinedYrs)

setwd(a) 


# -------------------------------------------------------------------------------------------------------------------------------



#  -------------------------------------  MERGE AND CONCATENATE B MFI DATA --------------------------------------------------
a <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/Analysis/medFI/")
files <-   list.files("./", pattern=').csv$', full=F)
mergeDataYr3 <- read.csv(file = "medianFI.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr3$BmfiBatch = "0"

for (i in 1:length(files))
{
  temp <- read.csv(files[i], stringsAsFactors = F, header=T, row.names=1)
  temp$BmfiBatch = i
  mergeDataYr3 <- rbind(mergeDataYr3, temp)
}
rownames(mergeDataYr3)[grep("-32_d21", rownames(mergeDataYr3))] <- "PBMCs_19616-032_d21_030.fcs"
rownames(mergeDataYr3)[grep("-31_d41", rownames(mergeDataYr3))] <- "PBMCs_19616-031_d41_003.fcs"


# *************************************************************************
# incomplete, need to rationally extract B cell data from year 2 and year 1
# *************************************************************************

# mergeData <- mergeDataYr3 

cleanColumnNames <- function (mergeData)
{
  mergeData <- mergeData[-grep(pattern="Mean", rownames(mergeData)),]
  mergeData <- mergeData[-grep(pattern="SD", rownames(mergeData)),]
  mergeData <- mergeData[,-which(colnames(mergeData) == "X.1")]
  
  temp <- rownames(mergeData); temp <- str_replace(temp, "PBMCs_", ""); 
  temp <- substr(temp, 1, nchar(temp)-8)   # clean up rownames a bit
  rownames(mergeData) <- temp
  # add columns for Cohort, Subject, TimePoint (numeric), TimeCategory (categorical), and Year
  mergeData$Label <- rownames(mergeData)
  mergeData$Subject <- word(mergeData$Label,1,sep = "\\_")  # take first chunk prior to _ character
  mergeData$Cohort <- 0;  
  mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "19")] <- "aPD1";  
  mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "FS")] <- "Healthy"; 
  mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "20")] <- "nonPD1"
  mergeData$TimePoint <- as.numeric(word(mergeData$Label,2,sep= "d"))  # parse out the number of days since vaccination
  mergeData$TimeCategory <- "a"
  mergeData$TimeCategory[which(mergeData$TimePoint == 0)] <- "baseline";  
  mergeData$TimeCategory[which(mergeData$TimePoint > 0 & mergeData$TimePoint <= 14)] <- "oneWeek"; 
  mergeData$TimeCategory[which(mergeData$TimePoint > 14)] <- "late"
  
  
  temp <- colnames(mergeData)
  temp <- str_replace(temp,"Median..", "medfi-")
  temp <- str_replace(temp,"Lymphocytes.Singlets.Live.CD14.CD4..CD19..", "CD19_")
  temp <- str_replace(temp,"CD19_NonnaiveB.IgD.CD71..ActBCells", "ActBCells")
  temp <- str_replace(temp,"CD19_NonnaiveB.IgD.CD71..CD20loCD71hi", "CD20loCD71hi")
  temp <- str_replace(temp,"CD19_NonnaiveB.IgD.CD71..", "IgDloCD71hi")
  
  
  colnames(mergeData) <- temp
  return(mergeData)  
}
Yr3data <- cleanColumnNames(mergeDataYr3);  Yr3data$Year <- 3
which(sapply(Yr3data, class) == "character")
# combinedYrs <- merge(x=Yr3data, y=Yr2data, all=T)

# write.csv(Yr3data, file = "../mergeData_Bmfi.csv")   # year 3 only
# write.csv(Yr3data, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_Bmfi_allyrs.csv")    # all years combined
rm(mergeDataYr3)
setwd(a) 

# -------------------------------------------------------------------------------------------------------------------------------





#  -------------------------------------  MERGE AND CONCATENATE T from B CELL PANEL --------------------------------------------------
a <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/Analysis/TfhAnalyses/")
files <-   list.files("./", pattern=').csv$', full=F)
mergeDataYr3 <- read.csv(file = "TfhAnalyses.csv", stringsAsFactors = F, header = T, row.names = 1)
mergeDataYr3$TfromBBatch = "0"

for (i in 1:length(files))
{
  temp <- read.csv(files[i], stringsAsFactors = F, header=T, row.names=1)
  temp$TfromBBatch = i
  mergeDataYr3 <- rbind(mergeDataYr3, temp)
}
rownames(mergeDataYr3)[grep("-32_d21", rownames(mergeDataYr3))] <- "PBMCs_19616-032_d21_030.fcs"
rownames(mergeDataYr3)[grep("-31_d41", rownames(mergeDataYr3))] <- "PBMCs_19616-031_d41_003.fcs"

# *************************************************************************
# incomplete, need to rationally extract B cell data from year 2 and year 1
# *************************************************************************

# mergeData <- mergeDataYr3 

cleanColumnNames <- function (mergeData)
{
  mergeData <- mergeData[-grep(pattern="Mean", rownames(mergeData)),]
  mergeData <- mergeData[-grep(pattern="SD", rownames(mergeData)),]
  mergeData <- mergeData[,-which(colnames(mergeData) == "X.1")]
  
  temp <- rownames(mergeData); temp <- str_replace(temp, "PBMCs_", ""); 
  temp <- substr(temp, 1, nchar(temp)-8)   # clean up rownames a bit
  rownames(mergeData) <- temp
  # add columns for Cohort, Subject, TimePoint (numeric), TimeCategory (categorical), and Year
  mergeData$Label <- rownames(mergeData)
  mergeData$Subject <- word(mergeData$Label,1,sep = "\\_")  # take first chunk prior to _ character
  mergeData$Cohort <- 0;  
  mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "19")] <- "aPD1";  
  mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "FS")] <- "Healthy"; 
  mergeData$Cohort[which(substr(mergeData$Label, 1, 2) == "20")] <- "nonPD1"
  mergeData$TimePoint <- as.numeric(word(mergeData$Label,2,sep= "d"))  # parse out the number of days since vaccination
  mergeData$TimeCategory <- "a"
  mergeData$TimeCategory[which(mergeData$TimePoint == 0)] <- "baseline";  
  mergeData$TimeCategory[which(mergeData$TimePoint > 0 & mergeData$TimePoint <= 14)] <- "oneWeek"; 
  mergeData$TimeCategory[which(mergeData$TimePoint > 14)] <- "late"
  
  
  temp <- colnames(mergeData)
  temp <- str_replace(temp,"Median..Comp.820_60.UV.A", "_medfi-CD20")
  temp <- str_replace(temp,"Median..Comp.800_30.Violet.A", "_medfi-Ki67")
  temp <- str_replace(temp,"Median..Comp.780_60.YG.A", "_medfi-Tbet")
  temp <- str_replace(temp,"Median..Comp.780_60.Red.A", "_medfi-CD80")
  temp <- str_replace(temp,"Median..Comp.780_60.Blue.A", "_medfi-CD23")
  temp <- str_replace(temp,"Median..Comp.750_30.Violet.A", "_medfi-CD4")
  temp <- str_replace(temp,"Median..Comp.740_35.UV.A", "_medfi-CD16")
  temp <- str_replace(temp,"Median..Comp.730_45.Red.A", "_medfi-CXCR5")
  temp <- str_replace(temp,"Median..Comp.710_50.YG.A", "_medfi-CD14")
  temp <- str_replace(temp,"Median..Comp.710_50.Violet.A", "_medfi-CD21")
  temp <- str_replace(temp,"Median..Comp.710_50.Blue.A", "_medfi-CD19")
  temp <- str_replace(temp,"Median..Comp.660_20.Red.A", "_medfi-Blimp1")
  temp <- str_replace(temp,"Median..Comp.660_40.YG.A", "_medfi-CXCR4")
  temp <- str_replace(temp,"Median..Comp.660_20.Violet.A", "_medfi-CD71")
  temp <- str_replace(temp,"Median..Comp.660_20.UV.A", "_medfi-CD11c")
  temp <- str_replace(temp,"Median..Comp.610_20.YG.A", "_medfi-CD38")
  temp <- str_replace(temp,"Median..Comp.610_20.Violet.A", "_medfi-CD138")
  temp <- str_replace(temp,"Median..Comp.586_15.YG.A", "_medfi-CD86")
  temp <- str_replace(temp,"Median..Comp.586_15.Violet.A", "_medfi-Dead")
  temp <- str_replace(temp,"Median..Comp.586_15.UV.A", "_medfi-CD25")
  temp <- str_replace(temp,"Median..Comp.515_30.UV.A", "_medfi-CD3")
  temp <- str_replace(temp,"Median..Comp.515_20.Violet.A", "_medfi-IgM")
  temp <- str_replace(temp,"Median..Comp.515_20.Blue.A", "_medfi-CD32")
  temp <- str_replace(temp,"Median..Comp.470_15.Violet.A", "_medfi-IgD")  
  temp <- str_replace(temp,"Median..Comp.450_50.Violet.A", "_medfi-Bcl6")  
  temp <- str_replace(temp,"Median..Comp.379_28.UV.A", "_medfi-CD27")  
    
  
  temp <- str_replace(temp,"Lymphocytes.Singlets.Live.CD4..CD3..CXCR5..CD38.Ki67....", "Ki67hiCD38hicTfh")
  colnames(mergeData) <- temp
  return(mergeData)  
}
Yr3data <- cleanColumnNames(mergeDataYr3);  Yr3data$Year <- 3

# combinedYrs <- merge(x=Yr3data, y=Yr2data, all=T)


# write.csv(Yr3data, file = "../mergeData_TfhAnalyses.csv")   # year 3 only
# write.csv(Yr3data, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_TmfiFromB_allyrs.csv")    # all years combined
rm(mergeDataYr3)
setwd(a) 

# -------------------------------------------------------------------------------------------------------------------------------




















