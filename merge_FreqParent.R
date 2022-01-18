library(reshape2)
library(stringr)
sessionInfo()

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
cleanColumnNames <- function (df)
{
  df <- df[-grep(pattern="Mean", rownames(df)),]
  df <- df[-grep(pattern="SD", rownames(df)),]
  df <- df[,-which(colnames(df) == "X.1")]
  
  temp <- rownames(df); temp <- str_replace(temp, "PBMCs_", ""); temp <- str_replace(temp, "PBMC_", ""); 
  temp <- str_replace(temp, "day", "d"); temp <- substr(temp, 1, nchar(temp)-8)   # clean up rownames a bit, lop off the last 8 characters _0xx.fcs
  rownames(df) <- temp
  # add columns for Cohort, Subject, TimePoint (numeric), TimeCategory (categorical), and Year
  df$Label <- rownames(df)
  df$Subject <- word(df$Label,1,sep = "\\_")  # take first chunk prior to _ character
  df$Cohort <- 0;  
  df$Cohort[which(substr(df$Label, 1, 2) == "19")] <- "aPD1";  
  df$Cohort[which(substr(df$Label, 1, 2) == "FS")] <- "Healthy"; 
  df$Cohort[which(substr(df$Label, 1, 2) == "20")] <- "nonPD1"
  df$TimePoint <- as.numeric(word(df$Label, 2, sep= "d"))  # parse out the number of days since vaccination
  df$TimeCategory <- 0;  
  df$TimeCategory[which(df$TimePoint == 0)] <- "baseline"  
  df$TimeCategory[which(df$TimePoint > 0 & df$TimePoint <= 14)] <- "oneWeek"; 
  df$TimeCategory[which(df$TimePoint > 14)] <- "late"

  temp <- colnames(df)
  temp <- str_replace(temp,"Freq..of.Parent....", "FreqParent")
  temp <- str_replace(temp,"Time.Lymphocytes.Singlets.Live.CD3..CD4..Nonnaive.CXCR5..", "cTfh_")
  temp <- str_replace(temp,"Time.Lymphocytes.Singlets.Live.CD19..", "CD19_")
  temp <- str_replace(temp,"Time.Lymphocytes.Singlets.Live.", "Live_")
  # temp[grep(pattern = "FreqParent.1", temp)] <- str_replace( temp[grep(pattern = "FreqParent.1", temp)],"cTfh_ICOS.CD38..","cTfh_ICOSloCD38lo_") 
  # temp <- str_replace( temp,"cTfh_ICOSloCD38lo.","cTfh_ICOSloCD38lo_") 
  # temp <- str_replace(temp,"cTfh_ICOS.CD38..","cTfh_ICOShiCD38hi_")
  temp[grep(pattern = "Live_CD3..CD4....FreqParent.1", temp)] <- str_replace( temp[grep(pattern = "Live_CD3..CD4....FreqParent.1", temp)], "Live_CD3..CD4.", "Live_CD3hi_CD4lo_" )
  temp <- str_replace(temp, paste0( "FreqParent.", seq(1:5), collapse="|"),"FreqParent")
  temp <- str_replace(temp,"cTfh\\.","")
  colnames(df) <- temp
  return(df)  
}

Yr3data <- cleanColumnNames(mergeDataYr3);  Yr3data$Year <- 3
Yr2data <- cleanColumnNames(mergeDataYr2);  Yr2data$Year <- 2
Yr1data <- cleanColumnNames(mergeDataYr1);  Yr1data$Year <- 1

names(Yr3data)[22:40] <- str_replace( names(Yr3data)[22:40],"cTfh_ICOS.CD38..","cTfh_ICOShiCD38hi_") 
names(Yr3data)[41:59] <- str_replace( names(Yr3data)[41:59],"cTfh_ICOS.CD38..","cTfh_ICOSloCD38lo_") 
names(Yr3data)[c(38,57)] <- str_replace(names(Yr3data)[c(38,57)], "CXCR3.CCR6.", "CXCR3lo_CCR6lo")
names(Yr3data)[c(39,58)] <- str_replace(names(Yr3data)[c(39,58)], "CXCR3.CCR6.", "CXCR3lo_CCR6hi")
names(Yr3data)[c(40,59)] <- str_replace(names(Yr3data)[c(40,59)], "CXCR3.CCR6.", "CXCR3hi_CCR6lo")

names(Yr2data)[22:37] <- str_replace( names(Yr2data)[22:37],"cTfh_ICOS.CD38..","cTfh_ICOShiCD38hi_") 
names(Yr2data)[38:53] <- str_replace( names(Yr2data)[38:53],"cTfh_ICOS.CD38..","cTfh_ICOSloCD38lo_") 

names(Yr1data)[19:30] <- str_replace( names(Yr1data)[19:30],"cTfh_ICOS.CD38..","cTfh_ICOShiCD38hi_") 
# names(Yr1data)[41:59] <- str_replace( names(Yr1data)[41:59],"cTfh_ICOS.CD38..","cTfh_ICOSloCD38lo_") 


combinedYrs <- merge(x=Yr3data, y=Yr2data, all=T)
combinedYrs <-  merge(x=combinedYrs, y=Yr1data, all=T)
rownames(combinedYrs) <- combinedYrs$Label


# write.csv(Yr3data, file = "../mergeData_freqParent.csv")   # year 3 only
write.csv(combinedYrs, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Tflow_freqParent_allyrs.csv")  # all years combined
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

# write.csv(Yr3data, file = "D:/Pembro-Fluvac/Analysis/mergedData/mergeData_Bflow_TmfiFromB_allyrs.csv")    # all years combined
rm(mergeDataYr3)
setwd(a) 

# -------------------------------------------------------------------------------------------------------------------------------

#  -------------------------------------  MERGE AND CONCATENATE BAA YEAR 3 DATA  --------------------------------------------------

week1 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week1/rherati_20130914_A1V1-091413_Senescence/Analysis/ICOShiCD38hi.csv")
week2 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week2/rherati_20130920_Week2_Senescence/Analysis/ICOShiCD38hi.csv")
week3 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week3/rherati_20130925_Week3_Senescence/Analysis/ICOShiCD38hi.csv")
week4 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week4/rherati_20131002_Week4_Senescence/Analysis/ICOShiCD38hi.csv")
week5 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week5/rherati_20131009_Week5_Senescence/Analysis/ICOShiCD38hi.csv")
week6 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week6/rherati_20131017_Week6_Senescence/Analysis/ICOShiCD38hi.csv")
week7 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week7/rherati_20131023_Week7_Senescence/Analysis/ICOShiCD38hi.csv")
  week7$X[10] <- "222-055-V1_0.fcs"   # mislabeling of the fcs file
week8 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week8/rherati_20131030_Week8_Senescence/Analysis/ICOShiCD38hi.csv")
week9 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week9/rherati_20131108_Week9_Senescence/Analysis/ICOShiCD38hi.csv")
week10 <- read.csv(file="D:/BAA project/Year3/Data files by week/Aim2 - Senescence and TransFactor/Week10/rherati_20131115_Week10_Senescence/Analysis/ICOShiCD38hi.csv")

BAAyear3long <- rbind(week1,week2,week3,week4,week5,week6,week7,week8,week9,week10); BAAyear3long$X.1 <- NULL; names(BAAyear3long) <- c("file","ICOShiCD38hi_FreqParent")
BAAyear3long <- BAAyear3long[-grep("222-053",BAAyear3long$file), ]        # exclude 222-053 because of poor quality CXCR5 stain at visit 1
BAAyear3long$Subject <- substr(BAAyear3long$file, 1, 7);  BAAyear3long <- BAAyear3long[-grep(paste0(c("Mean","SD"), collapse = "|"), BAAyear3long$file),]
BAAyear3long$Cohort <- "Young"
BAAyear3long$Cohort[grepl(pattern = "222",x = BAAyear3long$Subject)] <- "Elderly"
BAAyear3long$Visit <- substr(BAAyear3long$file, 9,10)
BAAyear3wide <- dcast(BAAyear3long, `Subject` + `Cohort` ~ `Visit`, value.var = c("ICOShiCD38hi_FreqParent"))
BAAyear3wide$FCtfh <- BAAyear3wide$V2 / BAAyear3wide$V1 

# write.csv(BAAyear3long, file = "D:/Pembro-Fluvac/Analysis/mergedData/BAAyear3long.csv")    
# write.csv(BAAyear3wide, file = "D:/Pembro-Fluvac/Analysis/mergedData/BAAyear3wide.csv")    


rm(week1,week2,week3,week4,week5,week6,week7,week8,week9,week10)




seroTemp <- read.csv(file= "D:/Pembro-Fluvac/Analysis/mergedData/20200711 SerologyMeasurements.csv", stringsAsFactors = F, header = T)
seroRock <- read.csv(file="D:/Pembro-Fluvac/Analysis/mergedData/20201217 Rockefeller_Serology.csv")
seroUpdate <- read.csv(file = "D:/Pembro-Fluvac/Analysis/Rockefeller/122721 Flu PD-1 Combined seasons.csv"); seroUpdate <- seroUpdate[,1:10]
names(seroUpdate)[1] <- "Subject";

seroMerge <- merge(x=seroTemp, y=seroUpdate, by = c("Subject","Updated.RU.identifier"), all.x = T)
plot(x=seroMerge$Plasma.CXCL13..pg.mL.,y=seroMerge$CXCL13.combined.repeat.2020)
# cor(x=seroMerge$Plasma.CXCL13..pg.mL.,y=seroMerge$CXCL13.combined.repeat.2020, use="complete.obs")
# 
# write.csv(seroMerge[, -which(names(seroMerge) %in% c("Day","Treatment","CXCL13.combined.repeat.2020","HAI.titer..H1N1pdm09."))],
          # file = "D:/Pembro-Fluvac/Analysis/mergedData/SerologyMeasurements.csv")

names(seroUpdate)[1] <- "Label";
seroMerge <- merge(x=seroRock, y=seroUpdate, by = c("Label","Updated.RU.identifier"), all.x = T)
plot(x=seroMerge$Plasma.CXCL13..pg.mL.,y=seroMerge$CXCL13.combined.repeat.2020)
# cor(x=seroMerge$Plasma.CXCL13..pg.mL.,y=seroMerge$CXCL13.combined.repeat.2020, use="complete.obs")
#
# write.csv(seroMerge[, -which(names(seroMerge) %in% c("Day","Treatment","CXCL13.combined.repeat.2020","HAI.titer..H1N1pdm09."))], 
#           file = "D:/Pembro-Fluvac/Analysis/mergedData/Rockefeller_Serology.csv")

rm(seroTemp, seroUpdate, seroMerge, seroRock)

