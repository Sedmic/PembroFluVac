library(cytofkit)
# cytofkit_GUI()

a <- getwd()
setwd(dir = "D:/Pembro-Fluvac/18-19season/Flow cytometry/Tcell/exportFCS/")
file <- list.files("./", pattern='.fcs$', full=TRUE)

file <- file[ -grep("19616-021_d0_004", file)]
# file <- file[ -grep("19616-021_d21_001", file)]
file <- file[ -grep("19616-027_d0_003", file)]
file <- file[ -grep("19616-028_d0_004", file)]
file <- file[ -grep("19616-029_d0_002", file)]
file <- file[ -grep("FS1819-117_d0_003", file)]
file <- file[ -grep("FS1819-126_d7_002", file)]
# file <- file[ -grep("FS1819-121_d28_029", file)] 


parameters <- list.files("./", pattern='.txt$', full=TRUE)
result <- cytofkit(fcsFiles = file, # [c(1:54,56:82)], 
                   markers = parameters,
                   projectName = "cytofkit", 
                   ifCompensation = FALSE,
                   transformMethod = c("logicle"),
                   mergeMethod = c("ceil"), 
                   fixedNum = 1000,
                   dimReductionMethod = c("tsne"),
                   clusterMethods = c("Rphenograph", "FlowSOM"),
                   visualizationMethods = c("tsne"),
                   progressionMethod = c("NULL"), 
                   FlowSOM_k = 50,
                   clusterSampleSize = 500, 
                   resultDir = "./run1/", 
                   saveResults = TRUE,
                   saveObject = TRUE, 
                   openShinyAPP = FALSE)

setwd() <- a
cytofkitShinyAPP()
