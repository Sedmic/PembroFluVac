
library("MASS")
library("glmnet")

mergedData <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/allMergedData.csv", header=T, stringsAsFactors = F)

## ------------------------------------------------- Systematic eval of differences  ----------------------------------------------------------

subsetData <- subset(mergedData, TimeCategory == "baseline")
summary( glm(Cohort ~ . , family = binomial(link='logit'), data = subsetData[ , c(grep("Cohort",colnames(subsetData),value = F), 150:160)]) )   # repeated in batches of ~10 columns
subsetData <- subset(mergedData, TimeCategory == "oneWeek")
summary( glm(Cohort ~ . , family = binomial(link='logit'), data = subsetData[ , c(grep("Cohort",colnames(subsetData),value = F), 150:160)]) )   # repeated in batches of ~10 columns


# need to do a regularization approach here because more variables than observations

#   BASELINE  differences between cohorts    CATEGORICAL

subsetData <- subset(mergedData, TimeCategory == "baseline" & Cohort != "nonPD1")
exclude <- grep (paste(c("Subject","TimeCategory","Year","Label","TimePoint.one","Sex","dummy", "Current.Immunotherapy","Cycle.of.Immunotherapy"), collapse="|"), colnames(subsetData), value = F)
subsetData <- subsetData[, -exclude] 
# subsetData <- subsetData[, -c(100:200)]
a <- sapply(subsetData, function (x) { sum(!is.na(x)) } )   # take away any columns that are purely NA
subsetData <- subsetData[, which(a > 2)]
a <- sapply(subsetData[,2:ncol(subsetData)], function (x) { sum(x, na.rm=T) } )  # take away all zero columns
subsetData <- subsetData[, - c(1+which(a == 0, arr.ind = T))]  # add one because the Cohort column was excluded, so a doesn't represent subsetData otherwise
a <- sapply(subsetData[,2:ncol(subsetData)], function (x) { sd(x, na.rm=T) } )   # take away all columns with 0 standard deviation
subsetData <- subsetData[, - c(1+which(a == 0, arr.ind = T))]  # add one because the Cohort column was excluded, so a doesn't represent subsetData otherwise
a <- subsetData[ , -grep(paste(c("Kd","Plasma","HAI"), collapse="|"),colnames(subsetData), value = F)]   # leave out the serologies columns because not enough paired data  
exclude <- which(is.na(rowSums(a[,-1], na.rm=T) ) ) 
if (length(exclude) >0)    {   a <- a[ -exclude, ] }  # exclude any completely blank rows 
set.seed(101)
train_rows <- sample(1:nrow(a), 0.66*nrow(a))
y.train <- a$Cohort[train_rows]
x.train <- as.matrix(a[train_rows, -grep("Cohort",colnames(a), value = F)])
ytest <- a$Cohort[-train_rows]
xtest <- as.matrix(a[-train_rows, -grep("Cohort",colnames(a), value = F)])
list.of.fits <- list()
for (i in 0:100) {   fit.name <- paste0("alpha",i/100); list.of.fits[[fit.name]] <- cv.glmnet(x.train, y.train, alpha=i/100, family="binomial", standardize=T, nfolds=10) }
results <- data.frame()
for (i in 0:100) {   
  fit.name <- paste0("alpha",i/100);   
  fit.predicted <- predict(list.of.fits[[fit.name]], s=list.of.fits[[fit.name]]$lambda.1se, newx = xtest, type="class");
  accuracy <- mean(ytest == fit.predicted) ;  
  temp <- data.frame(alpha=i/100, accuracy=accuracy);   
  results <- rbind(results, temp)    }
ggplot(results, aes(x=alpha, y=accuracy)) + geom_point() + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), title = element_text(size=24)) + 
  ggtitle("Accuracy for prediction of test data given model")
fit <- cv.glmnet(x.train, y.train, alpha=1, family = "binomial", standardize=T, nfolds = 10) 
coef(fit, fit$lambda.1se)






#   DAY 7 as prediction for CONTINUOUS              
runElasticNet <- function( data, testContinuous, oneWeekFlag="T")
{
  testColumn <- testContinuous
  resultObject <- list()
  set.seed(102)
  if (oneWeekFlag == 'T') {  subsetData <- subset(data, TimeCategory == "oneWeek" & Cohort != "nonPD1" & Year == "3") }
  if (oneWeekFlag != 'T') { subsetData <- subset(data, TimeCategory == "late" & Cohort != "nonPD1" ) }
  exclude <- grep (paste(c("Subject","TimeCategory","Year","Label","Cohort","TimePoint.one","Sex","dummy", "Current.Immunotherapy","Cycle.of.Immunotherapy"), collapse="|"), colnames(subsetData), value = F)
  subsetData <- subsetData[, -exclude] 
  a <- sapply(subsetData, function (x) { sum(!is.na(x)) } )   # take away any columns that are purely NA
  subsetData <- subsetData[, which(a > 2)]
  a <- sapply(subsetData , function (x) { sum(x, na.rm=T) } )  # take away all zero columns
  subsetData <- subsetData[, - c(which(a == 0, arr.ind = T))]  # 
  a <- sapply(subsetData , function (x) { sd(x, na.rm=T) } )   # take away all columns with 0 standard deviation
  subsetData <- subsetData[, - c(which(a == 0, arr.ind = T))]  # 
  exclude <- which(is.na(rowSums(subsetData, na.rm=T) ) )      # exclude any completely blank rows 
  if (length(exclude) >0)    {   subsetData <- subsetData[ -exclude, ] }  
  a <- subsetData;   
  train_rows <- sample(1:nrow(subsetData), 0.66*nrow(subsetData))
  train_rows <- train_rows[ !is.na(a[train_rows, testColumn] ) ]     # exclude rows where the y response variable is NA 
  y.train <- a[train_rows, testColumn]
  x.train <- as.matrix(a[train_rows, -grep(testColumn,colnames(a), value = F)])
  ytest <- a[-train_rows, testColumn]
  xtest <- as.matrix(a[-train_rows, -grep(testColumn,colnames(a), value = F)])
  list.of.fits <- list();  results <- data.frame()
  for (i in 0:50) {   
    fit.name <- paste0("alpha",i/50)
    list.of.fits[[fit.name]] <- cv.glmnet(x.train[complete.cases(x.train) , ], y.train[complete.cases(x.train)], alpha=i/50, standardize=T, nfolds=10)      }
  for (i in 0:50) {   
    fit.name <- paste0("alpha",i/50);   
    fit.predicted <- predict(list.of.fits[[fit.name]], s=list.of.fits[[fit.name]]$lambda.1se, newx = xtest[complete.cases(xtest), ]);
    MSE <- mean((fit.predicted - ytest[complete.cases(xtest)])^2)  
    temp <- data.frame(alpha=i/50, MSE=MSE);   
    results <- rbind(results, temp)    }
  resultObject[[1]] <- ggplot(results, aes(x=alpha, y=MSE)) + geom_point() + theme_bw() + theme(axis.title = element_text(size=20), axis.text = element_text(size=16), title = element_text(size=24)) + 
    ggtitle("Accuracy for prediction of test data given model")
  # now for the final model using the whole dataset
  z <- rownames(results)[which(results$MSE == min(results$MSE))]
  a <- subsetData
  fit <- cv.glmnet(x = as.matrix(a[complete.cases(a), -grep(testColumn,colnames(a), value = F)]), 
                   y = as.matrix(a[complete.cases(a), testColumn]), alpha= results[z[1],"alpha"] + 0.01, standardize=T, nfolds = 10) 
  resultObject[[2]] <- fit
  predictors <- as.matrix(coef(fit, fit$lambda.1se)); 
  resultObject[[3]] <- predictors[order(predictors, decreasing = F),]
  return(resultObject)
}
# very few biological replicates compared to the number of columns in the dataframe
# only terms with nonzero coefficients are 

# IgDloCD71hi..medfi.CD86.
# CD19_CD27.CD38....FreqParent
# CD27hiCD38hi_..FreqParent

# so nothing much new gained by looking at medfi or other subsets

runElasticNet(data = mergedData, testContinuous = "cTfh_ICOShiCD38hi_..FreqParent")[[3]]
runElasticNet(data = mergedData, testContinuous = "CD20loCD71hi...medfi.CD27.")[[3]]

subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1")
runElasticNet(data = subsetData, testContinuous = "FCtfh_oW", oneWeekFlag = "T")[[3]]
subsetData <- subset(mergedData, TimeCategory != "late" & Cohort != "nonPD1")
runElasticNet(data = subsetData, testContinuous = "FCPB_oW", oneWeekFlag = "T")[[3]]
subsetData <- subset(mergedData, TimeCategory != "oneWeek" & Cohort != "nonPD1")
runElasticNet(data = subsetData, testContinuous = "FChai_late", oneWeekFlag="F")[[3]]



