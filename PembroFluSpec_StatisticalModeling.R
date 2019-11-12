
library("MASS")
library("glmnet")

mergedData <- read.csv(file = "D:/Pembro-Fluvac/Analysis/mergedData/allMergedData.csv", header=T, stringsAsFactors = F)

## ------------------------------------------------- Systematic eval of differences  ----------------------------------------------------------

subsetData <- subset(mergedData, TimeCategory == "baseline")
summary( glm(Cohort ~ . , family = binomial(link='logit'), data = subsetData[ , c(grep("Cohort",colnames(subsetData),value = F), 150:160)]) )   # repeated in batches of ~10 columns
subsetData <- subset(mergedData, TimeCategory == "oneWeek")
summary( glm(Cohort ~ . , family = binomial(link='logit'), data = subsetData[ , c(grep("Cohort",colnames(subsetData),value = F), 150:160)]) )   # repeated in batches of ~10 columns


# need to do a regularization approach here because more variables than observations

#   BASELINE 

subsetData <- subset(mergedData, TimeCategory == "baseline" & Cohort != "nonPD1")
subsetData$Cohort <- factor(subsetData$Cohort, levels = c("Healthy", "aPD1"))
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
a <- a[ -which(is.na(rowSums(a[,-1]) ) ), ]
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

