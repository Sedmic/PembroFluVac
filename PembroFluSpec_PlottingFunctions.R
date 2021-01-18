library("grid")
library("ggplot2")
library("gplots")
library("viridis")
library("reshape2")
# library("ggpubr")
# library("ggcorrplot")
library("gridExtra")
library("RColorBrewer")
library("scales")
library("MASS")
library("ggrepel")
library("pheatmap")
library("e1071")

## ------------------------------------------- PLOTTING FUNCTIONS ------------------------------------------------------------------------------


twoSampleBox <- function (data, xData, yData, fillParam, title, yLabel, position = "left")
{
  fit <- lm(data[,yData] ~ data[,xData])
  pValue <- summary(fit)$coefficients[2,"Pr(>|t|)"];  CI <- confint(fit)[2,]; CI <- round(CI,2)

  if (pValue < 0.01)
  {
    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1), "\n", "Mean 95%CI: [",CI[1], ", ",CI[2], "]")
  }
  if (pValue >= 0.01)
  {
    annotationInfo <- paste0("P = ", round(pValue, 2),"\n", "Mean 95%CI: \n[",CI[1], ", ",CI[2], "]")
  }
  if(position == "left")  { my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=24)))   }
  if(position == "right")  { my_grob = grobTree(textGrob(annotationInfo, x=0.73,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=24)))   }
  
  return (
    ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) +  
      geom_boxplot(outlier.shape = NA) + geom_jitter(size=7, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
      ggtitle(title) + ylab(yLabel) +  
      theme(axis.text = element_text(size=28,hjust = 0.5), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=32,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")     )
}

twoSampleBarMelted <- function (data, xData, yData, fillParam, title, yLabel)
{
  set.seed(101)
  fit <- lm(data[,yData] ~ data[,xData])
  pValue <- summary(fit)$coefficients[2,"Pr(>|t|)"];  CI <- confint(fit)[2,]; CI <- round(CI,2)
  if (pValue < 0.01)
  {
    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1), "\n", "Mean 95%CI: [",CI[1], ", ",CI[2], "]")
  }
  if (pValue >= 0.01)
  {
    annotationInfo <- paste0("P = ", round(pValue, 2),"\n", "Mean 95%CI: \n[",CI[1], ", ",CI[2], "]")
  }
  my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=28)))
  
  if( ! is.factor(data[,xData])) {  data[,xData] <- factor(data[,xData])    }
  for (i in 1:length(levels(data[,xData])))
  {
    data[, paste(levels(data[,xData])[i])] <- mean(data[which(data[,xData] == levels(data[,xData])[i]), yData], na.rm=T)   # now mean calculated for each level of xData
  }
  
  overTime <- aggregate( x = data[,yData], by= list(variable = data$variable, TimeCategory = data$TimeCategory, Cohort = data$Cohort), FUN=mean, na.rm = T)
  names(overTime)[which(names(overTime) == 'x')] <- yData
  
  # ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) +  
  #  stat_summary(data=data, fun.y = 'mean', geom="bar")  + geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.1)) + 
  #  ggtitle(title) + ylab(yLabel) +  theme_bw() +  theme(axis.text = element_text(size=28,hjust = 0.5), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=32,hjust = 0.5)) + 
  #  annotation_custom(my_grob) + theme(legend.position = "none") 
  
  return (
    ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) +  
      geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity") + 
      geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.1)) + 
      ggtitle(title) + ylab(yLabel) +  theme_bw() +
      theme(axis.text = element_text(size=28,hjust = 0.5), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=32,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")
  )
}

twoSampleBar <- function (data, xData, yData, fillParam, title, yLabel, batch="none", position = "left", FCplot=F, confInt=F, ttest=T, nonparam=F)
{
  set.seed(102)
  if (ttest == T && batch == "none" && nonparam == F)
  {
    justforttest <- data[, c(xData,yData)]
    fit <- rstatix::t_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
    pValue <- fit$p
  }
  if (ttest == T && batch == "none" && nonparam == T)
  {
    justforttest <- data[, c(xData,yData)]
    fit <- rstatix::wilcox_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
    pValue <- fit$p
  }
  if (ttest == F)
  {
    if (batch == "none")  {
      fit <- wilcox_test(data = data, as.formula(paste(yData, xData, sep="~") ));  pValue <- fit$p  }
    if (batch != "none")  {  fit <- lm(data[,yData] ~ data[,xData] + data[,batch]); pValue <- summary(fit)$coefficients[2,"Pr(>|t|)"];  CI <- confint(fit)[2,]; CI <- round(CI,2)    }
  }
  if (! is.nan(pValue) && confInt == T)
  {
    if (pValue < 0.01)
      {    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1), "\n", "Mean 95%CI: [",CI[1], ", ",CI[2], "]")     }
    if (pValue >= 0.01)
    {    annotationInfo <- paste0("P = ", round(pValue, 2),"\n", "Mean 95%CI: \n[",CI[1], ", ",CI[2], "]")      }
    if(position == "left")  { my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))   }
    if(position == "right")  { my_grob = grobTree(textGrob(annotationInfo, x=0.65,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))   }
    if(position == "none") { my_grob = grobTree(textGrob(annotationInfo, x=10,  y=10, hjust=0, gp=gpar(col="black", fontsize=1)))   }
  }
  if (! is.nan(pValue) && confInt == F)
  {
    if (pValue < 0.01)
    {    annotationInfo <- paste0("P = ", formatC(pValue, format="e",digits=1))     }
    if (pValue >= 0.01)
    {    annotationInfo <- paste0("P = ", round(pValue, 2))      }
    if(position == "left")  { my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28)))   }
    if(position == "right")  { my_grob = grobTree(textGrob(annotationInfo, x=0.65,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28)))   }
    if(position == "none") { my_grob = grobTree(textGrob(annotationInfo, x=10,  y=10, hjust=0, gp=gpar(col="black", fontsize=1)))   }
    
  }

  if( ! is.factor(data[,xData])) {  data[,xData] <- factor(data[,xData])    }
  # for (i in 1:length(levels(data[,xData]))) { data[, paste(levels(data[,xData])[i])] <- mean(data[which(data[,xData] == levels(data[,xData])[i]), yData], na.rm=T)}   # now mean calculated for each level of xData
  
  overTime <- aggregate(x = data[,yData], by= list( data[,xData]), FUN=mean, na.rm = T)
  names(overTime)[which(names(overTime) == 'Group.1')] <- xData
  names(overTime)[which(names(overTime) == 'x')] <- yData
  
  # ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) +  
  #  stat_summary(data=data, fun.y = 'mean', geom="bar")  + geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.1)) + 
  #  ggtitle(title) + ylab(yLabel) +  theme_bw() +  theme(axis.text = element_text(size=28,hjust = 0.5), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=32,hjust = 0.5)) + 
  #  annotation_custom(my_grob) + theme(legend.position = "none") 
  
  if (FCplot == T)
  return (
    ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) + 
      geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5) + 
      geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity", color='black',size=0.1) + 
      geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.15)) + 
      ggtitle(title) + ylab(yLabel) +  theme_bw() +
      theme(axis.text = element_text(size=28,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), 
            plot.title = element_text(size=32,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")
  )
  if (FCplot == F)
    return (
      ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.6)) + scale_fill_manual(values = c("#7FAEDB", "#FFB18C")) + 
        geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity", color='black',size=0.1) + 
        geom_point(size=7, pch=21, fill="black", color="white", alpha=0.5, position = position_jitter(width=0.15)) + 
        ggtitle(title) + ylab(yLabel) +  theme_bw() +
        theme(axis.text = element_text(size=28,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), axis.title.x = element_blank(), 
              plot.title = element_text(size=32,hjust = 0.5)) + 
        annotation_custom(my_grob) + theme(legend.position = "none")
    )
}

twoSampleBarRock <- function(data, xData, yData, fillParam, title, yLabel, position="left")
{
  overTime <- aggregate(x = data[,yData], by= list( data[,xData]), FUN=mean, na.rm = T); 
  names(overTime)[which(names(overTime) == 'Group.1')] <- xData
  names(overTime)[which(names(overTime) == 'x')] <- yData
  justforttest <- data[, c(xData,yData)]; fit <- rstatix::wilcox_test(justforttest, formula = as.formula(paste(colnames(justforttest)[2], "~", paste(colnames(justforttest)[1]), sep = "") ))
  pValue <- fit$p;   annotationInfo <- paste0("P = ", round(pValue, 2))
  if(position == "left")  { my_grob = grobTree(textGrob(annotationInfo, x=0.03,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28)))   }
  if(position == "right")  { my_grob = grobTree(textGrob(annotationInfo, x=0.4,  y=0.93, hjust=0, gp=gpar(col="black", fontsize=28)))   }
  if(position == "none") { my_grob = grobTree(textGrob(annotationInfo, x=10,  y=10, hjust=0, gp=gpar(col="black", fontsize=1)))   }
  return(
    ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam, width=0.8)) + scale_fill_manual(values = c( "#d9eafb", "#ff9a6a"))  + 
      geom_bar(data=overTime, aes_string(x=xData, y=yData), position = position_dodge(), stat = "identity",color="black",size=0.1) + 
      geom_point(aes(shape = TreatmentFactored),size=4, fill="black", color="black", alpha=0.5, position = position_jitter(width=0.3, seed=58)) +   #48
      ggtitle(title) + ylab(yLabel) +  theme_bw() + scale_shape_manual(values=c(0:2,5:6)) + labs(col="Cohort",shape="Treatment") + 
      theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=24,hjust = 0.5), axis.title.x = element_blank(), 
            plot.title = element_text(size=28,hjust = 0.5),axis.text.x=element_text(angle=45,vjust=1,hjust=1)) + annotation_custom(my_grob) # + scale_y_continuous(breaks = seq(0,100,0.2),limits=c(0,1))
    )
}

univScatter <- function(data, xData, yData, fillParam, title, xLabel, yLabel, position = "left")
{
  pearson <- round(cor(data[,xData], data[,yData], method = "pearson", use = "complete.obs"), 2)
  pValue <- cor.test(data[,xData], data[,yData], method="pearson")
  if (pValue$p.value < 0.01)   {    annotationInfo <- paste0("Pearson r = ", pearson,"\n","P = ", formatC(pValue$p.value, format="e", digits=1))    }
  if (pValue$p.value >= 0.01) {    annotationInfo <- paste0("Pearson r = ", pearson,"\n","P = ", round(pValue$p.value,2))   }
  if(position == "left")  { my_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))   }
  if(position == "right")  { my_grob = grobTree(textGrob(annotationInfo, x=0.45,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))   }
  if(position == "none") { my_grob = grobTree(textGrob(annotationInfo, x=10,  y=10, hjust=0, gp=gpar(col="black", fontsize=1)))   }
 # my_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=28)))
  return (
    ggplot(data ) + 
      geom_smooth(aes_string(x=xData, y=yData, color=fillParam, fill=fillParam), method='lm') + 
      geom_point(aes_string(x=xData, y=yData, fill=fillParam), size=8, color="black", pch=21) + theme_bw() + 
      scale_color_manual(values=c("black")) + 
      scale_fill_manual(values=c("grey90")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=28,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")     )
}

bivScatter <- function(data1, data2, name1, name2, xData, yData, fillParam, title, xLabel, yLabel, statsOff = F)
{
  if (statsOff == F)
  {  
    pearson1 <- round(cor(data1[,xData], data1[,yData], method = "pearson", use = "complete.obs"), 2)
    pValue1 <- cor.test(data1[,xData], data1[,yData], method="pearson")
    pearson2 <- round(cor(data2[,xData], data2[,yData], method = "pearson", use = "complete.obs"), 2)
    pValue2 <- cor.test(data2[,xData], data2[,yData], method="pearson")
    if (pValue1$p.value < 0.01 | pValue2$p.value < 0.01)   {    
      annotationInfo1 <- paste0(name1, " Pearson r = ", pearson1,"\n","P = ", formatC(pValue1$p.value, format="e", digits=1) )
      annotationInfo2 <- paste0("\n", name2, " Pearson r = ", pearson2,"\n","P = ", formatC(pValue2$p.value, format="e", digits=1) )  }
    if (pValue1$p.value >= 0.01 & pValue2$p.value >= 0.01) {    
      annotationInfo1 <- paste0(name1, " Pearson r = ", pearson1,"\n","P = ", round(pValue1$p.value,2) )
      annotationInfo2 <- paste0("\n",name2," Pearson r = ", pearson2, "\n","P = ", round(pValue2$p.value,2))   }
    my_grob1 = grobTree(textGrob(annotationInfo1, x=0.05,  y=0.88, hjust=0, gp=gpar(col="#7FAEDB", fontsize=28)))
    my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.75, hjust=0, gp=gpar(col="#FFB18C", fontsize=28)))
    return (
      ggplot() + 
        geom_point(data = data1, aes_string(x=xData, y=yData, fill= fillParam), size=8, color="black", pch=21) + theme_bw() + 
        geom_point(data = data2, aes_string(x=xData, y=yData, fill=fillParam), size=8, color="black", pch=21) + theme_bw() + 
        geom_smooth(data = data1, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
        geom_smooth(data = data2, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
        scale_color_manual(values=c("#FFB18C", "#7FAEDB")) + 
        scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=28,hjust = 0.5,color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5)) + 
        annotation_custom(my_grob1) + annotation_custom(my_grob2) + 
        theme(legend.position = "none")     )
  }
  if (statsOff == T)
  {  
    return (
      ggplot() + 
        geom_point(data = data1, aes_string(x=xData, y=yData, fill= fillParam), size=8, color="black", pch=21) + theme_bw() + 
        geom_point(data = data2, aes_string(x=xData, y=yData, fill=fillParam), size=8, color="black", pch=21) + theme_bw() + 
        geom_smooth(data = data1, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
        geom_smooth(data = data2, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
        scale_color_manual(values=c("#FFB18C", "#7FAEDB")) + 
        scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=28,hjust = 0.5,color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5)) + 
        #annotation_custom(my_grob1) + annotation_custom(my_grob2) + 
        theme(legend.position = "none")     )
  }
}


prePostTime <- function(data, xData, yData, fillParam, groupby, title, xLabel, yLabel)
{
  targets <- which(table(data$Subject) > 0)         # change to >1 if you want to show only the repeated measures values
  subsetData <- data[ which( data$Subject %in% names(targets)   ), ]
  subsetData <- subsetData[order(subsetData$Subject, subsetData$TimePoint, decreasing = F),]
  return(
    ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
      geom_bar(stat = "summary", aes_string(fill=fillParam),color="black",size=0.1) + 
      geom_path(aes_string(group=groupby), color="grey60") + 
      geom_point(size = 1, color="grey60") + facet_wrap(fillParam ) +   # , scales='free'
      scale_color_manual(values=c("#7FAEDB", "#FFB18C")) + 
      scale_fill_manual(values=c("#7FAEDB", "#FFB18C")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=18,hjust = 0.5, color="black"), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
            legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white"), 
            axis.title.x = element_blank(), axis.text.x = element_text(angle=45,hjust=1,vjust=1)) 
  )
}

prePostTimeAveraged <- function(data, title, xLabel, yLabel)
{
  subsetData <- data
  subsetData <- subsetData[ which( !is.na(subsetData$value) ),]
  overTime <- aggregate( subsetData$value, by= list(subsetData$variable, subsetData$TimeCategory, subsetData$Cohort), FUN=mean, na.rm = T); names(overTime)[4] <- "Mean"
  overTimeSD <- aggregate( subsetData$value, by= list(subsetData$variable, subsetData$TimeCategory, subsetData$Cohort), FUN=sd, na.rm = T); names(overTimeSD)[4] <- "SD"
  temp <- merge(overTime, overTimeSD, by = c("Group.1","Group.2", "Group.3"))
  #overTimeN <- t(as.matrix(table(subsetData$Cohort, subsetData$TimeCategory)))
  #overTimeSD$SE <- overTimeSD$x / sqrt(as.vector(overTimeN[1:nrow(overTime)]))              
  #overTime$SE <- overTimeSD$SE
  temp2 <- as.data.frame(unlist(t(as.matrix(table(subsetData$Cohort, subsetData$TimeCategory)))));  names(temp2)[3] <- "N"
  temp3 <- merge(temp, temp2, by.x = c("Group.2","Group.3"), by.y = c("Var1","Var2"))
  overTime <- temp3
  overTime$SE <- overTime$SD / sqrt(overTime$N)
  # temp <- colnames(overTime); temp[which(temp == "x")] <- "Ave"; colnames(overTime) <- temp
  annotationInfo1 <- paste0(unique(overTime$Group.3)[1])  # aPD1 group should get the peach FFB18C color
  my_grob1 = grobTree(textGrob(annotationInfo1, x=0.05,  y=0.93, hjust=0, gp=gpar(col="#FFB18C", fontsize=28)))
  annotationInfo2 <- paste0(unique(overTime$Group.3)[2])  # healthy group should get the blue 7FAEDB color
  my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.83, hjust=0, gp=gpar(col="#7FAEDB", fontsize=28)))
  
  return(
    ggplot(data=overTime, aes(x=Group.2, y=Mean, group = Group.3, color=Group.3)) + theme_bw() +   # group.1 = variable, Group.2 = TimeCategory, Group.3 = Cohort
      geom_ribbon(aes(x=Group.2, ymin=Mean-SE, ymax=Mean+SE, fill=Group.3), alpha=0.1) + 
      geom_line(aes(size=3)) +  
      geom_point(aes(size=4, fill=Group.3), pch=21) + 
      scale_color_manual(values=c("#7FAEDB", "#FFB18C")) +     
      scale_fill_manual(values=c("#7FAEDB", "#FFB18C")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=24,hjust = 0.5, color="black"), axis.title = element_text(size=28,hjust = 0.5), plot.title = element_text(size=36,hjust = 0.5), 
            axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none")  + annotation_custom(my_grob1) + annotation_custom(my_grob2)
  )
}



prePostTimeGene <- function(data, xData, yData, fillParam, groupby, title, xLabel, yLabel, paired=F)
{
  if( paired == T)
  {
    targets <- which(table(data$Subject) > 1)    # strip down to just subjects with more than one value (e.g. paired measurements)
    subsetData <- data[ which( data$Subject %in% names(targets)   ), ]
    subsetData <- subsetData[order(subsetData$Subject, subsetData$TimeCategory, decreasing = F),]
    return(
      ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
        geom_bar(stat = "summary", aes_string(fill=fillParam), color="black") + 
        geom_path(aes_string(group=groupby), color="grey60") + 
        geom_point(size = 1, color="grey60") + facet_wrap(fillParam ) + 
        scale_color_manual(values=c("#FFB18C", "#7FAEDB")) + 
        scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
              legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) 
    )
  }
  if( paired == F)   # use all data, not just the stripped down list
  {
    return(
      ggplot(data=data, aes_string(x=xData, y=yData) ) + theme_bw() + 
        geom_bar(stat = "summary", aes_string(fill=fillParam), color="black") + 
        geom_path(aes_string(group=groupby), color="grey60") + 
        geom_point(size = 1, color="grey60") + facet_wrap(fillParam ) + 
        scale_color_manual(values=c("#FFB18C", "#7FAEDB")) + 
        scale_fill_manual(values=c("#FFB18C", "#7FAEDB")) + 
        ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
        theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
              legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) 
    )
  }
}


prePostTimeAveragedGene <- function(data, title, xLabel, yLabel)
{
  subsetData <- data
  overTime <- aggregate( subsetData$value, by= list(subsetData$TimeCategory, subsetData$Cohort), FUN=mean, na.rm = T)
  overTimeSD <- aggregate( subsetData$value, by= list(subsetData$TimeCategory, subsetData$Cohort), FUN=sd, na.rm = T)
  overTimeN <- t(as.matrix(table(subsetData$Cohort, subsetData$TimeCategory)))
  overTimeSD$SE <- overTimeSD$x / sqrt(as.vector(overTimeN[1:4]))
  overTime$SE <- overTimeSD$SE
  temp <- colnames(overTime); temp[which(temp == "x")] <- "Ave"; colnames(overTime) <- temp
  overTime$Group.2 <- factor(overTime$Group.2, levels = c("Healthy","aPD1"))
  overTime <- overTime[order(overTime$Group.2, decreasing = T),]
  annotationInfo1 <- "Healthy"  
  my_grob1 = grobTree(textGrob(annotationInfo1, x=0.05,  y=0.90, hjust=0, gp=gpar(col="#7FAEDB", fontsize=24)))
  annotationInfo2 <- "aPD1"   
  my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.82, hjust=0, gp=gpar(col="#FFB18C", fontsize=24)))
  return(
    ggplot(data=overTime, aes(x=Group.1, y=Ave, group = Group.2, color=as.factor(Group.2))) + theme_bw() +     # Group.1 = TimeCategory, Group.2 = Cohort
      geom_ribbon(aes(x=Group.1, ymin=Ave-SE, ymax=Ave+SE, fill=Group.2), alpha=0.1) + 
      geom_line(aes(size=3)) +   
      geom_point(aes(size=4, fill=Group.2), pch=21) + 
      scale_color_manual(values=c("#7FAEDB", "#FFB18C")) +     
      scale_fill_manual(values=c("#7FAEDB", "#FFB18C")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
            legend.position = "none")   + annotation_custom(my_grob1) + annotation_custom(my_grob2) 
  )
}



volcanoPlot <- function( diffExprData, repelThresh, title, leftLabel, rightLabel )
{
  data <- diffExprData[which(diffExprData$lfcSE<2),]     # filter out very high SE which may be outliers
  left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.03,hjust=0, gp=gpar(col="black", fontsize=18)))
  right_grob <- grobTree(textGrob(rightLabel, x=0.9,y=0.03,hjust=0, gp=gpar(col="black", fontsize=18)))
  return(
    ggplot(data, aes(x=log2FoldChange, y=-log10(padj), label=row.names(data))) +
    geom_point(alpha=0.2, size=3, colour="black") + theme_bw() +  theme(legend.position="none") + 
    geom_text_repel(data=subset(data, padj< repelThresh), aes(label=row.names(subset(data, padj<repelThresh)))) + 
    ggtitle(title) + xlab("log2 fold change") + ylab("-log10 p-adj") + 
    theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5)) +
    annotation_custom(left_grob) + annotation_custom(right_grob)
  )
}


plotGSEAlollipop <- function( mergeResults, title, leftLabel, rightLabel, sizebyFDR=F, colorRight="#FFB18C", colorLeft="#7FAEDB")
{
  mergeResults$NAME <- toTitleCase(tolower(substr(mergeResults$NAME,start =10, stop=50)))
  mergeResults$NAME <- factor(mergeResults$NAME, levels = mergeResults$NAME[order(mergeResults$NES, decreasing = T)])
  left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.03,hjust=0, gp=gpar(col="black", fontsize=12)))
  right_grob <- grobTree(textGrob(rightLabel, x=0.8,y=0.03,hjust=0, gp=gpar(col="black", fontsize=12)))
  mergeResults <- subset(mergeResults, `FDR.q.val` < 0.05)
  if(sizebyFDR == F){
    return(
      ggplot(data=mergeResults) + 
      geom_bar( data = subset(mergeResults, `NES` > 0), aes(x=NAME, y=NES) , stat="Identity", width=0.15, fill=colorRight, size=0.01, color="black") +
      geom_bar( data = subset(mergeResults, `NES` < 0), aes(x=NAME, y=NES) , stat="Identity", width=0.15, fill=colorLeft, size=0.01, color="black") +
      geom_point(aes(x=NAME, y=NES), size=6) + 
      coord_flip() + theme_bw() + ggtitle(title) + ylab("Normalized Enrichment Score") + xlab(NULL) + 
      theme(axis.title.x = element_text(size=18), axis.text = element_text(size=14), title = element_text(size=18)) +  scale_y_continuous(minor_breaks = seq(-4,4,2)) + 
      annotation_custom(left_grob) + annotation_custom(right_grob) 
    )}
  if(sizebyFDR == T){
    return(
      ggplot(data=mergeResults) + 
        geom_bar( data = subset(mergeResults, `NES` > 0), aes(x=NAME, y=NES) , stat="Identity", width=0.15, fill=colorRight, size=0.01, color="black") +
        geom_bar( data = subset(mergeResults, `NES` < 0), aes(x=NAME, y=NES) , stat="Identity", width=0.15, fill=colorLeft, size=0.01, color="black") +        
        geom_point(aes(x=NAME, y=NES, size=FDR.q.val)) + 
        coord_flip() + theme_bw() + ggtitle(title) + ylab("Normalized Enrichment Score") + xlab(NULL) + 
        theme(axis.title.x = element_text(size=18), axis.text = element_text(size=14), title = element_text(size=18), legend.title = element_text(size=12)) + 
        scale_size( range=c(8,3),name = "False\nDiscovery\nRate") +  scale_y_continuous(minor_breaks = seq(-4,4,2)) + 
        annotation_custom(left_grob) + annotation_custom(right_grob)
    )}
}


plotGSEAtraditional <- function (mergeResults, singlePathway, pathwayName, title, leftLabel, rightLabel, cohortCompare, legendUpperRight)
{
  if (cohortCompare == F)     {     colorGrob1 <- "black"; colorGrob2 <- "black"     }
  if (cohortCompare == T)   { colorGrob1 <- "#FFB18C"; colorGrob2 <- "#7FAEDB"}
  annotationInfo <- paste0("NES: ", round(mergeResults$NES[grep(pathwayName,mergeResults$NAME)],2), "\n", "FDR: ", formatC(mergeResults$FDR.q.val[grep(pathwayName,mergeResults$NAME)], format = "e", digits = 1))
  if (legendUpperRight == T)   {   main_grob = grobTree(textGrob(annotationInfo, x=0.58,  y=0.82, hjust=0, gp=gpar(col="black", fontsize=28)))     }
  if (legendUpperRight == F)   {   main_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.15, hjust=0, gp=gpar(col="black", fontsize=28)))     }
  left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.97,hjust=0, gp=gpar(col=colorGrob1, fontsize=18)))
  right_grob <- grobTree(textGrob(rightLabel, x=0.8,y=0.97,hjust=0, gp=gpar(col=colorGrob2, fontsize=18))) 
  return(
    ggplot(data=singlePathway, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES) ) + geom_line(color="black", size=1) + 
    geom_rug(sides="b", size=0.75, alpha=0.5) + theme_bw() +
    ggtitle(title) + ylab("Enrichment score") + xlab("Rank in gene list") + 
    theme(axis.text = element_text(size=24,hjust = 0.5), axis.title = element_text(size=24,hjust = 0.5), plot.title = element_text(size=32,hjust = 0.5))+
    annotation_custom(main_grob) + geom_hline(yintercept = 0) + 
      annotation_custom(left_grob) + annotation_custom(right_grob)
  )

}



plotIPAlollipop <- function( pathways, title, leftLabel, rightLabel)
{
  pathways$NAME <- factor(pathways[,1], levels = pathways[order(pathways$z.score, decreasing = T), 1])  # factor names by zscore
  left_grob <- grobTree(textGrob(leftLabel, x=0.05,y=0.03,hjust=0, gp=gpar(col="black", fontsize=12)))
  right_grob <- grobTree(textGrob(rightLabel, x=0.8,y=0.03,hjust=0, gp=gpar(col="black", fontsize=12)))
  return(
    ggplot(data=pathways) + geom_point(aes(x=NAME, y=z.score), size=6) + 
      geom_bar( data = subset(pathways, `z.score` > 0), aes(x=NAME, y=z.score) , stat="Identity", width=0.01, color="#FFB18C", size=1) +
      geom_bar( data = subset(pathways, `z.score` < 0), aes(x=NAME, y=z.score) , stat="Identity", width=0.01, color="#7FAEDB", size=1) +
      coord_flip() + theme_bw() + ggtitle(title) + ylab("z score") + xlab(NULL) + 
      theme(axis.title.x = element_text(size=18), axis.text = element_text(size=14), title = element_text(size=18)) + 
      annotation_custom(left_grob) + annotation_custom(right_grob)
  )
}


