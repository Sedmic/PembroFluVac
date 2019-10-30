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
# library("shadowtext")


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
  my_grob = grobTree(textGrob(annotationInfo, x=0.07,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
  return (
    ggplot(data=data, aes_string(x=xData, y=yData, fill=fillParam)) + scale_fill_manual(values = c("#abcdef", "#ffcab1")) +  
      geom_boxplot(outlier.shape = NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black", color="white", stroke=1) + theme_bw() + 
      ggtitle(title) + ylab(yLabel) +  
      theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")     )
}


univScatter <- function(data, xData, yData, fillParam, title, xLabel, yLabel)
{
  pearson <- round(cor(data[,xData], data[,yData], method = "pearson", use = "complete.obs"), 2)
  pValue <- cor.test(data[,xData], data[,yData], method="pearson")
  if (pValue$p.value < 0.01)   {    annotationInfo <- paste0("Pearson r = ", pearson,"\n","P = ", formatC(pValue$p.value, format="e", digits=1))    }
  if (pValue$p.value >= 0.01) {    annotationInfo <- paste0("Pearson r = ", pearson,"\n","P = ", round(pValue$p.value,2))   }
  my_grob = grobTree(textGrob(annotationInfo, x=0.05,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
  return (
    ggplot(data ) + 
      geom_smooth(aes_string(x=xData, y=yData, color=fillParam, fill=fillParam), method='lm') + 
      geom_point(aes_string(x=xData, y=yData, fill=fillParam), size=7, color="black", pch=21) + theme_bw() + 
      scale_color_manual(values=c("black")) + 
      scale_fill_manual(values=c("grey90")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5)) + 
      annotation_custom(my_grob) + theme(legend.position = "none")     )
}

bivScatter <- function(data1, data2, name1, name2, xData, yData, fillParam, title, xLabel, yLabel)
{
  pearson1 <- round(cor(data1[,xData], data1[,yData], method = "pearson", use = "complete.obs"), 2)
  pValue1 <- cor.test(data1[,xData], data1[,yData], method="pearson")
  pearson2 <- round(cor(data2[,xData], data2[,yData], method = "pearson", use = "complete.obs"), 2)
  pValue2 <- cor.test(data2[,xData], data2[,yData], method="pearson")
  if (pValue1$p.value < 0.01 | pValue2$p.value < 0.01)   {    
    annotationInfo1 <- paste0(name1, " Pearson r = ", pearson1,"\n","P = ", formatC(pValue1$p.value, format="e", digits=1) )
    annotationInfo2 <- paste0("\n", name2, " Pearson r = ", pearson2,"\n","P = ", formatC(pValue2$p.value, format="e", digits=1) )  }
  if (pValue2$p.value >= 0.01 | pValue2$p.value >= 0.01) {    
    annotationInfo1 <- paste0("Pearson r = ", pearson1,"\n","P = ", round(pValue1$p.value,2) )
    annotationInfo2 <- paste0("\n","Pearson r = ", pearson2, "\n","P fr f = ", round(pValue2$p.value,2))   }
  my_grob1 = grobTree(textGrob(annotationInfo1, x=0.05,  y=0.85, hjust=0, gp=gpar(col="#abcdef", fontsize=24)))
  # my_grob2 = shadowtext::shadowtextGrob(label=annotationInfo2, x=0.15, y=0.75, bg.r=0.01,  gp=gpar(col="#ffcab1",cex=2))
  my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.75, hjust=0, gp=gpar(col="#ffcab1", fontsize=24)))
  return (
    ggplot() + 
      geom_point(data = data1, aes_string(x=xData, y=yData, fill= fillParam), size=7, color="black", pch=21) + theme_bw() + 
      geom_point(data = data2, aes_string(x=xData, y=yData, fill=fillParam), size=7, color="black", pch=21) + theme_bw() + 
      geom_smooth(data = data1, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
      geom_smooth(data = data2, aes_string(x=xData, y=yData, color=fillParam, fill=NA), method='lm', fullrange=T) +
      scale_color_manual(values=c("#abcdef", "#ffcab1")) + 
      scale_fill_manual(values=c("#abcdef", "#ffcab1")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5)) + 
      annotation_custom(my_grob1) + annotation_custom(my_grob2) + theme(legend.position = "none")     )
}


prePostTime <- function(data, xData, yData, fillParam, groupby, title, xLabel, yLabel)
{
  targets <- which(table(data$Subject) > 1)
  subsetData <- data[ which( data$Subject %in% names(targets)   ), ]
  subsetData <- subsetData[order(subsetData$Subject, subsetData$TimePoint, decreasing = F),]
  return(
    ggplot(data=subsetData, aes_string(x=xData, y=yData) ) + theme_bw() + 
      geom_bar(stat = "summary", aes_string(fill=fillParam), color="black") + 
      geom_path(aes_string(group=groupby), color="grey60") + 
      geom_point(size = 1, color="grey60") + facet_wrap(fillParam ) + 
      scale_color_manual(values=c("#abcdef", "#ffcab1")) + 
      scale_fill_manual(values=c("#abcdef", "#ffcab1")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
            legend.position = "none", strip.text = element_text(size = 24, color="black"), strip.background = element_rect(fill="white")) 
      
  )
}

prePostTimeAveraged <- function(data, title, xLabel, yLabel)
{
  subsetData <- data
  overTime <- aggregate( subsetData$value, by= list(subsetData$variable, subsetData$TimeCategory, subsetData$Cohort), FUN=mean, na.rm = T)
  overTimeSD <- aggregate( subsetData$value, by= list(subsetData$variable, subsetData$TimeCategory, subsetData$Cohort), FUN=sd, na.rm = T)
  overTimeN <- t(as.matrix(table(subsetData$Cohort, subsetData$TimeCategory)))
  overTimeSD$SE <- overTimeSD$x / sqrt(as.vector(overTimeN[1:6]))
  overTime$SE <- overTimeSD$SE
  temp <- colnames(overTime); temp[which(temp == "x")] <- "Ave"; colnames(overTime) <- temp
  annotationInfo1 <- paste0(unique(overTime$Group.3)[1])   
  my_grob1 = grobTree(textGrob(annotationInfo1, x=0.05,  y=0.90, hjust=0, gp=gpar(col="#abcdef", fontsize=24)))
  annotationInfo2 <- paste0(unique(overTime$Group.3)[2])   
  my_grob2 = grobTree(textGrob(annotationInfo2, x=0.05,  y=0.82, hjust=0, gp=gpar(col="#ffcab1", fontsize=24)))
  
  return(
    ggplot(data=overTime, aes(x=Group.2, y=Ave, group = Group.3, color=as.factor(Group.3))) + theme_bw() + 
      geom_ribbon(aes(x=Group.2, ymin=Ave-SE, ymax=Ave+SE, fill=Group.3), alpha=0.1) + 
      geom_line(aes(size=3)) +  
      geom_point(aes(size=4, fill=Group.3), pch=21) + 
      scale_color_manual(values=c("#abcdef", "#ffcab1")) +     
      scale_fill_manual(values=c("#abcdef", "#ffcab1")) + 
      ggtitle(title) + ylab(yLabel) + xlab(xLabel)  +
      theme(axis.text = element_text(size=18,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), plot.title = element_text(size=28,hjust = 0.5), 
            legend.position = "none")  + annotation_custom(my_grob1) + annotation_custom(my_grob2)
  )

}



















