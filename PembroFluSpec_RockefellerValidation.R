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
  queryName="Rockefeller_freqParent", 
  viewName="demo_freqParent", 
  colFilter=NULL, 
  containerFilter=NULL
)

dataset <- labkey.data
dataset$`Visit` <- factor(dataset$`Visit`, levels=c("1","2","3"))


FC_Tfhresponse <- dcast(dataset, `Participant ID`+`Trx group`~`Visit`, value.var = c("ICOS+CD38+ c Tfh freq")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`2`/FC_Tfhresponse$`1`
#  exclude subject 87 because day 7 timepoint had unreasonably few CD4 recorded to be confident of the point estimate
FC_Tfhresponse <- FC_Tfhresponse[-which(FC_Tfhresponse$`Participant ID` == 86),]
FC_Tfhresponse$`Trx group`<- factor(FC_Tfhresponse$`Trx group`, levels=c("Chemotherapy","PD-L1 only", "TKI only", "PD1 only", "PD1 + TKI"))
fit <- lm(FC_Tfhresponse$FC~FC_Tfhresponse$`Trx group`); summary(fit)
# summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
# annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=FC_Tfhresponse, aes(x=`Trx group`, y=FC, fill=`Trx group`)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black") + 
  ggtitle("Rockefeller: ICOS+CD38+ cTfh responses") + ylab("ICOS+CD38+ cTfh fold-change") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  # annotation_custom(my_grob1) + 
  theme(legend.position = "none")
  

FC_Tfhresponse <- dcast(dataset, `Participant ID`+`Any A PD1`~`Visit`, value.var = c("ICOS+CD38+ c Tfh freq")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`2`/FC_Tfhresponse$`1`
#  exclude subject 87 because day 7 timepoint had unreasonably few CD4 recorded to be confident of the point estimate
FC_Tfhresponse <- FC_Tfhresponse[-which(FC_Tfhresponse$`Participant ID` == 86),]
FC_Tfhresponse$`Any A PD1`<- factor(FC_Tfhresponse$`Any A PD1`, levels=c("F","T"))
fit <- lm(FC_Tfhresponse$FC~FC_Tfhresponse$`Any A PD1`); summary(fit)
# summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
# annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=FC_Tfhresponse, aes(x=`Any A PD1`, y=FC, fill=`Any A PD1`)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black") + 
  ggtitle("Rockefeller: ICOS+CD38+ cTfh responses") + ylab("ICOS+CD38+ cTfh fold-change") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), plot.title = element_text(size=28,hjust = 0.5)) + 
  # annotation_custom(my_grob1) + 
  theme(legend.position = "none")


FC_Tfhresponse <- dcast(dataset, `Participant ID`+`Treatment`~`Visit`, value.var = c("ICOS+CD38+ c Tfh freq")); 
FC_Tfhresponse$FC <- FC_Tfhresponse$`2`/FC_Tfhresponse$`1`
#  exclude subject 87 because day 7 timepoint had unreasonably few CD4 recorded to be confident of the point estimate
FC_Tfhresponse <- FC_Tfhresponse[-which(FC_Tfhresponse$`Participant ID` == 86),]
FC_Tfhresponse$`Treatment`<- factor(FC_Tfhresponse$`Treatment`, levels=c("Carbo Gem","axitinib","cabozantinib","pazopanib",
                                                                         "Atezolizumab","Pembrolizumab","Nivolumab",
                                                                         "Pembrolizumab/lenvatinib","Nivolumab/cabozantinib"))
fit <- lm(FC_Tfhresponse$FC~FC_Tfhresponse$`Treatment`); summary(fit)
# summary(fit); pValue <- round(summary(fit)$coefficients[2,"Pr(>|t|)"],2);   CI <- confint(fit)[2,]; CI <- round(CI,2)
# annotationInfo <- paste0("p = ", pValue,"\n","95%CI: [",CI[1],",",CI[2],"]"); my_grob1 = grobTree(textGrob(annotationInfo, x=0.63,  y=0.92, hjust=0, gp=gpar(col="black", fontsize=18)))
ggplot(data=FC_Tfhresponse, aes(x=`Treatment`, y=FC, fill=`Treatment`)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(size=5, pch=21, width=0.05, fill="black") + 
  ggtitle("Rockefeller: ICOS+CD38+ cTfh responses") + ylab("ICOS+CD38+ cTfh fold-change") + # scale_y_continuous(breaks=seq(1:3), limits=c(0,3)) +   
  theme(axis.text = element_text(size=22,hjust = 0.5), axis.title = element_text(size=22,hjust = 0.5), axis.title.x = element_blank(), 
        plot.title = element_text(size=28,hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1)) + 
  # annotation_custom(my_grob1) + 
  theme(legend.position = "none") 





