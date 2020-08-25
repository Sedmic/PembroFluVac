library(cytofkit2)
library(flowCore)
library(flowStats)
library(uwot)
library(FlowSOM)
library(flowCore)
library(Biobase)
library(ggplot2)
library(tidyverse)
library(mem)
library(uwot)
library(RColorBrewer)
library(reshape2)
library(tibble)
library(pheatmap)




# with help from https://github.com/cytolab/roe_phospho_flow/blob/master/01_data_analysis.Rmd


COFACTOR_FOR_ARCSINH_SCALE = 5          # set cofactor for all channels
NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE = 1000     # choose number of cells to sample per file 
# (if you want proportional sampling, enter a number larger than the max number of cells in a single file)
TARGET_NUMBER_OF_CLUSTERS = 10          # choose target number of clusters for FlowSOM
SEED = 1                                # set seed for reproducible results (can choose any number)
a <- getwd()
setwd(dir = "D:/Pembro-Fluvac/18-19season/Flow cytometry/Tcell/exportNonNavCD4/")  
fcs.files <- list.files(path = '.', pattern='.fcs$', full=TRUE, ignore.case=T)
data.lists <- lapply(lapply(fcs.files, read.FCS), exprs)
my.data = as.data.frame(do.call(rbind, mapply(cbind, data.lists, "File_ID" = c(1:length(data.lists)), SIMPLIFY = F)))
orig_names = c(colnames(my.data),"UMAP1","UMAP2","cluster")
colnames(my.data)[1:length(my.data) - 1] <- as.character(read.FCS(fcs.files[[1]])@parameters@data[["desc"]])

# equally sample the data 
files.to.sample = split(my.data,my.data$`File_ID`)
sampled.data <- list()
for (i in 1: length(files.to.sample)){
  if (nrow(files.to.sample[[i]])>NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE){            
    sample.df =  files.to.sample[[i]]
    sampled.data[[i]] = as.data.frame(sample.df[sample(nrow(sample.df), NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE), ])}        
  else{
    sampled.data[[i]] = files.to.sample[[i]]}}
my.sampled.data = as.data.frame(do.call(rbind, sampled.data))   

# Run UMAP on chosen markers
# select all channels to use in UMAP by OPENING CONSOLE below
set.seed(SEED)
chosen.markers = as.data.frame(as.data.frame(my.sampled.data)[,c(choose.markers(my.sampled.data))])        # ************** INTERACTIVE STEP ********
umap.markers <- chosen.markers %>%
  mutate_all(function(x)
    asinh(x / COFACTOR_FOR_ARCSINH_SCALE))
myumap <- umap(umap.markers, ret_model = TRUE,verbose = TRUE)
umap.data = as.data.frame(myumap$embedding)
colnames(umap.data) <- c("UMAP1", "UMAP2")

# create flowFrame from UMAP data
umap.matrix <- as.matrix(umap.data)
UMAP.metadata <- data.frame(name = dimnames(umap.matrix)[[2]], desc = paste('UMAP', dimnames(umap.matrix)[[2]]))
UMAP.metadata$range <- apply(apply(umap.matrix, 2, range), 2, diff)
UMAP.metadata$minRange <- apply(umap.matrix, 2, min)
UMAP.metadata$maxRange <- apply(umap.matrix, 2, max)
umap.flowframe <- new("flowFrame", exprs=umap.matrix,parameters = AnnotatedDataFrame(UMAP.metadata))
# run FlowSOM on UMAP axes
fSOM.umap <- FlowSOM(umap.flowframe, compensate = FALSE, transform = FALSE, toTransform=c(1:2), scale = TRUE, colsToUse = c(1:2), nClus = TARGET_NUMBER_OF_CLUSTERS, seed = SEED)
FlowSOM.clusters <- as.numeric(as.vector(as.matrix(fSOM.umap[[2]][fSOM.umap[[1]]$map$mapping[,1]])))
analysis.data = as.data.frame(cbind(my.sampled.data,umap.data,FlowSOM.clusters))
colnames(analysis.data)[ncol(analysis.data)]<-"cluster"

# plot UMAP axes and FlowSOM per file used in analysis as well as for the concatenated data (last plot)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = col_vector[-c(4,17,19,27,29:45)]  # take away certain olors
values = sample(col_vector)
range <- apply(apply(umap.data, 2, range), 2, diff)
graphical.ratio <- (range[1]/range[2])

# plot FlowSOM clusters on UMAP axes for concatenated data
ggplot(data.frame(x = analysis.data$UMAP1, y = analysis.data$UMAP2, col = as.factor(analysis.data$cluster))) + coord_fixed(ratio=graphical.ratio) + 
  geom_point(aes(x=x, y=y, color=col),cex = 1.5) + guides(colour = guide_legend(override.aes = list(size=5), nrow = 13)) +
  labs(x = "UMAP 1", y = "UMAP 2",title = "FlowSOM Clustering on UMAP Axes (concatenated)", color = "FlowSOM Cluster") + theme_bw() + scale_color_manual(values = values)

# plot FlowSOM clusters on UMAP axes for each file used in analysis
separate.fcs.files = split(analysis.data,analysis.data$`File_ID`)
for (i in 1:length(separate.fcs.files)){
  newname  = str_remove(fcs.files[i], ".fcs")
  plot <- data.frame(x = separate.fcs.files[[i]][["UMAP1"]], y = separate.fcs.files[[i]][["UMAP2"]], col = as.factor(separate.fcs.files[[i]][["cluster"]]))
  print(ggplot(plot) + geom_point(aes(x=x, y=y, col = col)) +
          coord_fixed(ratio=graphical.ratio)+ labs(color = "FlowSOM Cluster", x = "UMAP1", y = "UMAP2", title = "FlowSOM Clustering on UMAP Axes",caption = newname) + 
          scale_color_manual(values = values) + guides(colour = guide_legend(override.aes = list(size=5)))+ theme_bw() + theme(plot.caption = element_text(size = 6)))}


# export new FCS files with UMAP1, UMAP2, and FlowSOM cluster ID as added channels 
# markers used for the UMAP will also have a (u) next to their channel name
my.new.data  = analysis.data
for (s in 1:ncol(my.new.data)){
  colnames(my.new.data)[which(colnames(my.new.data) == colnames(umap.markers)[s])] <- paste0(colnames(analysis.data)[which(colnames(analysis.data) == colnames(umap.markers)[s])]," (u)")}
new_desc = colnames(my.new.data)
colnames(analysis.data)<-orig_names
separate.files = split(analysis.data,as.factor(analysis.data$`File_ID`))
for (i in 1:length(separate.files)){
  single.file = separate.files[[i]]
  remove.ID  = single.file[-c(ncol(my.sampled.data))]
  mat <- as.matrix(single.file)
  metadata <- data.frame(name = dimnames(mat)[[2]], desc = new_desc)
  metadata$range <- apply(apply(mat, 2, range), 2, diff)
  metadata$minRange <- apply(mat, 2, min)
  metadata$maxRange <- apply(mat, 2, max)
  export.flowframe <- new("flowFrame", exprs = mat, parameters = AnnotatedDataFrame(metadata))
  newname  = str_remove(fcs.files[i], ".fcs")
  filename = paste0(newname,"_UMAP_FlowSOM.fcs")
  write.FCS(export.flowframe,filename = filename)
  print(i)80 
  }







COFACTOR_FOR_ARCSINH_SCALE = 5          # set cofactor for all channels
NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE = 10000     # choose number of cells to sample per file 
# (if you want proportional sampling, enter a number larger than the max number of cells in a single file)
TARGET_NUMBER_OF_CLUSTERS = 20          # choose target number of clusters for FlowSOM
SEED = 1                                # set seed for reproducible results (can choose any number)
origDir <- getwd()
setwd(dir = "D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/exportNonnavB/d0")
fcs.files <- list.files(path = '.', pattern='.fcs$', full=TRUE, ignore.case=T)
data.lists <- lapply(lapply(fcs.files, read.FCS), exprs)                          # expression data, one list element per file
data.lists <- lapply(data.lists, function(df) tibble::rownames_to_column(as.data.frame(df),var="cellID"))
my.data = as.data.frame(do.call(rbind, mapply(cbind, data.lists, "File_ID" = c(1:length(data.lists)), SIMPLIFY = F)))
orig_names = c(colnames(my.data),"UMAP1","UMAP2","cluster")
colnames(my.data)[1:length(my.data) - 1] <- as.character(read.FCS(fcs.files[[1]])@parameters@data[["desc"]])

# equally sample the data 
files.to.sample = split(my.data,my.data$`File_ID`)
sampled.data <- list()
for (i in 1: length(files.to.sample)){
  if (nrow(files.to.sample[[i]])>NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE){            
    sample.df =  files.to.sample[[i]]
    sampled.data[[i]] = as.data.frame(sample.df[sample(nrow(sample.df), NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE), ])}        
  else{
    sampled.data[[i]] = files.to.sample[[i]]}}
my.sampled.data = as.data.frame(do.call(rbind, sampled.data))   

# Run UMAP on chosen markers
# select all channels to use in UMAP by OPENING CONSOLE below
set.seed(SEED)
chosen.markers = as.data.frame(as.data.frame(my.sampled.data)[,c(choose.markers(my.sampled.data))])        # ************** INTERACTIVE STEP *****************
#  chosen markers:  CD27, Bcl6, IgD, CD32, IgM, CD86, CD138, CD38, Blimp1, CD11c, CD71, CXCR4, CD21, CD14, CXCR5, CD16, CD23, CD80, Tbet, Ki67, CD20         7:11,15:32
umap.markers <- chosen.markers %>%   mutate_all(function(x)    asinh(x / COFACTOR_FOR_ARCSINH_SCALE))
myumap <- umap(umap.markers, ret_model = TRUE,verbose = TRUE)
umap.data = as.data.frame(myumap$embedding)
colnames(umap.data) <- c("UMAP1", "UMAP2")

# create flowFrame from UMAP data
umap.matrix <- as.matrix(umap.data)
UMAP.metadata <- data.frame(name = dimnames(umap.matrix)[[2]], desc = paste('UMAP', dimnames(umap.matrix)[[2]]))
UMAP.metadata$range <- apply(apply(umap.matrix, 2, range), 2, diff)
UMAP.metadata$minRange <- apply(umap.matrix, 2, min)
UMAP.metadata$maxRange <- apply(umap.matrix, 2, max)
umap.flowframe <- new("flowFrame", exprs=umap.matrix,parameters = AnnotatedDataFrame(UMAP.metadata))
# run FlowSOM on UMAP axes
fSOM.umap <- FlowSOM(umap.flowframe, compensate = FALSE, transform = FALSE, toTransform=c(1:2), scale = TRUE, colsToUse = c(1:2), nClus = TARGET_NUMBER_OF_CLUSTERS, seed = SEED)
FlowSOM.clusters <- as.numeric(as.vector(as.matrix(fSOM.umap[[2]][fSOM.umap[[1]]$map$mapping[,1]])))
analysis.data = as.data.frame(cbind(my.sampled.data,umap.data,FlowSOM.clusters))
colnames(analysis.data)[ncol(analysis.data)]<-"cluster"

# plot UMAP axes and FlowSOM per file used in analysis as well as for the concatenated data (last plot)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = col_vector[-c(4,17,19,27,29:45)]  # take away certain olors
values = sample(col_vector)
range <- apply(apply(umap.data, 2, range), 2, diff)
graphical.ratio <- (range[1]/range[2])

# plot FlowSOM clusters on UMAP axes for concatenated data                                          # pch=21, color='black',stroke=0.05,
subsampled <- data.frame(UMAP1 = analysis.data$UMAP1, UMAP2 = analysis.data$UMAP2, cluster = as.factor(analysis.data$cluster))
subsampled <- subsampled[sample(1:nrow(subsampled), 5000),]
# ggplot(data.frame(x = analysis.data$UMAP1, y = analysis.data$UMAP2, col = as.factor(analysis.data$cluster))) + coord_fixed(ratio=graphical.ratio) + 
ggplot(subsampled, aes(x = UMAP1, y = UMAP2, col = cluster)) + coord_fixed(ratio=graphical.ratio) + 
  geom_point(size = 3) + guides(colour = guide_legend(override.aes = list(size=5), nrow = 13)) + 
  labs(x = "UMAP 1", y = "UMAP 2",title = "FlowSOM Clustering for Nonnaive B cells at baseline", color = "FlowSOM Cluster") + theme_bw() + scale_color_manual(values = values)
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_UMAP_nonnavB_d0_10kcells_20clusters.pdf")

# separateUMAPplots(analysis.data,fcs.files)


propClusters <- analysis.data[,-c(1:6,14,33:34,36:37)] %>% group_by(File_ID, cluster)  %>% summarize(n=n())  %>% mutate(freq = n/sum(n))
propClusters$n <- NULL
propClusters$File_ID <- sapply(propClusters$File_ID, function(x) { fcs.files[x] } )
propClusters$Cohort <- "HC"; propClusters$Cohort[grep("19616",propClusters$File_ID, value = F)] <- "aPD1"
propClusters$cluster <- factor(propClusters$cluster)
ggplot(propClusters, aes(x=cluster, y=freq, group=Cohort, color=Cohort)) +  stat_smooth(method="loess",span=0.1,se=TRUE, aes(fill=Cohort), alpha=0.1) + theme_bw() + 
  ggtitle("Nonnaive B cells at baseline") + ylab("Frequency (% nonnaive B") + xlab("FlowSOM cluster") + 
  scale_fill_manual(values = c("#FFB18C", "#7FAEDB")) + scale_color_manual(values = c("#FFB18C", "#7FAEDB" )) + 
  theme(axis.text.x = element_text(angle=45,vjust = 0.9, hjust=1), axis.text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=24),
        legend.text = element_text(size=12), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_clusterComparison_byCohort_nonnavB_d0.pdf")
clusterSummary <- analysis.data[,-c(1:6,14,33:37)] %>% group_by(cluster)  %>% # mutate_all(function(x)  asinh(x / COFACTOR_FOR_ARCSINH_SCALE))%>% 
  summarize_all(list(mean)) 
clusterSummary <- clusterSummary %>% column_to_rownames("cluster")
clusterSummary <- clusterSummary[,-grep(paste0(c("CD4","CD3"), collapse='|'), colnames(clusterSummary),value=F)]
draw_colnames_65 <- function (coln, gaps, ...) {              # https://slowkow.com/notes/pheatmap-tutorial/
  coord <- pheatmap:::find_coordinates(length(coln), gaps);  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(  coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.75, hjust = 1, rot = 90, gp = grid::gpar(...)   )
  return(res)     }
assignInNamespace(   x = "draw_colnames",   value = "draw_colnames_65",   ns = asNamespace("pheatmap") )

pheatmap(clusterSummary, scale = "column", cluster_rows = F, color=viridis::inferno(100), main="Cluster expression for all subjects", fontsize=14, 
         border_color = "grey30", cellwidth = 19, angle_col = 90,
          filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_clusterComparison_byCohort_nonnavB_d0_heatmap.pdf",
         )  
dev.off()

propClusters %>% group_by(cluster) %>% rstatix::t_test(data= ., formula =  freq ~ Cohort) %>% rstatix::adjust_pvalue(method = 'fdr')
propClusters %>% group_by(cluster==1) %>% rstatix::t_test(data= ., formula =  freq ~ Cohort)
propClusters %>% group_by(cluster==18) %>% rstatix::t_test(data= ., formula =  freq ~ Cohort)

writeFCSfiles( analysis.data, umap.markers, orig_names)


##  *********************************************************************************************************************************************************
##  *********************************************************************************************************************************************************
##  *********************************************************************************************************************************************************


COFACTOR_FOR_ARCSINH_SCALE = 5          # set cofactor for all channels
NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE = 10000     # choose number of cells to sample per file 
TARGET_NUMBER_OF_CLUSTERS = 20          # choose target number of clusters for FlowSOM
SEED = 1                                # set seed for reproducible results (can choose any number)
origDir <- getwd()
setwd(dir = "D:/Pembro-Fluvac/18-19season/Flow cytometry/Tcell/exportNonNavCD4/")
fcs.files <- list.files(path = '.', pattern='.fcs$', full=TRUE, ignore.case=T)
data.lists <- lapply(lapply(fcs.files, read.FCS), exprs)                          # expression data, one list element per file
data.lists <- lapply(data.lists, function(df) tibble::rownames_to_column(as.data.frame(df),var="cellID"))
my.data = as.data.frame(do.call(rbind, mapply(cbind, data.lists, "File_ID" = c(1:length(data.lists)), SIMPLIFY = F)))
orig_names = c(colnames(my.data),"UMAP1","UMAP2","cluster")
colnames(my.data)[1:length(my.data) - 1] <- as.character(read.FCS(fcs.files[[1]])@parameters@data[["desc"]])

# equally sample the data 
files.to.sample = split(my.data,my.data$`File_ID`)
sampled.data <- list()
for (i in 1: length(files.to.sample)){
  if (nrow(files.to.sample[[i]])>NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE){            
    sample.df =  files.to.sample[[i]]
    sampled.data[[i]] = as.data.frame(sample.df[sample(nrow(sample.df), NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE), ])}        
  else{
    sampled.data[[i]] = files.to.sample[[i]]}}
my.sampled.data = as.data.frame(do.call(rbind, sampled.data))   

# Run UMAP on chosen markers
# select all channels to use in UMAP by OPENING CONSOLE below
set.seed(SEED)
chosen.markers = as.data.frame(as.data.frame(my.sampled.data)[,c(choose.markers(my.sampled.data))])        # ************** INTERACTIVE STEP *****************
#  chosen markers:  CD27, Bcl6, IgD, CD32, IgM, CD86, CD138, CD38, Blimp1, CD11c, CD71, CXCR4, CD21, CD14, CXCR5, CD16, CD23, CD80, Tbet, Ki67, CD20         7:11,15:32
umap.markers <- chosen.markers %>%   mutate_all(function(x)    asinh(x / COFACTOR_FOR_ARCSINH_SCALE))
myumap <- umap(umap.markers, ret_model = TRUE,verbose = TRUE)
umap.data = as.data.frame(myumap$embedding)
colnames(umap.data) <- c("UMAP1", "UMAP2")

# create flowFrame from UMAP data
umap.matrix <- as.matrix(umap.data)
UMAP.metadata <- data.frame(name = dimnames(umap.matrix)[[2]], desc = paste('UMAP', dimnames(umap.matrix)[[2]]))
UMAP.metadata$range <- apply(apply(umap.matrix, 2, range), 2, diff)
UMAP.metadata$minRange <- apply(umap.matrix, 2, min)
UMAP.metadata$maxRange <- apply(umap.matrix, 2, max)
umap.flowframe <- new("flowFrame", exprs=umap.matrix,parameters = AnnotatedDataFrame(UMAP.metadata))
# run FlowSOM on UMAP axes
fSOM.umap <- FlowSOM(umap.flowframe, compensate = FALSE, transform = FALSE, toTransform=c(1:2), scale = TRUE, colsToUse = c(1:2), nClus = TARGET_NUMBER_OF_CLUSTERS, seed = SEED)
FlowSOM.clusters <- as.numeric(as.vector(as.matrix(fSOM.umap[[2]][fSOM.umap[[1]]$map$mapping[,1]])))
analysis.data = as.data.frame(cbind(my.sampled.data,umap.data,FlowSOM.clusters))
colnames(analysis.data)[ncol(analysis.data)]<-"cluster"

# plot UMAP axes and FlowSOM per file used in analysis as well as for the concatenated data (last plot)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = col_vector[-c(4,17,19,27,29:45)]  # take away certain olors
values = sample(col_vector)
range <- apply(apply(umap.data, 2, range), 2, diff)
graphical.ratio <- (range[1]/range[2])

# plot FlowSOM clusters on UMAP axes for concatenated data                                          # pch=21, color='black',stroke=0.05,
ggplot(data.frame(x = analysis.data$UMAP1, y = analysis.data$UMAP2, col = as.factor(analysis.data$cluster))) + coord_fixed(ratio=graphical.ratio) + 
  geom_point(aes(x=x, y=y, color=col),size = 3, alpha=0.75) + guides(colour = guide_legend(override.aes = list(size=5), nrow = 13)) + 
  labs(x = "UMAP 1", y = "UMAP 2",title = "FlowSOM Clustering for Nonnaive B cells at baseline", color = "FlowSOM Cluster") + theme_bw() + scale_color_manual(values = values)
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_UMAP_nonnavB_d0_10kcells_20clusters.pdf")


# plot FlowSOM clusters on UMAP axes for each file used in analysis
separate.fcs.files = split(analysis.data,analysis.data$`File_ID`)
for (i in 1:1) { # length(separate.fcs.files)){            # remove the comment sign # to apply this to all files in the set
  newname  = str_remove(fcs.files[i], ".fcs")
  plot <- data.frame(x = separate.fcs.files[[i]][["UMAP1"]], y = separate.fcs.files[[i]][["UMAP2"]], col = as.factor(separate.fcs.files[[i]][["cluster"]]))
  print(ggplot(plot) + geom_point(aes(x=x, y=y, col = col)) +
          coord_fixed(ratio=graphical.ratio)+ labs(color = "FlowSOM Cluster", x = "UMAP1", y = "UMAP2", title = "FlowSOM Clustering on UMAP Axes",caption = newname) + 
          scale_color_manual(values = values) + guides(colour = guide_legend(override.aes = list(size=5)))+ theme_bw() + theme(plot.caption = element_text(size = 6)))}


propClusters <- analysis.data[,-c(1:6,14,33:34,36:37)] %>% group_by(File_ID, cluster)  %>% summarize(n=n())  %>% mutate(freq = n/sum(n))
propClusters$n <- NULL
propClusters$File_ID <- sapply(propClusters$File_ID, function(x) { fcs.files[x] } )
propClusters$Cohort <- "HC"; propClusters$Cohort[grep("19616",propClusters$File_ID, value = F)] <- "aPD1"
propClusters$cluster <- factor(propClusters$cluster)
ggplot(propClusters, aes(x=cluster, y=freq, group=Cohort, color=Cohort)) +  stat_smooth(method="loess",span=0.1,se=TRUE, aes(fill=Cohort), alpha=0.1) + theme_bw() + 
  ggtitle("Nonnaive B cells at baseline") + ylab("Frequency (% nonnaive B") + xlab("FlowSOM cluster") + 
  theme(axis.text.x = element_text(angle=45,vjust = 0.9, hjust=1), axis.text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=24),
        legend.text = element_text(size=12), legend.title = element_text(size=14))
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_clusterComparison_byCohort_nonnavB_d0.pdf")
clusterSummary <- analysis.data[,-c(1:6,14,33:37)] %>% group_by(cluster)  %>% # mutate_all(function(x)  asinh(x / COFACTOR_FOR_ARCSINH_SCALE))%>% 
  summarize_all(list(mean)) 
clusterSummary <- clusterSummary %>% column_to_rownames("cluster")

draw_colnames_45 <- function (coln, gaps, ...) {              # https://slowkow.com/notes/pheatmap-tutorial/
  coord <- pheatmap:::find_coordinates(length(coln), gaps);  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(  coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)   )
  return(res)     }
assignInNamespace(   x = "draw_colnames",   value = "draw_colnames_45",   ns = asNamespace("pheatmap") )

pheatmap(clusterSummary, scale = "column", cluster_rows = F, color=viridis::inferno(100), main="Cluster expression for all subjects", fontsize=12, border_color = "grey30",
         #filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_clusterComparison_byCohort_nonnavB_d0_heatmap.pdf",
);  dev.off()



writeFCSfiles( analysis.data, umap.markers, orig_names)



# placeholder










# function to choose markers to use in your analysis (do not need to change anything here)
choose.markers <- function(exp_data) {
  print("Numbered column names, in order they appear in file: ")
  print(paste(c(1:(ncol(exp_data))), ": ", 
              colnames(exp_data[, c(1:(ncol(exp_data)))]), sep = ""))
  markers = readline("Enter column numbers to include (e.g. 1:5,6,8:10).\n")
  sep_vals = unlist(strsplit(markers, ","))
  list_vals = vector()
  for (i in 1:length(sep_vals)) {
    val = sep_vals[i]
    if (length(unlist(strsplit(val, ":"))) > 1) {
      new_val = as.numeric(unlist(strsplit(val, ":"))[1]):
        as.numeric(unlist(strsplit(val, ":"))[2])
    } else{
      new_val = as.numeric(sep_vals[i])
    }
    list_vals = c(list_vals, new_val)
  }
  markerList = c(list_vals)
  return(markerList)
}

writeFCSfiles <- function( analysis.data, umap.markers, orig_names)
{
  # export new FCS files with UMAP1, UMAP2, and FlowSOM cluster ID as added channels 
  # markers used for the UMAP will also have a (u) next to their channel name
  my.new.data  = analysis.data
  for (s in 1:ncol(my.new.data)){
    colnames(my.new.data)[which(colnames(my.new.data) == colnames(umap.markers)[s])] <- paste0(colnames(analysis.data)[which(colnames(analysis.data) == colnames(umap.markers)[s])]," (u)")}
  new_desc = colnames(my.new.data) #[-c(35)]
  colnames(analysis.data)<-orig_names
  separate.files = split(analysis.data,as.factor(analysis.data$`File_ID`))          # separate the analysis.data into respective file IDs
  dir.create(a <- str_replace_all( Sys.time(), pattern=":",replacement=""))
  setwd(dir = a)
  for (i in  1:length(separate.files)){  # 1:3) {#
    single.file = separate.files[[i]]                               # take list elements one at a time
    single.file = single.file[c("cellID","UMAP1","UMAP2","cluster")]
    single.file = merge(x = data.lists[[i]], y = single.file, by.x = 'cellID', by.y = 'cellID' )
    single.file$cellID <- NULL
    mat <- as.matrix(single.file)                                   # reformat data as matrix 
    metadata <- data.frame(name = dimnames(mat)[[2]], desc = new_desc[-c(34:35)])
    metadata$range <- apply(apply(mat, 2, range), 2, diff)
    metadata$minRange <- apply(mat, 2, min)
    metadata$maxRange <- apply(mat, 2, max)
    export.flowframe <- new("flowFrame", exprs = mat, parameters = AnnotatedDataFrame(metadata))          # object with all of the key data for the FCS file
    newname  = str_remove(fcs.files[i], ".fcs")
    filename = paste0(newname,"_UMAP_FlowSOM.fcs")
    write.FCS(export.flowframe,filename = filename)                 # write the final FCS file
    print(i)
  }
  setwd(origDir)
}

separateUMAPplots <- function (analysis.data, fcs.files)
{
  # plot FlowSOM clusters on UMAP axes for each file used in analysis
  separate.fcs.files = split(analysis.data,analysis.data$`File_ID`)
  for (i in 1:length(separate.fcs.files)){            # remove the comment sign # to apply this to all files in the set
    newname  = str_remove(fcs.files[i], ".fcs")
    plot <- data.frame(x = separate.fcs.files[[i]][["UMAP1"]], y = separate.fcs.files[[i]][["UMAP2"]], col = as.factor(separate.fcs.files[[i]][["cluster"]]))
    print(ggplot(plot) + geom_point(aes(x=x, y=y, col = col)) +
            coord_fixed(ratio=graphical.ratio)+ labs(color = "FlowSOM Cluster", x = "UMAP1", y = "UMAP2", title = "FlowSOM Clustering on UMAP Axes",caption = newname) + 
            scale_color_manual(values = values) + guides(colour = guide_legend(override.aes = list(size=5)))+ theme_bw() + theme(plot.caption = element_text(size = 6)))}
  ggsave(filename = paste0("PerSubjectUMAP_number_",i,".pdf"))
}








# with help from https://github.com/cytolab/roe_phospho_flow/blob/master/01_data_analysis.Rmd


# COFACTOR_FOR_ARCSINH_SCALE = 5          # set cofactor for all channels
# NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE = 1000     # choose number of cells to sample per file 
## (if you want proportional sampling, enter a number larger than the max number of cells in a single file)
# TARGET_NUMBER_OF_CLUSTERS = 10          # choose target number of clusters for FlowSOM
# SEED = 1                                # set seed for reproducible results (can choose any number)

# fcs.files <- list.files(path = '.', pattern='.fcs$', full=TRUE, ignore.case=T)
# data.lists <- lapply(lapply(fcs.files, read.FCS), exprs)
# my.data = as.data.frame(do.call(rbind, mapply(cbind, data.lists, "File_ID" = c(1:length(data.lists)), SIMPLIFY = F)))
# orig_names = c(colnames(my.data),"UMAP1","UMAP2","cluster")
# colnames(my.data)[1:length(my.data) - 1] <- as.character(read.FCS(fcs.files[[1]])@parameters@data[["desc"]])

## equally sample the data 
# files.to.sample = split(my.data,my.data$`File_ID`)
# sampled.data <- list()
# for (i in 1: length(files.to.sample)){
#   if (nrow(files.to.sample[[i]])>NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE){            
#     sample.df =  files.to.sample[[i]]
#     sampled.data[[i]] = as.data.frame(sample.df[sample(nrow(sample.df), NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE), ])}        
#   else{
#     sampled.data[[i]] = files.to.sample[[i]]}}
# my.sampled.data = as.data.frame(do.call(rbind, sampled.data))   

## Run UMAP on chosen markers
## select all channels to use in UMAP by OPENING CONSOLE below
# set.seed(SEED)
# chosen.markers = as.data.frame(as.data.frame(my.sampled.data)[,c(choose.markers(my.sampled.data))])        # ************** INTERACTIVE STEP ********
# umap.markers <- chosen.markers %>%
#   mutate_all(function(x)
#     asinh(x / COFACTOR_FOR_ARCSINH_SCALE))
# myumap <- umap(umap.markers, ret_model = TRUE,verbose = TRUE)
# umap.data = as.data.frame(myumap$embedding)
# colnames(umap.data) <- c("UMAP1", "UMAP2")

## create flowFrame from UMAP data
# umap.matrix <- as.matrix(umap.data)
# UMAP.metadata <- data.frame(name = dimnames(umap.matrix)[[2]], desc = paste('UMAP', dimnames(umap.matrix)[[2]]))
# UMAP.metadata$range <- apply(apply(umap.matrix, 2, range), 2, diff)
# UMAP.metadata$minRange <- apply(umap.matrix, 2, min)
# UMAP.metadata$maxRange <- apply(umap.matrix, 2, max)
# umap.flowframe <- new("flowFrame", exprs=umap.matrix,parameters = AnnotatedDataFrame(UMAP.metadata))
## run FlowSOM on UMAP axes
# fSOM.umap <- FlowSOM(umap.flowframe, compensate = FALSE, transform = FALSE, toTransform=c(1:2), scale = TRUE, colsToUse = c(1:2), nClus = TARGET_NUMBER_OF_CLUSTERS, seed = SEED)
# FlowSOM.clusters <- as.numeric(as.vector(as.matrix(fSOM.umap[[2]][fSOM.umap[[1]]$map$mapping[,1]])))
# analysis.data = as.data.frame(cbind(my.sampled.data,umap.data,FlowSOM.clusters))
# colnames(analysis.data)[ncol(analysis.data)]<-"cluster"

## plot UMAP axes and FlowSOM per file used in analysis as well as for the concatenated data (last plot)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# col_vector = col_vector[-c(4,17,19,27,29:45)]  # take away certain olors
# values = sample(col_vector)
# range <- apply(apply(umap.data, 2, range), 2, diff)
# graphical.ratio <- (range[1]/range[2])

## plot FlowSOM clusters on UMAP axes for concatenated data
# ggplot(data.frame(x = analysis.data$UMAP1, y = analysis.data$UMAP2, col = as.factor(analysis.data$cluster))) + coord_fixed(ratio=graphical.ratio) + 
#   geom_point(aes(x=x, y=y, color=col),cex = 1.5) + guides(colour = guide_legend(override.aes = list(size=5), nrow = 13)) +
#   labs(x = "UMAP 1", y = "UMAP 2",title = "FlowSOM Clustering on UMAP Axes (concatenated)", color = "FlowSOM Cluster") + theme_bw() + scale_color_manual(values = values)

## plot FlowSOM clusters on UMAP axes for each file used in analysis
# separate.fcs.files = split(analysis.data,analysis.data$`File_ID`)
# for (i in 1:length(separate.fcs.files)){
#   newname  = str_remove(fcs.files[i], ".fcs")
#   plot <- data.frame(x = separate.fcs.files[[i]][["UMAP1"]], y = separate.fcs.files[[i]][["UMAP2"]], col = as.factor(separate.fcs.files[[i]][["cluster"]]))
#   print(ggplot(plot) + geom_point(aes(x=x, y=y, col = col)) +
#           coord_fixed(ratio=graphical.ratio)+ labs(color = "FlowSOM Cluster", x = "UMAP1", y = "UMAP2", title = "FlowSOM Clustering on UMAP Axes",caption = newname) + 
#           scale_color_manual(values = values) + guides(colour = guide_legend(override.aes = list(size=5)))+ theme_bw() + theme(plot.caption = element_text(size = 6)))}


## export new FCS files with UMAP1, UMAP2, and FlowSOM cluster ID as added channels 
## markers used for the UMAP will also have a (u) next to their channel name
# my.new.data  = analysis.data
# for (s in 1:ncol(my.new.data)){
#   colnames(my.new.data)[which(colnames(my.new.data) == colnames(umap.markers)[s])] <- paste0(colnames(analysis.data)[which(colnames(analysis.data) == colnames(umap.markers)[s])]," (u)")}
# new_desc = colnames(my.new.data)
# colnames(analysis.data)<-orig_names
# separate.files = split(analysis.data,as.factor(analysis.data$`File_ID`))
# for (i in 1:length(separate.files)){
#   single.file = separate.files[[i]]
#   remove.ID  = single.file[-c(ncol(my.sampled.data))]
#   mat <- as.matrix(single.file)
#   metadata <- data.frame(name = dimnames(mat)[[2]], desc = new_desc)
#   metadata$range <- apply(apply(mat, 2, range), 2, diff)
#   metadata$minRange <- apply(mat, 2, min)
#   metadata$maxRange <- apply(mat, 2, max)
#   export.flowframe <- new("flowFrame", exprs = mat, parameters = AnnotatedDataFrame(metadata))
#   newname  = str_remove(fcs.files[i], ".fcs")
#   filename = paste0(newname,"_UMAP_FlowSOM.fcs")
#   write.FCS(export.flowframe,filename = filename)
#   print(i)80 
# }