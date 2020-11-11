# library(cytofkit2)
Sys.setenv("MC_CORES" = 8L)
library("parallel")
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
Spectre::package.check()
Spectre::package.load()
library(dplyr)
library(rstatix)


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


COFACTOR_FOR_ARCSINH_SCALE <- 1000          # set cofactor for all channels
NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE <- 10000     # choose number of cells to sample per file 
# (if you want proportional sampling, enter a number larger than the max number of cells in a single file)
TARGET_NUMBER_OF_CLUSTERS <- 20          # choose target number of clusters for FlowSOM
SEED <- 1                                # set seed for reproducible results (can choose any number)
origDir <- getwd()
setwd(dir = "D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/exportNonnavB/d0")
fcs.files <- list.files(path = '.', pattern='.fcs$', full=TRUE, ignore.case=T)
data.lists <- lapply(lapply(fcs.files, read.FCS), exprs)                          # expression data, one list element per file
data.lists <- lapply(data.lists, function(df) tibble::rownames_to_column(as.data.frame(df),var="cellID"))
my.data <- as.data.frame(do.call(rbind, mapply(cbind, data.lists, "File_ID" = c(1:length(data.lists)), SIMPLIFY = F)))
orig_names <- c(colnames(my.data),"UMAP1","UMAP2","cluster")
colnames(my.data)[1:length(my.data) - 1] <- as.character(read.FCS(fcs.files[[1]])@parameters@data[["desc"]])

# equally sample the data 
files.to.sample <- split(my.data,my.data$`File_ID`)
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
chosen.markers <- as.data.frame(as.data.frame(my.sampled.data)[,c(choose.markers(my.sampled.data))])        # ************** INTERACTIVE STEP *****************
#  chosen markers:  CD27, Bcl6, IgD, CD32, IgM, CD86, CD138, CD38, Blimp1, CD11c, CD71, CXCR4, CD21, CD14, CXCR5, CD16, CD23, CD80, Tbet, Ki67, CD20         7:11,15:32
umap.markers <- chosen.markers %>%   mutate_all(function(x)    asinh(x / COFACTOR_FOR_ARCSINH_SCALE))
myumap <- umap(umap.markers, ret_model = TRUE,verbose = TRUE)
umap.data <- as.data.frame(myumap$embedding)
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
analysis.data <- as.data.frame(cbind(my.sampled.data,umap.data,FlowSOM.clusters))
colnames(analysis.data)[ncol(analysis.data)]<-"cluster"

# plot UMAP axes and FlowSOM per file used in analysis as well as for the concatenated data (last plot)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[-c(4,17,19,27,29:45)]  # take away certain olors
values <- sample(col_vector)
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
         # filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_clusterComparison_byCohort_nonnavB_d0_heatmap.pdf",
         )  
dev.off()

propClusters %>% group_by(cluster) %>% rstatix::t_test(data= ., formula =  freq ~ Cohort) %>% rstatix::adjust_pvalue(method = 'fdr')
# propClusters %>% group_by(cluster==1) %>% rstatix::t_test(data= ., formula =  freq ~ Cohort)
# propClusters %>% group_by(cluster==18) %>% rstatix::t_test(data= ., formula =  freq ~ Cohort)

writeFCSfiles( analysis.data, umap.markers, orig_names)


##  *********************************************************************************************************************************************************
##  *********************************************************************************************************************************************************
##  *********************************************************************************************************************************************************


COFACTOR_FOR_ARCSINH_SCALE <- 5          # set cofactor for all channels
NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE <- 20000     # choose number of cells to sample per file 
TARGET_NUMBER_OF_CLUSTERS <- 25          # choose target number of clusters for FlowSOM
SEED <- 1 ; set.seed(SEED)                               # set seed for reproducible results (can choose any number)
origDir <- getwd()
setwd(dir = "D:/Pembro-Fluvac/18-19season/Flow cytometry/Tcell/exportNonNavCD4/")
fcs.files <- list.files(path = '.', pattern='.fcs$', full=TRUE, ignore.case=T)
data.lists <- lapply(lapply(fcs.files, read.FCS), exprs)                          # expression data, one list element per file
data.lists <- lapply(data.lists, function(df) tibble::rownames_to_column(as.data.frame(df),var="cellID"))
data.lists <- lapply(data.lists, function(df) df[,c(2:ncol(df), 1)])              # nove cellID column to the end 
my.data <- as.data.frame(do.call(rbind, mapply(cbind, data.lists, "File_ID" = c(1:length(data.lists)), SIMPLIFY = F)))
orig_names <- c(colnames(my.data),"UMAP1","UMAP2","cluster")
colnames(my.data)[1:length(my.data) - 1] <- as.character(read.FCS(fcs.files[[1]])@parameters@data[["desc"]])        # had to change this from 1:len... to 3:len...

# equally sample the data 
files.to.sample <- split(my.data,my.data$`File_ID`)
sampled.data <- list()
for (i in 1: length(files.to.sample)){
  if (nrow(files.to.sample[[i]])>NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE){            
    sample.df =  files.to.sample[[i]]
    sampled.data[[i]] = as.data.frame(sample.df[sample(nrow(sample.df), NUMBER_OF_CELLS_TO_SAMPLE_PER_FILE), ])}        
  else{
    sampled.data[[i]] = files.to.sample[[i]]}}
my.sampled.data <- as.data.frame(do.call(rbind, sampled.data))   

# Run UMAP on chosen markers
# select all channels to use in UMAP by OPENING CONSOLE below
set.seed(SEED)
chosen.markers <- as.data.frame(as.data.frame(my.sampled.data)[,c(choose.markers(my.sampled.data))])        # ************** INTERACTIVE STEP *****************
#  chosen markers:  CD27, CTLA4, SLAM, CD45RA, CD28, Tcf1, CD38, CD25, CXCR4, CD127, CXCR3, Foxp3, CXCR5, CCR6, CD36, Ki67, ICOS        7:8,10,16:26,28:30
umap.markers <- chosen.markers %>%   mutate_all(function(x)    asinh(x / COFACTOR_FOR_ARCSINH_SCALE))
myumap <- umap(umap.markers, ret_model = TRUE,verbose = TRUE)
umap.data <- as.data.frame(myumap$embedding)
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
analysis.data <- as.data.frame(cbind(my.sampled.data,umap.data,FlowSOM.clusters))
colnames(analysis.data)[ncol(analysis.data)]<-"cluster"

# plot UMAP axes and FlowSOM per file used in analysis as well as for the concatenated data (last plot)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[-c(4,17,19,27,29:45)]  # take away certain colors
values <- sample(col_vector)
range <- apply(apply(umap.data, 2, range), 2, diff)
graphical.ratio <- (range[1]/range[2])

# plot FlowSOM clusters on UMAP axes for concatenated data                                          # pch=21, color='black',stroke=0.05,
ggplot(data.frame(x = analysis.data$UMAP1, y = analysis.data$UMAP2, col = as.factor(analysis.data$cluster))) + coord_fixed(ratio=graphical.ratio) + 
  geom_point(aes(x=x, y=y, color=col),size = 3, alpha=0.75) + guides(colour = guide_legend(override.aes = list(size=5), nrow = 13)) + 
  labs(x = "UMAP 1", y = "UMAP 2",title = "FlowSOM Clustering for Nonnaive T cells at baseline", color = "FlowSOM Cluster") + theme_bw() + scale_color_manual(values = values)
# ggsave( filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_UMAP_nonnavT_d0_20kcells_25clusters.pdf")

# separateUMAPplots(analysis.data, fcs.files)


propClusters <- analysis.data[,-c(1:6, 32:33)] %>% group_by(File_ID, cluster)  %>% summarize(n=n())  %>% mutate(freq = n/sum(n))
propClusters$n <- NULL
propClusters$File_ID <- sapply(propClusters$File_ID, function(x) { fcs.files[x] } )
propClusters$Cohort <- "HC"; propClusters$Cohort[grep("19616",propClusters$File_ID, value = F)] <- "aPD1"
propClusters$cluster <- factor(propClusters$cluster)
a <- ggplot(propClusters, aes(x=cluster, y=freq, group=Cohort, color=Cohort)) +  stat_smooth(method="loess",span=0.1,se=T, aes(fill=Cohort), alpha=0.1) +  
  ggtitle("Nonnaive T cells at baseline") + ylab("Frequency (% nonnaive T)") + xlab("FlowSOM cluster") + theme_bw() +
  theme(axis.text.x = element_text(angle=45,vjust = 0.9, hjust=1), axis.text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=24),
        legend.text = element_text(size=12), legend.title = element_text(size=14)); a
# ggsave(filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_clusterComparison_byCohort_nonnavT_d0.pdf")
clusterSummary <- analysis.data[,-c(1:6,32:36)] %>% group_by(cluster)  %>% # mutate_all(function(x)  asinh(x / COFACTOR_FOR_ARCSINH_SCALE))%>% 
  summarize_all(list(mean)) 
clusterSummary <- clusterSummary %>% column_to_rownames("cluster")
clusterSummary <- clusterSummary[,-grep(paste0(c("CD19","dead","CD20"), collapse='|'), colnames(clusterSummary),value=F)]

# assignInNamespace(   x = "draw_colnames",   value = "draw_colnames_90",   ns = asNamespace("pheatmap") )

b <- pheatmap(clusterSummary, scale = "column", cluster_rows = F, color=viridis::inferno(100), main="Cluster expression for all subjects", fontsize=12, 
         border_color = "grey30", angle_col = 90,
         #filename = "D:/Pembro-Fluvac/Analysis/Images/FlowSOM_clusterComparison_byCohort_nonnavT_d0_heatmap.pdf",
); b
dev.off()

res <- propClusters %>% group_by(cluster) %>% rstatix::t_test(data= ., formula =  freq ~ Cohort) %>% rstatix::adjust_pvalue(method = 'fdr')
res <- propClusters %>% group_by(cluster, Cohort) %>% rstatix::get_summary_stats(data= .) # %>% rstatix::adjust_pvalue(method = 'fdr')

writeFCSfiles( analysis.data, umap.markers, orig_names)

# rm(my.sampled.data, my.sampled.data2, my.data)
# placeholder


##  *********************************************************************************************************************************************************
##  *********************************************************************************************************************************************************
##  *********************************************************************************************************************************************************

# library('devtools')
# install_github("immunedynamics/spectre")

PrimaryDirectory <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/exportNonnavB/d0/spectre/data");  InputDirectory <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/exportNonnavB/d0/spectre/metadata");  MetaDirectory <- getwd()
setwd("D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/exportNonnavB/d0/spectre/output");  OutputDirectory <- getwd()
data.list <- Spectre::read.files(file.loc = InputDirectory, file.type = '.fcs', do.embed.file.names = T)
do.list.summary(data.list)
exclude <- grep(pattern = paste(c("029_d0"),collapse='|'), names(data.list))
minCells <- min(unlist(do.list.summary(data.list[-exclude])$nrow.check))  # need to downsample so that there is equal representation of each sample
temp <- lapply(data.list[-exclude], function(x) do.subsample(x, minCells, seed=1))

cell.dat <- Spectre::do.merge.files(dat = temp)          
setwd(InputDirectory);  oneFCS <- flowCore::read.FCS(filename = "export_PBMCs_19616-017_d0_001_NonnaiveB.fcs")
namesFCS <- oneFCS@parameters$desc
namesFCS[1:6] <- names(cell.dat)[1:6]         # manually confirmed these are the same parameters and in the same order
names(cell.dat)[1:32] <- namesFCS[1:32]       # move the human-readable labels over to cell.dat after manually confirming channels are in filter + alphabetical order
names(cell.dat)
cell.dat <- do.asinh(dat = cell.dat, use.cols = names(cell.dat)[c(7:32)], cofactor = 1000)
transformed.cols <- paste0(names(cell.dat)[c(7:32)], "_asinh")
# for(i in transformed.cols){   make.colour.plot(do.subsample(cell.dat, 20000), i, "Bcl6") }
setwd(MetaDirectory)
meta.data <- fread(input = "sample.details_fcs.csv")
cell.dat <- do.add.cols(dat = cell.dat, base.col = "FileName", add.dat = meta.data, add.by = "Filename", rmv.ext = TRUE)
cellular.cols <- names(cell.dat)[c(36:40,42,44:55,57:61)]            # choose the arcsinh-transformed columns as the source data
#  chosen markers:  CD27, Bcl6, IgD, CD32, IgM, CD86, CD138, CD38, Blimp1, CD11c, CD71, CXCR4, CD21, CXCR5, CD16, CD23, CD80, Tbet, Ki67, CD20        
# cluster.cols <- names(cell.dat)[c(36:40,44:52,54:55,57:61)]             # specify which columns to cluster
cluster.cols <- names(cell.dat)[grep(pattern = paste0(c("CD86_asinh", "CD71_asinh","CD38_asinh", "Tbet_asinh",
                                                        "CD20_asinh","Ki67_asinh","CXCR5_asinh","CXCR4_asinh",
                                                        "CD21_asinh","CD138_asinh","IgD_asinh", "CD11c_asinh"), collapse = '|'), names(cell.dat), value = F)]  
exp.name <- "CNS experiment"; sample.col <- "Sample"; group.col <- "Group"; batch.col <- "Batch"
data.frame(table(cell.dat[[group.col]]))
cell.flowSOM <- run.flowsom(cell.dat, cluster.cols, meta.k = 15, meta.seed=1, clust.seed = 1)       # renaming now to cell.flowSOM in case of re-clustering later
sub.targets <- c(5000, 5000)
cell.sub <- do.subsample(dat = cell.flowSOM, sub.targets, group.col, seed = 1)
cell.sub <- run.umap(cell.sub, cluster.cols)
setwd(OutputDirectory)
make.colour.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE, filename = "UMAP_FlowSOM_k20_allSubj.pdf")
make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", cluster.cols, filename="UMAP_FlowSOM_perChannel_overlay.pdf")

make.multi.plot(cell.sub, "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')
make.colour.plot(subset(cell.sub, Group=="Healthy"), "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE, filename = "UMAP_FlowSOM_k20_Healthy.pdf")
make.colour.plot(subset(cell.sub, Group=="aPD1"), "UMAP_X", "UMAP_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE, filename = "UMAP_FlowSOM_k20_aPD1.pdf")

exp <- do.aggregate(cell.flowSOM, cellular.cols, by = "FlowSOM_metacluster")
make.pheatmap(exp, "FlowSOM_metacluster", cellular.cols, standard.colours = "inferno", file.name = "Pheatmap_flowSOM_metacluster.pdf" )
dev.off()
# write.sumtables(dat = cell.flowSOM, sample.col = sample.col,  measure.col = cellular.cols, annot.col = c(group.col, batch.col), group.col = group.col, do.proportions = TRUE, do.mfi.per.sample = FALSE, do.mfi.per.marker = TRUE)

propClust.Spect <-cell.flowSOM[,-c(1:33)] %>% group_by(FileName, FlowSOM_metacluster) %>% summarize(n=n()) %>% mutate(freq = n/sum(n))
# propClusters <- analysis.data[,-c(1:6,14,33:34,36:37)] %>% group_by(File_ID, cluster)  %>% summarize(n=n())  %>% mutate(freq = n/sum(n))
propClust.Spect$n <- NULL
# propClusters$File_ID <- sapply(propClusters$File_ID, function(x) { fcs.files[x] } )
propClust.Spect$Cohort <- "HC"; propClust.Spect$Cohort[grep("19616",propClust.Spect$FileName, value = F)] <- "aPD1"
propClust.Spect$FlowSOM_metacluster <- factor(propClust.Spect$FlowSOM_metacluster)
ggplot(propClust.Spect, aes(x=FlowSOM_metacluster, y=freq, group=Cohort, color=Cohort)) + # stat_smooth(method="loess",span=0.1,se=TRUE, aes(fill=Cohort), alpha=0.1) + 
  theme_bw() + stat_summary(fun="median", geom="bar", position = position_dodge(width=1), aes(fill=Cohort)) + geom_jitter(position = position_dodge(width=1), color='black') + 
  #geom_bar(stat = 'summary', position=position_dodge(width=1), aes(fill=Cohort)) + geom_jitter(position = position_dodge(width=1), color='black') + 
  ggtitle("Nonnaive B cells at baseline") + ylab("Frequency (% nonnaive B") + xlab("FlowSOM cluster") + 
  scale_fill_manual(values = c("#FFB18C", "#7FAEDB")) + scale_color_manual(values = c("#FFB18C", "#7FAEDB" )) + 
  theme(axis.text.x = element_text(angle=45,vjust = 0.9, hjust=1), axis.text = element_text(size=12), axis.title = element_text(size=14), title = element_text(size=24),
        legend.text = element_text(size=12), legend.title = element_text(size=14))
ggsave(filename = "D:/Pembro-Fluvac/18-19season/Flow cytometry/Bcell/exportNonnavB/d0/spectre/output/NonnavB_baseline_flowSOM.pdf", width=8)

tableOfResults <- propClust.Spect %>% group_by(FlowSOM_metacluster, Cohort) %>% rstatix::get_summary_stats() #  %>% print(n=50)
unstabClusters <- which(colSums(as.matrix(table(tableOfResults$Cohort, tableOfResults$FlowSOM_metacluster))) != 2)  # where is there complete absence of cluster in cohort
ifelse(length(unstabClusters) < 1 , 
       forStats <- propClust.Spect, 
       forStats <- propClust.Spect[-which(propClust.Spect$FlowSOM_metacluster == unstabClusters), ])
forStats  %>% group_by(FlowSOM_metacluster) %>% rstatix::wilcox_test(data= ., formula =  freq ~ Cohort) %>% rstatix::adjust_pvalue(method = 'fdr') %>% print(n=50)







##  *********************************************************************************************************************************************************
##  *********************************************************************************************************************************************************
##  *********************************************************************************************************************************************************



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