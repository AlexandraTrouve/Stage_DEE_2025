#####* ROC analysis for SVM models *#####
##### function used for ROC analysis #####
library("data.table")
library("ROCR")
library(tools)
library(ggplot2)
library(stringr)
rocPrArea<-function(positive, negative) {

  positive$state<-rep(1,times=nrow(positive))
  negative$state<-rep(0,times=nrow(negative))
  aucResult<-c()
  pr_plotValueX<-list()
  pr_plotValueY<-list()
  pos <- data.frame(V1 = positive[[2]], state = positive$state)
  neg <- data.frame(V1 = negative[[2]], state = negative$state)
  tempData <- rbind(pos, neg)
  pred <- prediction(tempData$V1,tempData$state) 
  perf <- performance( pred, "tpr", "fpr" )
  
  auc_temp <-performance( pred, measure = "auc")
  aucResult[1] <- round(unlist(slot(auc_temp, "y.values")), digits = 3)
    
  pr_temp<-data.frame(unlist(perf@x.values),unlist(perf@y.values))
  names(pr_temp)<-c("x","y")
  pr_temp<-na.omit(pr_temp)
  pr_plotValueX[[1]]<-pr_temp$x
  pr_plotValueY[[1]]<-pr_temp$y
  return(list(pr_plotValueX,pr_plotValueY,aucResult))}

args = commandArgs(trailingOnly=TRUE)
pos_files <- args[1:2]
neg_files <- args[3:4]
output_file <- args[5]

filtered_pos_files <- c()
filtered_neg_files <- c()
plotLegend = list()
file_name <- tools::file_path_sans_ext(basename(pos_files[[1]]))
model <- str_extract(file_name, "(?<=using_)[^.]+(?=_model)")
for (file_path in pos_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  cluster <- stringr::str_extract(file_name, "^(.*?)(?=_pos)")
  if (cluster != model) {filtered_pos_files <- c(filtered_pos_files, file_path)
  plotLegend <- append(plotLegend, cluster)}}
pos_files <- filtered_pos_files
for (file_path in neg_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  cluster <- stringr::str_extract(file_name, "^(.*?)(?=_neg)")
  if (!is.na(cluster) && cluster != model) {filtered_neg_files <- c(filtered_neg_files, file_path)}}
neg_files <- filtered_neg_files

##### take T1 enhancers as an example #####
## T1 prediction based on T1 model 
pos1<-fread(pos_files[[1]])
neg1<-fread(neg_files[[1]])

stage1<-rocPrArea(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
pdf(file = output_file, width = 8, height = 7)
plotColors<-c("pink","lightblue")
names(plotColors) <- c("high_variability","low_variability")

par(mar=c(5,5,2,2))
plot(unlist(aucOthers[[1]][1]),unlist(aucOthers[[2]][1]),col=plotColors[plotLegend[[1]]], lwd = 3,cex=1.5,main=paste0("zebrafish ",model," promoters"),xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.2,0.2,paste(plotLegend[[1]]," ","(AUC = ",round(aucOthers[[3]][1],digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[plotLegend[[1]]]) 

dev.off()
