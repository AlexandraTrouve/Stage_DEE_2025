#####* ROC analysis for SVM models *#####
##### function used for ROC analysis #####
library("data.table")
library("ROCR")
library(tools)
library(ggplot2)
library(stringr)
rocPrArea<-function(positive, negative) {
  
  positive$state<-rep(1,tiems=nrow(positive))
  negative$state<-rep(0,tiems=nrow(negative))
  aucResult<-c()
  pr_plotValueX<-list()
  pr_plotValueY<-list()
  
  for (i in 2:4) {
    pos<-cbind(positive[[i]],positive[,5])
    neg<-cbind(negative[[i]],negative[,5])
    tempData<-rbind(pos,neg)
    pred <- prediction(tempData$V1,tempData$state) 
    perf <- performance( pred, "tpr", "fpr" )
    
    auc_temp <-performance( pred, measure = "auc")
    aucResult[i-1]<-round(unlist(slot(auc_temp, "y.values")),digits = 3)
    
    pr_temp<-data.frame(unlist(perf@x.values),unlist(perf@y.values))
    names(pr_temp)<-c("x","y")
    pr_temp<-na.omit(pr_temp)
    pr_plotValueX[[i-1]]<-pr_temp$x
    pr_plotValueY[[i-1]]<-pr_temp$y
  } 
  return(list(pr_plotValueX,pr_plotValueY,aucResult))
}

args = commandArgs(trailingOnly=TRUE)
pos_files <- args[1:3]
neg_files <- args[4:6]
output_file <- args[7]

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
## T2 prediction based on T1 model 
pos2<-fread(pos_files[[2]])
neg2<-fread(neg_files[[2]])

pos<-cbind(pos1,pos2$V2)
neg<-cbind(neg1,neg2$V2)
stage1<-rocPrArea(pos,neg)
aucOthers<-list(stage1[[1]],stage1[[2]],stage1[[3]])

## plot
pdf(file = output_file, width = 8, height = 7)
plotColors<-c("lightseagreen","lightcoral","cornflowerblue")
names(plotColors) <- c("sample1","sample1","sample3")

par(mar=c(5,5,2,2))
plot(unlist(aucOthers[[1]][1]),unlist(aucOthers[[2]][1]),col=plotColors[plotLegend[[1]]], lwd = 3,cex=1.5,main=paste0("ROC curves for model train on ",model," promoters"),xlab="False positive rate",ylab="Ture positive rate",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5,type="l")
legend(0.3,0.15,paste(plotLegend[[1]]," ","(AUC = ",round(aucOthers[[3]][1],digits = 3),")",sep=""),lwd = 3,bty="n",cex=1.4,col = plotColors[plotLegend[[1]]]) 

for (i in 2:3) {
  lines(unlist(aucOthers[[1]][i]),unlist(aucOthers[[2]][i]),col=plotColors[plotLegend[[i]]],lwd=3)}

legendY<-c(0.1, 0.05)
for (i in 2:3) {
  legendName<-paste(plotLegend[[i]]," ","(AUC = ",round(aucOthers[[3]][i],digits = 3),")",sep="")
  legend(0.3,legendY[i-1],legendName,lwd = 3,bty="n",cex=1.4,col = plotColors[plotLegend[[i]]])}

dev.off()
