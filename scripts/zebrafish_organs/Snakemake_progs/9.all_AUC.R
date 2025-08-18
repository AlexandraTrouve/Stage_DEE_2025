library("data.table")
library("ROCR")
library(tools)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

input_files <- args[1:8]
output_file <- args[9]

perfs <- list()
AUC_values <- list()
samples <- c()

for (file_path in input_files) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  cluster <- str_extract(file_name, "^(.*?)(?=.cvpred)")
  samples <- c(samples, cluster)
  cv <- fread(file_path, header = FALSE, sep = "\t")
  colnames(cv) <- c("position", "prediction", "real_state", "cv_number")
  pred <- prediction(cv$prediction, cv$real_state)
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, measure = "auc")
  AUC_val <- signif(unlist(slot(auc, "y.values")), digits = 3)
  perfs[[cluster]] <- perf
  AUC_values[[cluster]] <- AUC_val}

pdf(file = output_file, width = 8, height = 7)
plotColors<-c("thistle2","pink","lightblue","darkseagreen3","lightsalmon","snow3","brown2","lightgoldenrod2")
names(plotColors) <- c("brain","ovary","testis","eye","liver","pectoral_fin","gills","intestine")

plot(perfs[[1]]@x.values[[1]], perfs[[1]]@y.values[[1]], lwd = 1, cex = 0.05, col = plotColors[samples[1]],
     xlab = "False positive rate", ylab = "True positive rate",
     xlim = c(0, 1), ylim = c(0, 1))
title(main = "ROC Curves", cex.main = 1.8)
abline(0, 1, col = "red", lty = 2, lwd = 1)

AUC_values_for_legend <- c(paste0(samples[1], " (AUC = ", round(AUC_values[[1]], 3), ")"))
for (i in 1:length(perfs)) {
  if (i > 1) {
    lines(perfs[[i]]@x.values[[1]], perfs[[i]]@y.values[[1]],
          col = plotColors[samples[i]], lwd = 1)}
  AUC_values_for_legend[i] <- paste0(samples[i], " (AUC = ", round(AUC_values[[i]], 3), ")")}

legend("bottomright",
       legend = AUC_values_for_legend,
       col = plotColors[1:length(perfs)],
       lwd = 1,
       cex = 1.2,
       bty = "n")

dev.off()
