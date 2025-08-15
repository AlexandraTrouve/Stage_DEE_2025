library(readxl)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(tibble)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

scores_file = args[1]
SI_file = args[2]
output_file = args[3]

scores <- read.csv(scores_file, header = TRUE, sep = ',')
SI <- read.table(SI_file, header = TRUE)

SI <- merge(SI, scores[, c("gene_name", "SVM_score")], by = "gene_name")
SI <- subset(SI, select = -c(chromosome, start, end))

pdf(file = output_file, width = 8, height = 7)
ggplot(data = SI, aes(x = shape_index, y = SVM_score, color = Label)) +
  scale_color_manual(values = c("hourglass" = "pink", "inverse_hourglass" = "lightsalmon",
                               "early_constraint" = "lightblue", "late_constraint" = "thistle2",
                               "uniform_high" = "darkseagreen3", "uniform_low" = "snow3",
                               'other_genes' = 'darkgoldenrod1')) +
  geom_point() +
  labs(
    title = "SVM score in terms of promoter shape index",
    x = "Promoter shape index",
    y = "Promoter SVM score",
    color = "Cluster") + 
  theme_minimal()
dev.off()
