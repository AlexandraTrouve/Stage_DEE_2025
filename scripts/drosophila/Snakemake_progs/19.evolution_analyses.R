library(tidyverse)
library(rstatix)
library(ggplot2)
library(tibble)
library(dplyr)
library(ggpubr)
library(car)
args = commandArgs(trailingOnly=TRUE)

SI_file = args[1]
phastcons_file = args[2]
phylop_file = args[3]
output_file = args[4]

phastcons_scores <- read.csv(phastcons_file, header=TRUE, sep=',')
colnames(phastcons_scores) <- c("gene_name", "phastCons_score")
phylop_scores <- read.csv(phylop_file, header=TRUE, sep=',')
colnames(phylop_scores) <- c("gene_name", "phyloP_score")

shape_index <- read.csv(SI_file, header=TRUE, sep='')
shape_index <- subset(shape_index, select=-c(chromosome, start, end))
shape_index <- merge(shape_index, phastcons_scores[,c("gene_name", "phastCons_score")], by="gene_name",all.y=TRUE)
shape_index <- merge(shape_index, phylop_scores[,c("gene_name", "phyloP_score")], by="gene_name",all.y=TRUE)
shape_index <- shape_index[!is.na(shape_index$Label), ]

compare_means(
  formula = phastCons_score ~ Label,       
  data = shape_index,               
  method = "wilcox.test",
  ref.group = "other_genes")
compare_means(
  formula = phyloP_score ~ Label,       
  data = shape_index,               
  method = "wilcox.test",
  ref.group = "other_genes")
shape_index %>% wilcox_effsize(phastCons_score ~ Label)
shape_index %>% wilcox_effsize(phyloP_score ~ Label)

shape_index$Label <- factor(shape_index$Label, levels = c("other_genes", "hourglass", "inverse_hourglass", 
                                                          "early_constraint", "late_constraint", 
                                                          "uniform_high", "uniform_low"))

pdf(file = output_file, width = 10, height = 9)
ggplot(data = shape_index, aes(x = Label, y = phastCons_score, fill = Label)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, notch = TRUE) +
  scale_fill_manual(values = c("other_genes" = "darksalmon", 
                               "hourglass" = "lightgrey", 
                               "inverse_hourglass" = "lightgrey", 
                               "early_constraint" = "lightgrey", 
                               "late_constraint" = "lightgrey", 
                               "uniform_high" = "lightgrey", 
                               "uniform_low" = "lightgrey")) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), hide.ns = TRUE, ref.group = "other_genes") +
  labs(title = "Promoter phastCons score in terms of cluster",
    x = "Clusters",
    y = "Promoter phastCons score",
    fill = "Cluster") +
  theme_minimal()

ggplot(data = shape_index, aes(x = Label, y = phyloP_score, fill = Label)) +
  geom_violin() +
  geom_boxplot(width = 0.2, outlier.shape = NA, notch = TRUE) +
  scale_fill_manual(values = c("other_genes" = "darksalmon", 
                               "hourglass" = "lightgrey", 
                               "inverse_hourglass" = "lightgrey", 
                               "early_constraint" = "lightgrey", 
                               "late_constraint" = "lightgrey", 
                               "uniform_high" = "lightgrey", 
                               "uniform_low" = "lightgrey")) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), hide.ns = TRUE, ref.group = "other_genes") +
  labs(title = "Promoter phyloP score in terms of cluster",
    x = "Clusters",
    y = "Promoter phyloP score",
    fill = "Cluster") + 
  theme_minimal()
dev.off()


