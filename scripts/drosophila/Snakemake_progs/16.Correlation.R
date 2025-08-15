library(tidyverse)
library(rstatix)
library(ggplot2)
library(tibble)
library(dplyr)
library(ggpubr)
args = commandArgs(trailingOnly=TRUE)

SI_file = args[1]
promoters_scores_file = args[2]
table_promoters_scores_file = args[3]
boxplot_output_file = args[4]
correlation_output_file = args[5]

SI <- read.table(SI_file, header = TRUE)

compare_means(
  formula = shape_index ~ Label,       
  data = SI,               
  method = "wilcox.test",
  ref.group = "other_genes")
SI %>% wilcox_effsize(shape_index ~ Label)
SI_lim <- -1

pdf(file = boxplot_output_file, width = 10, height = 9)
ggplot(data = SI, aes(x = Label, y = shape_index, fill = Label)) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), hide.ns = TRUE, ref.group = "other_genes") +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, notch = TRUE) +
  scale_fill_manual(values = c("other_genes" = "darksalmon", 
                               "hourglass" = "lightgrey", 
                               "inverse_hourglass" = "lightgrey", 
                               "early_constraint" = "lightgrey", 
                               "late_constraint" = "lightgrey", 
                               "uniform_high" = "lightgrey", 
                               "uniform_low" = "lightgrey")) +
  geom_hline(yintercept = SI_lim,
             linetype = "dotted",
             color = "black",
             size = 1) +
  labs(title = "Promoter shape index in terms of cluster",
    x = "Clusters",
    y = "Promoter shape index") +
  theme_minimal()

dev.off()

SI <- subset(SI, select=-c(chromosome, start, end))
SI <- subset(SI, Label != "other_genes")

promoters <- read.csv(promoters_scores_file, header = TRUE, sep = ',')
colnames(promoters) <- c("gene_name", "promoter_scores")
table <- read.csv(table_promoters_scores_file, header = TRUE, sep = ',' )
colnames(table) <- c("gene_name", "table_promoter_scores")
SI <- merge(SI, promoters[, c("gene_name", "promoter_scores")], by = "gene_name", all.y = TRUE)
SI <- merge(SI, table[, c("gene_name", "table_promoter_scores")], by = "gene_name", all.y = TRUE)

pdf(correlation_output_file, width = 8, height = 7)
ggplot(SI, aes(x = promoter_scores, y = table_promoter_scores)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, col = "purple", fill = "pink") +
  stat_regline_equation(label.x = min(SI$promoter_scores) + 5, label.y = max(SI$table_promoter_scores) - 5) + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = min(SI$promoter_scores) + 5, label.y = max(SI$table_promoter_scores) - 50) +
  labs(title = "Linear regression between table promoters scores and our promoters scores",
       x = "Our promoters scores",
       y = "Table promoters scores") +
  theme_minimal()

ggplot(SI, aes(x = table_promoter_scores, y = shape_index)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, col = "black", fill = "grey") +
  stat_regline_equation(label.x = -200, label.y = -7.3) + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = -200, label.y = -7) +
  labs(title = "Linear regression between table promoters scores and promoters shape index",
       x = "Table promoters scores",
       y = "Promoters shape index") +
  theme_minimal()

ggplot(SI, aes(x = promoter_scores, y = shape_index)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, col = "orange", fill = "yellow") +
  stat_regline_equation(label.x = -100, label.y = -5.5) + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = -100, label.y = -5.2) +
  labs(title = "Linear regression between table promoters scores and promoters shape index",
       x = "Promoters scores",
       y = "Promoters shape index") +
  theme_minimal()
dev.off()
