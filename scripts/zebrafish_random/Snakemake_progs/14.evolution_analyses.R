library(tidyverse)
library(rstatix)
library(ggplot2)
library(tibble)
library(dplyr)
library(ggpubr)
library(car)
args = commandArgs(trailingOnly=TRUE)

modules_files = args[1:2]
phastcons_file = args[3]
phylop_file = args[4]
output_file = args[5]

#cluster1 = '/Users/children/Desktop/Stage_2025/zebrafish/Data/DRE_bias_genes_high_variability.txt'
#cluster2 = '/Users/children/Desktop/Stage_2025/zebrafish/Data/DRE_bias_genes_low_variability.txt'
#modules_files <- list(cluster1, cluster2)
#phastcons_file = '/Users/children/Desktop/Stage_2025/zebrafish/Results/analyses/evolution_scores/phastCons/mean_scores.csv'
#phylop_file = '/Users/children/Desktop/Stage_2025/zebrafish/Results/analyses/evolution_scores/phyloP/mean_scores.txt'
#output_file = '/Users/children/Desktop/Stage_2025/zebrafish/Results/analyses/evolution_scores/analyses1.pdf'

phastcons_scores <- read.csv(phastcons_file, header=TRUE, sep=',')
colnames(phastcons_scores) <- c("gene_name", "phastCons_score")
phylop_scores <- read.csv(phylop_file, header=TRUE, sep=',')
colnames(phylop_scores) <- c("gene_name", "phyloP_score")
scores <- merge(phastcons_scores, phylop_scores[,c("gene_name", "phyloP_score")], by="gene_name",all.y=TRUE)
  
cluster1 <- read.table(modules_files[[1]], header=FALSE)
file_name <- tools::file_path_sans_ext(basename(modules_files[[1]]))
pieces <- strsplit(file_name, "_")[[1]]
cluster_name <- paste(pieces[4], pieces[5], sep = "_")
cluster1$Label<-rep(cluster_name,times=nrow(cluster1))
colnames(cluster1) <- c('gene_name', 'Label')
cluster2 <- read.table(modules_files[[2]], header=FALSE)
file_name <- tools::file_path_sans_ext(basename(modules_files[[2]]))
pieces <- strsplit(file_name, "_")[[1]]
cluster_name <- paste(pieces[4], pieces[5], sep = "_")
cluster2$Label<-rep(cluster_name,times=nrow(cluster2))
colnames(cluster2) <- c('gene_name', 'Label')
clusters <- rbind(cluster1, cluster2)

scores <- merge(scores, clusters[, c("gene_name", "Label")], by="gene_name", all.x = TRUE)
scores$Label[is.na(scores$Label)] <- "other_genes"
#scores <- subset(scores, phyloP_score != 0)

leveneTest(phastCons_score ~ Label, data = scores)
leveneTest(phyloP_score ~ Label, data = scores)
pc_aov_result <- aov(phastCons_score ~ Label, data = scores)
pp_aov_result <- aov(phyloP_score ~ Label, data = scores)
TukeyHSD(pc_aov_result)
TukeyHSD(pp_aov_result)

pdf(file = output_file, width = 8, height = 7)
ggplot(data = scores, aes(x = Label, y = phyloP_score)) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), hide.ns = TRUE, ref.group = "other_genes") +
  geom_violin() +
  geom_boxplot(width = 0.2, outliers=F, notch=T) +
  coord_cartesian(ylim=c(0, 1)) +
  labs(
    title = "Promoter phyloP score in terms of cluster",
    x = "Clusters",
    y = "Promoter phyloP score")
ggplot(data = scores, aes(x = Label, y = phastCons_score)) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), hide.ns = TRUE, ref.group = "other_genes") +
  geom_violin() +
  geom_boxplot(width = 0.2, outliers=F, notch=T) +
  labs(
    title = "Promoter phastCons score in terms of cluster",
    x = "Clusters",
    y = "Promoter phastCons score")
dev.off()


