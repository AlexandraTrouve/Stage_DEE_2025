library(readxl)
library(tidyverse)
library(rstatix)
library(tibble)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

coord_file = args[1]
studied_file = args[2]
results_file = args[3]

promoters = read.table(coord_file)
promoters <- promoters %>%
  distinct()
studied = read.table(studied_file, header = TRUE)
colnames(promoters) <- c("chromosome", "start", "end", "gene_name")
studied <- rownames_to_column(studied, var = "gene_name")
promoters <- merge(promoters, studied[, c("gene_name", "Label")], by = "gene_name", all.x = TRUE)
promoters$Label[is.na(promoters$Label)] <- "other_genes"
promoters <- promoters %>%
  relocate(gene_name, .after = end)
promoters <- promoters %>%
  filter(!is.na(start))

promoters_hourglass <- promoters[promoters$Label == "hourglass", ]
promoters_hourglass <- subset(promoters_hourglass, select=-Label)
promoters_inverse_hourglass <- promoters[promoters$Label == "inverse_hourglass", ]
promoters_inverse_hourglass <- subset(promoters_inverse_hourglass, select=-Label)
promoters_uniform_high <- promoters[promoters$Label == "uniform_high", ]
promoters_uniform_high <- subset(promoters_uniform_high, select=-Label)
promoters_uniform_low <- promoters[promoters$Label == "uniform_low", ]
promoters_uniform_low <- subset(promoters_uniform_low, select=-Label)
promoters_late_constraint <- promoters[promoters$Label == "late_constraint", ]
promoters_late_constraint <- subset(promoters_late_constraint, select=-Label)
promoters_early_constraint <- promoters[promoters$Label == "early_constraint", ]
promoters_early_constraint <- subset(promoters_early_constraint, select=-Label)
promoters_other_genes <- promoters[promoters$Label == "other_genes", ]
promoters_other_genes <- subset(promoters_other_genes, select=-Label)


write.table(promoters_hourglass,
          file = paste0(results_file,'hourglass.bed'),
          sep = "\t",          
          quote = FALSE,       
          col.names = FALSE,   
          row.names = FALSE)

write.table(promoters_inverse_hourglass,
            file = paste0(results_file,'inverse_hourglass.bed'),
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)

write.table(promoters_early_constraint,
            file = paste0(results_file,'early_constraint.bed'),
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)

write.table(promoters_late_constraint,
            file = paste0(results_file,'late_constraint.bed'),
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)

write.table(promoters_uniform_high,
            file = paste0(results_file,'uniform_high.bed'),
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)

write.table(promoters_uniform_low,
            file = paste0(results_file,'uniform_low.bed'),
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)

write.table(promoters_other_genes,
            file = paste0(results_file,'other_genes.bed'),
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)
