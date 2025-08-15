library(readxl)
library(tidyverse)
library(rstatix)
library(tibble)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

coord_file = args[1]
modules_file = args[2]
output_file = args[3]

promoters = read.table(coord_file)
colnames(promoters) <- c("chromosome", "start", "end", "gene_name")
studied = read.table(modules_file, header = FALSE)
colnames(studied) <- c("gene_name")
promoters <- subset(promoters, gene_name %in% studied$gene_name)

promoters <- promoters %>%
  filter(!is.na(start))

write.table(promoters,
          file = output_file,
          sep = "\t",          
          quote = FALSE,       
          col.names = FALSE,   
          row.names = FALSE)
