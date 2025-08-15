library(readxl)
library(tidyverse)
library(rstatix)
library(ggResidpanel)
library(tibble)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

ucsc_file = args[1]
correspondance_file = args[2]
ensembl_file = args[3]

correspondance <- read.table(correspondance_file, header = TRUE, sep = ',')
genes <- read.table(ucsc_file)
colnames(genes) <- c("chromosome", "start", "end", "gene_id")
colnames(correspondance) <- c("gene_name","gene_id")
genes <- merge(genes, correspondance[, c("gene_name", "gene_id")], by = "gene_id")
genes <- subset(genes, select=-gene_id)

write.table(genes,
            file = ensembl_file,
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)
