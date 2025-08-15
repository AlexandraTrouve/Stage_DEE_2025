# library iportation
library(readxl)
library(tidyverse)
library(rstatix)
library(tibble)
library(dplyr)
library(ggResidpanel)
args = commandArgs(trailingOnly=TRUE)

overlap_file = args[1]
table_file = args[2]
reference_file = args[3]
output_file = args[4]

# table cleaning
table <- read.table(table_file)
colnames(table) <- c("chromosome", "start", "end", "gene_id", "shape_index")

# overlap cleaning
overlap <- read.table(overlap_file)
overlap <- subset(overlap, select=-c(V1, V2, V3))
overlap <- subset(overlap, V9 != 0)
overlap <- overlap %>%
  distinct()
colnames(overlap) <- c("gene_name","chromosome","start","end","gene_id","overlapping")
overlap <- merge(overlap, table[, c("start","shape_index")], by = "start")
overlap <- overlap %>%
  group_by(gene_name) %>%
  filter(overlapping == max(overlapping)) %>%
  ungroup()
overlap <- overlap %>%
  group_by(gene_id) %>%
  filter(overlapping == max(overlapping)) %>%
  ungroup()
overlap <- subset(overlap, select=-c(overlapping,gene_id))

# addition of the cluster
reference <- read.table(reference_file)
reference <- rownames_to_column(reference, var = "gene_name")
overlap <- merge(overlap, reference[, c("gene_name","Label")], by = "gene_name", all.x = TRUE)
overlap$Label[is.na(overlap$Label)] <- "other_genes"
overlap <- overlap %>% relocate(chromosome, .before = start)
overlap <- overlap %>% relocate(gene_name, .after = end)

write.table(overlap,
            file = output_file,
            sep = "\t",          
            quote = FALSE,       
            col.names = TRUE,   
            row.names = FALSE)

