library(readxl)
library(tidyverse)
library(rstatix)
library(ggResidpanel)
library(tibble)
library(dplyr)
library(BSgenome)
library(BiocManager)
args = commandArgs(trailingOnly=TRUE)

plength = as.numeric(args[1])
genes_file = args[2]
studied_file = args[3]
table_file = args[4]
conversion_file = args[5]
genes_output = args[6]
ucsc_output = args[7]
table_to_liftOver = args[8]

### selection of first and last TSS per gene
genes_promoters = read.table(genes_file, sep=',', header=TRUE)
genes_promoters <- distinct(genes_promoters)
colnames(genes_promoters) <- c("gene_name","Ensembl_chrs","gene_start","gene_end","strand","TSS")
tss_min_max <- data.frame(genes_promoters %>%
                            group_by(gene_name) %>%
                            summarise(Minimum = min(TSS), Maximum = max(TSS)))
genes_promoters <- subset(genes_promoters, select=-c(TSS, gene_start, gene_end))
genes_promoters <- merge(genes_promoters, tss_min_max[, c("gene_name", "Minimum","Maximum")], by = "gene_name")
genes_promoters <- distinct(genes_promoters)
conversion <- read.table(conversion_file)
colnames(conversion) <- c("Ensembl_chrs", "chromosome", "length")
genes_promoters <- merge(genes_promoters, conversion[, c("Ensembl_chrs", "chromosome")], by = "Ensembl_chrs")
genes_promoters <- subset(genes_promoters, select=-Ensembl_chrs)

### obtaining promoter coordinates and check start and end positions
genes_promoters$Maximum <- as.numeric(genes_promoters$Maximum)
genes_promoters$Minimum <- as.numeric(genes_promoters$Minimum)
for (i in 1:nrow(genes_promoters)){
  if(genes_promoters$strand[i] == 1){genes_promoters$coord[i] <- genes_promoters$Maximum[i] - plength}
  else{genes_promoters$coord[i] <- genes_promoters$Minimum[i] + plength}}

genes_pstrand <- subset(genes_promoters, strand != "-1")
genes_pstrand <- genes_pstrand %>%
  mutate(start = pmin(coord, Maximum), end = pmax(coord, Maximum))
genes_pstrand <- subset(genes_pstrand, select=-c(Minimum, Maximum, coord))

genes_nstrand <- subset(genes_promoters, strand != "1")
genes_nstrand <- genes_nstrand %>%
  mutate(start = pmin(coord, Minimum), end = pmax(coord, Minimum))
genes_nstrand <- subset(genes_nstrand, select=-c(Minimum, Maximum, coord))

genes_pstrand$start[genes_pstrand$start < 1] <- 1
genome <- getBSgenome("BSgenome.Dmelanogaster.UCSC.dm6")
chrs_lengths <- seqlengths(genome)
genes_nstrand$end <- pmin(genes_nstrand$end, chrs_lengths[genes_nstrand$chromosome])

promoters <- rbind(genes_nstrand, genes_pstrand)
promoters <- promoters %>%
  mutate(start = pmin(start, end), end = pmax(start, end))
promoters <- subset(promoters, select=-strand)
promoters <- promoters %>%
  relocate(gene_name, .after = end)

### table cleaning
table <- read.table(table_file)
table <- subset(table, select=-c(V2,V3,V6,V7,V8))
table <- table %>%
  separate(col = V9, into = c("gene_name", "gene_id", "Window","total_reads","shape_index"), sep = ";")
table <- subset(table, select=-c(gene_name, Window, total_reads))
table <- table %>%
  separate(col = gene_id, into = c("delete","geneid"), sep = "=")
table <- table %>%
  separate(col = shape_index, into = c("delete","shape_index"), sep = "=")
table <- subset(table, select=-delete)
colnames(table) <- c("chromosome","start","end","gene_name", "shape_index")

write.table(table,
            file = table_to_liftOver,
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,   
            row.names = FALSE)

table <- subset(table, select=-shape_index)
table <- distinct(table)

### creation of coordinates bed files
write.table(promoters,
            file = genes_output,
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)

write.table(table,
            file = ucsc_output,
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)