library(tidyverse)
library(tibble)
library(dplyr)
library(Biostrings)
library(BSgenome)
args = commandArgs(trailingOnly=TRUE)

plength = as.numeric(args[1])
genes_file = args[2]
conversion_file =args[3]
output_file = args[4]

### slection of first & last TSS per gene
genes_promoters = read.table(genes_file, sep='\t', header=TRUE)
genes_promoters <- distinct(genes_promoters)
colnames(genes_promoters) <- c("gene_name","Ensembl_chrs","gene_start","gene_end","strand","TSS")
tss_min_max <- data.frame(genes_promoters %>%
                            group_by(gene_name) %>%
                            summarise(Minimum = min(TSS), Maximum = max(TSS)))
genes_promoters <- subset(genes_promoters, select=-c(TSS, gene_start, gene_end))
genes_promoters <- merge(genes_promoters, tss_min_max[, c("gene_name", "Minimum","Maximum")], by = "gene_name", all.y = TRUE)
genes_promoters <- distinct(genes_promoters)
conversion <- read.table(conversion_file)
colnames(conversion) <- c("Ensembl_chrs", "chromosome")
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
genome <- getBSgenome("BSgenome.Drerio.UCSC.danRer11")
chrs_lengths <- seqlengths(genome)
genes_nstrand$end <- pmin(genes_nstrand$end, chrs_lengths[genes_nstrand$chromosome])

promoters <- rbind(genes_nstrand, genes_pstrand)
promoters <- promoters %>%
  mutate(start = pmin(start, end), end = pmax(start, end))
promoters <- subset(promoters, select=-strand)
promoters <- promoters %>%
  relocate(gene_name, .after = end)

### creation of coordinates bed file
write.table(promoters,
            file = output_file,
            sep = "\t",          
            quote = FALSE,       
            col.names = FALSE,   
            row.names = FALSE)
