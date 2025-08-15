library(readxl)
library(tidyverse)
library(rstatix)
library(tibble)
library(dplyr)
library(ggResidpanel)
args = commandArgs(trailingOnly=TRUE)

file = args[1]
output = args[2]

cluster <- read.table(file)
colnames(cluster) <- c('chromosome','start','end','gene_name')
cluster$length <- cluster$end - cluster$start

cluster <- subset(cluster, length > 10)
cluster <- subset(cluster, select=-length)

write.table(cluster,
       file = output,
       sep = "\t",
       quote = FALSE,
       col.names = FALSE,
       row.names = FALSE)