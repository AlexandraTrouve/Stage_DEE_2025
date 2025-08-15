# libraries importation
library(tidyverse)
library(tibble)
library(dplyr)
library(tools)
args = commandArgs(trailingOnly=TRUE)

input_files <- args[1:2]
output_file <- args[3]

all_data <- list()
for (file_path in input_files) {
  cluster <- tools::file_path_sans_ext(basename(file_path))
  df <- read.table(file_path)
  colnames(df) <- c("Kmer", paste0(cluster))
  all_data[[cluster]] <- df}

df_names <- names(all_data)
final_df <- all_data[[df_names[1]]]
for (i in 2:length(all_data)) {
  final_df <- full_join(final_df, all_data[[df_names[i]]], by = "Kmer")}

write.csv(final_df,
          output_file,
          row.names = FALSE)

