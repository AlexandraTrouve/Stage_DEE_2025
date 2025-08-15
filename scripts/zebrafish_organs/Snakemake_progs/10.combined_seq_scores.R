# libraries importation
library(tidyverse)
library(tibble)
library(dplyr)
library(tools)
args = commandArgs(trailingOnly=TRUE)

input_files <- args[1:8]
output_file <- args[9]

all_data_list <- list()
for (file_path in input_files) {
  df <- read.csv(file_path, header = TRUE, sep = ',')
  colnames(df) <- c("gene_name", "SVM_score")
  all_data_list[[length(all_data_list) + 1]] <- df}

final_df <- do.call(rbind, all_data_list)

write.csv(final_df,
          output_file,
          row.names = FALSE)

