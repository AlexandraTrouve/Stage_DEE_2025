library(readxl)
library(tidyverse)
library(rstatix)
library(ggResidpanel)
library(tibble)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

input_files = args[1:6]
modules_file = args[7]
output_file =args[8]

all_data_list <- list()
for (files in input_files) {
  df <- read.csv(files, header = TRUE)
  all_data_list[[length(all_data_list) + 1]] <- df}
final_df <- do.call(rbind, all_data_list)
final_df <- subset(final_df, select = c(ID, AlphaPurif, p_value_null_purif, p_value_purif_pos, p_value_null_pos, Conclusion))
final_df$corrected_p_value_null_purif <- p.adjust(final_df$p_value_null_purif, method = "fdr")
final_df$corrected_p_value_purif_pos <- p.adjust(final_df$p_value_purif_pos, method = "fdr")
final_df$corrected_p_value_null_pos <- p.adjust(final_df$p_value_null_pos, method = "fdr")
final_df <- subset(final_df, select = -c(p_value_null_purif, p_value_purif_pos, p_value_null_pos))

final_df <- final_df %>%
  mutate(corrected_conclusion = case_when(
    corrected_p_value_purif_pos < 0.01 & corrected_p_value_null_pos ~ "Directional (+)",
    corrected_p_value_null_purif < 0.01 & corrected_p_value_purif_pos < 0.01 & corrected_p_value_null_pos ~ "Neutral model",
    corrected_p_value_null_purif < 0.01 & AlphaPurif > 1 ~ "Stabilizing",
    corrected_p_value_null_purif < 0.01 & AlphaPurif < 1 ~ "Disruptive",
    TRUE ~ Conclusion
  ))

final_df <- final_df %>%
  separate(ID, into = c("chromosome", "start", "end", "gene_name"), sep = ":")
final_df <- subset(final_df, select=-c(chromosome, start, end))
modules <- read.table(modules_file, header = TRUE)
modules <- rownames_to_column(modules, var = "gene_name")
final_df <- merge(final_df, modules[, c("gene_name", "Label")], by = "gene_name", all.x = TRUE)
final_df$Label[is.na(final_df$Label)] <- "other_genes"

df_prop1 <- final_df %>%
  group_by(corrected_conclusion, Label) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Label) %>%
  mutate(proportion = count / sum(count))

pdf(file = output_file, width = 10, height = 9)
ggplot(data = df_prop1, aes(x = Label, y= proportion, fill = corrected_conclusion)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(proportion, 2)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 3) +
  ylim(0, 1) +
  scale_fill_manual(values = c("Directional (+)" = "lightcoral", "Disruptive" = "darkgoldenrod1",
                               "Neutral model" = "darkgrey", "Stabilizing" = "mediumaquamarine",
                               "Directional (-)" = "lightsteelblue")) +
  labs(title = "Proportion of selection pressure in terms of Cluster",
       x = "Cluster",
       y = "Proportion of promoters") +
  theme_minimal()

final_df <- final_df[!(final_df$Label %in% "other_genes"),]
df_prop <- final_df %>%
  group_by(corrected_conclusion, Label) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(corrected_conclusion) %>%
  mutate(proportion = count / sum(count))
ggplot(data = df_prop, aes(x = corrected_conclusion, y= proportion, fill = Label)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(proportion, 2)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 3) +
  ylim(0, 0.35) +
  scale_fill_manual(values = c("hourglass" = "pink", "inverse_hourglass" = "lightsalmon",
                               "early_constraint" = "lightblue", "late_constraint" = "thistle2",
                               "uniform_high" = "darkseagreen3", "uniform_low" = "snow3")) +
  labs(title = "Proportion of promoters form each cluster in terms of selection pressure",
       x = "Selection pressure",
       y = "Proportion of promoters") +
  theme_minimal()

dev.off()

