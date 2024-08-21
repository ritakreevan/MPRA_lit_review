library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(GenomciRanges)
library(biomartr)
library(motifbreakR)
library(cowplot)
library(readr)

#--------------------------------------------------
#STEP 1
##This is for reading in positions tables
##Position tables are from those studies
##that report chromosome-position-tested alleles
#--------------------------------------------------
folder_path = "/Users/rita/Desktop/MPRA_meta_analysis/Collected_variants/Positions"

# Get a list of all txt files in the folder, so I dont have to do 1-by-1
file_paths = list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)

col_names = c("chrom", "pos", "ref", "alt", "author", "study_type")

process_file = function(file_path) {
  df <- read_delim(file_path, delim = "\t", col_names = FALSE, show_col_types = FALSE)
  print(paste("Processing file:", file_path))
  print(head(df))  # Print first few rows to check data
  print(dim(df))   # Print dimensions to check number of columns
  colnames(df) <- col_names
  df <- df %>%
    mutate(chrom = as.character(chrom)) %>%
    mutate(chrom = ifelse(chrom == "X", "23", chrom))

  df <- df %>%
    mutate(across(c(pos, ref, alt), as.character))  
  df <- df %>%
    mutate(tested_variant = paste(chrom, pos, ref, alt, sep = "_")) %>%
    mutate(tested_variant = gsub("_", "/", tested_variant, fixed = TRUE))
  df_unique <- df %>%
    distinct(tested_variant, .keep_all = TRUE)
  
  return(df_unique)
}

list_of_dfs <- map(file_paths, ~ process_file(.x))

# Combine all dataframes into one master table with all variants
tested_variants_all <- bind_rows(list_of_dfs) %>%
  distinct(tested_variant, .keep_all = TRUE) %>%
  arrange(tested_variant)

#How many lines total in all the tables
total_lines <- map_int(list_of_dfs, nrow) %>% sum()

#how many unique values in "tested_variants" column
unique_count <- tested_variants_all %>%
  pull(tested_variant) %>%
  n_distinct()
