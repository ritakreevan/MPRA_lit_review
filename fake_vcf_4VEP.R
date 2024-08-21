library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(GenomciRanges)
library(biomartr)
library(motifbreakR)
library(cowplot)
library(readr)

#---------------------------
##STEP 3
##This is for making VEP
##suitable table
#--------------------------
unique_variants$tested_variant = NULL

vcf_table <- data.frame(
  `#CHROM` = unique_variants$chrom,
  POS = unique_variants$pos,
  ID = ".",
  REF = unique_variants$ref,
  ALT = unique_variants$alt,
  QUAL = ".",
  FILTER = "PASS",
  INFO = ".",
  FORMAT = "GT",
  FUNCTIONAL = "1/1",
  stringsAsFactors = FALSE
)

write.table(vcf_table, file = "tested_variats_vcf.tsv", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)