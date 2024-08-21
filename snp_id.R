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
#STEP 2
##This is for merging the tables that reported
## alleles and rsID-s
##hg38 positions for the rsID-s was pulled with
##bcftools
#--------------------------------------------------
##Ajore et al., 2022 , n=816, 4 lost from bcftools
ajore = merge(ajore_hg38_out, ajore_input, by.x = "V3", by.y = "V1")
ajore <- ajore[, c("V1", "V2.x", "V2.y", "V3.y")]
ajore$V1 <- sub("chr", "", ajore$V1)
ajore$author <- "Ajore et al., 2022"
ajore$study_type <- "Fine-mapping"
ajore$tested_variant <- paste(ajore$V1, ajore$V2.x, paste(ajore$V2.y, ajore$V3.y, sep = "/"), sep = "_")
names(ajore) <- c("chrom", "pos", "ref", "alt", "author", "study_type", "tested_variant")

#Castaldi et al.,2019 , n=296
castaldi = merge(castaldi_hg38_out, castaldi_input, by.x = "V3", by.y = "V1")
castaldi = castaldi[, c("V1", "V2.x", "V2.y", "V3.y")]
castaldi$V1 = sub("chr", "", castaldi$V1)
castaldi$author = "Castaldi et al., 2019"
castaldi$study_type = "Fine-mapping"
castaldi$tested_variant = paste(castaldi$V1, castaldi$V2.x, paste(castaldi$V2.y, castaldi$V3.y, sep = "/"), sep = "_")
names(castaldi) = c("chrom", "pos", "ref", "alt", "author", "study_type", "tested_variant")

#Joslin et al., 2021
joslin = merge(joslin_hg38_out, joslin_input, by.x = "V3", by.y = "V1")
joslin = joslin[, c("V1", "V2.x", "V2.y", "V3.y")]
joslin$V1 = sub("chr", "", joslin$V1)
joslin$author = "joslin et al., 2021"
joslin$study_type = "Fine-mapping"
joslin$tested_variant = paste(joslin$V1, joslin$V2.x, paste(joslin$V2.y, joslin$V3.y, sep = "/"), sep = "_")
names(joslin) = c("chrom", "pos", "ref", "alt", "author", "study_type", "tested_variant")

#Klein, et al., 2019
klein = merge(klein_hg38_out, klein_input, by.x = "V3", by.y = "V1")
klein = klein[, c("V1", "V2.x", "V2.y", "V3.y")]
klein$V1 = sub("chr", "", klein$V1)
klein$author = "Klein et al., 2019"
klein$study_type = "Fine-mapping"
klein$tested_variant = paste(klein$V1, klein$V2.x, paste(klein$V2.y, klein$V3.y, sep = "/"), sep = "_")
names(klein) = c("chrom", "pos", "ref", "alt", "author", "study_type", "tested_variant")

#Soemedi et al. 2017
soemedi = merge(soemedi_hg38_out, soemedi_input, by.x = "V3", by.y = "V1")
soemedi = soemedi[, c("V1", "V2.x", "V2.y", "V3.y")]
soemedi$V1 = sub("chr", "", soemedi$V1)
soemedi$author = "Soemedi et al., 2019"
soemedi$study_type = "Fine-mapping"
soemedi$tested_variant = paste(soemedi$V1, soemedi$V2.x, paste(soemedi$V2.y, soemedi$V3.y, sep = "/"), sep = "_")
names(soemedi) = c("chrom", "pos", "ref", "alt", "author", "study_type", "tested_variant")

#Madan et al., 2019
madan = madan_hg38_out %>%
  mutate(V1 = gsub("^chr", "", V1)) %>%
  select(-V5) %>%
  mutate(V4 = sapply(strsplit(V4, ","), `[`, 1)) %>%
  mutate(author = "Madan et al., 2019") %>%
  mutate(study_type = "Fine-mapping") %>%
  mutate(tested_variant = paste(V1, V2, V3, sep = "_") %>%
           paste(., V4, sep = "/")) %>%
  rename(
    chrom = V1,
    pos = V2,
    ref = V3,
    alt = V4
  )

#Liu et al., 2017
liu = liu_hg38_out %>%
  mutate(V1 = gsub("^chr", "", V1)) %>%
  select(-V5) %>%
  mutate(V4 = sapply(strsplit(V4, ","), `[`, 1)) %>%
  mutate(author = "Liu et al.,2017") %>%
  mutate(study_type = "Fine-mapping") %>%
  mutate(tested_variant = paste(V1, V2, V3, sep = "_") %>%
           paste(., V4, sep = "/")) %>%
  rename(
    chrom = V1,
    pos = V2,
    ref = V3,
    alt = V4
  )


#Tewhey et al., 2016
tewhey = tewhey_hg38_out %>%
  mutate(V1 = gsub("^chr", "", V1)) %>%
  select(-V5) %>%
  mutate(V4 = sapply(strsplit(V4, ","), `[`, 1)) %>%
  mutate(author = "Tewhey et al., 2016") %>%
  mutate(study_type = "Fine-mapping") %>%
  mutate(tested_variant = paste(V1, V2, V3, sep = "_") %>%
           paste(., V4, sep = "/")) %>%
  rename(
    chrom = V1,
    pos = V2,
    ref = V3,
    alt = V4
  )

#Bhattarai et al., 2023
bhattarai = bhattarai_hg38_out %>%
  mutate(V1 = gsub("^chr", "", V1)) %>%
  select(-V5) %>%
  mutate(V4 = sapply(strsplit(V4, ","), `[`, 1)) %>%
  mutate(author = "Bhattarai et al., 2023") %>%
  mutate(study_type = "Fine-mapping") %>%
  mutate(tested_variant = paste(V1, V2, V3, sep = "_") %>%
           paste(., V4, sep = "/")) %>%
  rename(
    chrom = V1,
    pos = V2,
    ref = V3,
    alt = V4
  )


bcf_out_combined = bind_rows(ajore, bhattarai, castaldi, joslin,
                             klein, liu, madan, soemedi, tewhey)

all_variants = rbind(tested_variants_all, bcf_out_combined)

unique_variants <- all_variants %>%
  distinct(tested_variant, .keep_all = TRUE)


write.table(unique_variants, file = "tested_variats_all_hg38.tsv", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
