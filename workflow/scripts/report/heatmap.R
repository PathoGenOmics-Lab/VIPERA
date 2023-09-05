# LIBRERIAS #######
library(tidyverse)

# DATOS

vcf <- read_tsv(snakemake@input[["vcf"]])
metadata <- read.csv(snakemake@params[["metadata"]])

# Orden temporal 

date_order <- metadata %>%
  arrange(CollectionDate) %>%
  pull(ID) %>%
  unique()

vcf <- vcf %>% 
  mutate(GFF_FEATURE = gsub(":.*","",GFF_FEATURE),
         SNP = case_when(!is.na(REF_AA) ~ paste(GFF_FEATURE,":",REF_AA,POS_AA,ALT_AA, sep = ""),
                         T ~ paste(REF,POS,ALT, sep = ""))) %>%
  unique() %>%
  dplyr::select(SNP,REGION,ALT_FREQ, POS)

vcf <- select(vcf,SNP,REGION,ALT_FREQ) %>% 
  pivot_wider(names_from = SNP, values_from = ALT_FREQ, values_fill = 0) %>%
  arrange(factor(REGION, levels = date_order)) %>%
  column_to_rownames(var = "REGION")

write.csv(vcf, snakemake@output[["table"]])
