# Jordi Sevilla

# LIBRERIAS ####
library(tidyverse)

# DATOS ####

metadata <- read.csv(snakemake@params[["metadata"]])

pango_report <- read.csv(snakemake@input[["report"]])

# SUMMARY ####

metadata <- select(metadata,ID,CollectionDate) %>%
            filter(ID %in% pango_report$taxon) %>%
            left_join(select(pango_report,taxon,lineage), by = c("ID" = "taxon")) %>% 
            rename(Sample = ID,
            Collection_Date = CollectionDate,
            Lineage = lineage)

write.csv(metadata, snakemake@output[["table"]], row.names = F)
