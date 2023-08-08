# Jordi Sevilla

# LIBRERIAS ####
library(tidyverse)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

# DATOS ####

metadata <- read.csv(snakemake@params[["metadata"]])

pango_report <- read.csv(snakemake@input[["report"]])

date_order <- read_csv(snakemake@params[["metadata"]]) %>%
arrange(CollectionDate) %>%
filter(ID %in% pango_report$taxon) %>%
pull(ID) %>%
unique()

# SUMMARY ####

index  <- data.frame(Sample = date_order, 
                    Index = seq(1,length(date_order),1))

metadata <- select(metadata,ID,CollectionDate) %>%
            filter(ID %in% pango_report$taxon) %>%
            left_join(select(pango_report,taxon,lineage), by = c("ID" = "taxon")) %>% 
            rename(Sample = ID,
            Collection_Date = CollectionDate,
            Lineage = lineage)

table <- left_join(index,metadata)

write.csv(table, snakemake@output[["table"]], row.names = F)
