# Jordi Sevilla

# Cargar librerias ####
library(tidyverse)

# Leer variables de snakemake ####
tsv <- snakemake@input[["tsv"]]
output <- snakemake@output[["filtered_tsv"]]

# Leer el tsv
data <- read_tsv(tsv)

#filtrar por p-valor 

data <- filter(data, as.logical(PASS))

# filtrar por más de dos lecturas por cadena

is_deletion <- str_detect(data$ALT, "^[A-Z]", negate = T)
inBothStrands <- data$ALT_RV > 2 & data$ALT_DP > 2 & (data$ALT_RV + data$ALT_DP) >= 20
filter <- is_deletion | inBothStrands

data <- filter(data, filter)


#Añadir algunas variables de interés

data <- mutate(data,
                synonimous = case_when(REF_AA == ALT_AA ~ "yes",
                                        T ~"No"))



write_tsv(data,output)