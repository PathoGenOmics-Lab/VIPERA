# Jordi Sevilla

# Cargar librerias ####
library(tidyverse)

# Leer variables de snakemake ####
vcf <- snakemake@input[["vcf"]]
output <- snakemake@output[["vaf_file"]]

# Leer el archivo 
data <- read_delim(vcf, comment = "##")

# Nuevo df donde calcular VAF
new_df <- transmute(data, ID = `#CHROM`,
                         DP4 = str_extract(INFO, "DP4=([^;]*)" ,group = 1)) %>%
  rowwise() %>%
  mutate(VAF = (as.numeric(str_split(DP4,",")[[1]][1]) + as.numeric(str_split(DP4,",")[[1]][2]))/sum(as.numeric(str_split(DP4,",")[[1]])) ) %>%
  ungroup()

# introducir VAF en el df original 

data["VAF"] <- new_df$VAF

# Añadir otras variables interesantes de la anotación
data <- data %>%
    rowwise() %>%
    mutate(gen = str_extract(INFO,"Gene=([^;]*)", group = 1),
           IsSynonymus = str_extract(INFO,"IsSynonymous=([^;]*)", group = 1),
           COV = gsub(":\\.","",COV))

write_tsv(data,output)
