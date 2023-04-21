# Jordi Sevilla

# Cargar librerias ####
library(tidyverse)

# Leer variables de snakemake ####
directories <- snakemake@input[["directories"]]
output_path <- snakemake@params[["demix_folder"]]

# Resumir datos de freyja ####

# Crear dataframe vacío donde contener los datos
demix <- data.frame("lineages" = NA, "abundances" = NA, "sample" = NA) %>%
  filter(!is.na(sample))

# Para cada archivo extraer la info e introducirla en el df
for (directory in directories){

    COV <- str_match(directory,"COV......")

    path <- sprintf("%s/%s_demixed.tsv", directory, COV)

    data <- read_tsv(path, col_names = c("variable", "valor")) %>%
        filter(row_number() %in% c(3, 4)) %>%
        pivot_wider(names_from = variable, values_from = valor) %>%
        separate_rows(lineages, abundances, sep = " ")

    data["sample"] <- COV

    demix <- rbind(demix, data)
}

write.csv(demix, snakemake@output[["summary_df"]], row.names = FALSE)