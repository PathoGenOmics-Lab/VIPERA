# LIBRERIAS #######
library(tidyverse)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")


# DISEÑO DE PLOTS ####
source(snakemake@params[["design"]])

# DATOS ####
vcf <- read_delim(snakemake@input[["vcf"]])
metadata <- read_delim(snakemake@params[["metadata"]])
N_S_position <- read_delim(snakemake@input[["N_S"]])

# ANÁLISIS ####
# Modificacion de datos de variantes
vcf <- vcf %>% 
  mutate(SNP = paste(REF,POS,ALT, sep = "-")) %>% # SNP al que da lugar cada variante
  dplyr::select(SNP,REGION,ALT_FREQ, GFF_FEATURE, synonimous) %>%
  rowwise() %>%
  mutate(POS = strsplit(SNP,"-")[[1]][2]) %>% # Volver a obtener la posición del SNP
  ungroup() 

# Modificar metadatos
metadata <- metadata %>%
  mutate(interval = as.numeric(as.Date(CollectionDate) - min(as.Date(CollectionDate)))) %>% # Tiempo en dias entre muestras
  select(ID,interval) %>%
  rename(REGION = ID)

# Combinar bases de datos
vcf <- left_join(vcf,metadata)


# FIGURA ####

plot <- vcf %>%
              group_by(REGION,synonimous) %>%
              summarise(Freq = sum(ALT_FREQ, na.rm = T)) %>% # Cácular la frecuencia de mutaciones sinónimas y no sinónimas por genoma
              pivot_wider( names_from = synonimous, values_from = Freq, values_fill = 0 )  %>%
              transmute(dn = No/sum(N_S_position$S), # Mutaciones no sinónimas por sitio no sinónimo
                        ds = yes/sum(N_S_position$S)) %>% # Mutaciones sinónimas por sitio no sinónimo
                        ungroup() %>%
              pivot_longer(c("dn","ds"), values_to = "value", names_to = "d") %>%
              left_join(unique(select(vcf,REGION,interval))) %>% # Recuperar info de tiempo
              ggplot() + 
              aes(x = interval, y = value, color = d) + 
              geom_point() +
              geom_line() + 
              scale_color_hue(labels = c("dN","dS")) + 
              labs(y = "", x = "Time since first sample", color = "Parameter")

ggsave(filename = snakemake@output[["plot"]], 
        plot = plot, width=159.2, 
        height=119.4, units="mm", 
        dpi=250)