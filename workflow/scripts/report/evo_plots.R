# LIBRERIAS #######
library(tidyverse)


# DISEÑO DE PLOTS ####
source(snakemake@params[["design"]])

# DATOS ####
vcf <- read_delim(snakemake@input[["vcf"]])

vcf <- vcf %>% 
  mutate(SNP = paste(REF,POS,ALT, sep = "-")) %>%
  dplyr::select(SNP,REGION,ALT_FREQ, GFF_FEATURE, synonimous) %>%
  rowwise() %>%
  mutate(POS = strsplit(SNP,"-")[[1]][2]) %>%
  ungroup() 

metadata <- read_delim(snakemake@params[["metadata"]]) %>%
  mutate(interval = as.numeric(as.Date(CollectionDate) - min(as.Date(CollectionDate)))) %>%
  select(ID,interval) %>%
  rename(REGION = ID)

vcf <- left_join(vcf,metadata)


N_S_position <- read_delim(snakemake@input[["N_S"]])

# PLOT

plot <- vcf %>%
              filter(ALT_FREQ > 0.05) %>%
              group_by(REGION,synonimous) %>%
              summarise(Freq = sum(ALT_FREQ, na.rm = T)) %>%
              pivot_wider( names_from = synonimous, values_from = Freq, values_fill = 0 )  %>%
              transmute(dn = No/sum(N_S_position$S),
                        ds = yes/sum(N_S_position$S)) %>%
                        ungroup() %>%
              pivot_longer(c("dn","ds"), values_to = "value", names_to = "d") %>%
              left_join(unique(select(vcf,REGION,interval))) %>%
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