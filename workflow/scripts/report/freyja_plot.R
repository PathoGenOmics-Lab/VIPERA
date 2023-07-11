# LIBRERIAS #######
library(tidyverse)

# DISEÑO DE PLOTS ####
source(snakemake@params[["design"]])

# DATOS ####

demix <- read_csv(snakemake@input[["summary_demixing"]])


# ANÁLISIS ####

# Muestras a destacar 
main_lineages <- demix %>%
  group_by(sample) %>%
  top_n(1,abundances) %>%
  ungroup() %>%
  pull(lineages) %>%
  unique()

# date order
date_order <- read_csv(snakemake@params[["metadata"]]) %>%
arrange(CollectionDate) %>%
filter(ID %in% demix$sample) %>%
pull(ID) %>%
unique()

# FIGURA #####

demix_plot <- demix %>%
  mutate(lin_2 = case_when(lineages %in% main_lineages ~ lineages,
                            T ~ "Other")) %>%
  ggplot() + 
  aes(x = factor(sample,date_order), y = as.numeric(abundances), fill = lin_2) + 
  scale_fill_viridis_d(option = "Mako") + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom") + 
  labs(x = "Sample", y = "Relative abundance", fill = "Lineage")


ggsave(filename = snakemake@output[["fig"]], 
        plot = demix_plot, width=159.2, 
        height=119.4, units="mm", 
        dpi=250)