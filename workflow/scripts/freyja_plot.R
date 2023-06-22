# LIBRERIAS #######
library(pacman)
p_load("tidyverse",
       "stringi",
       "flextable",
       "ggpubr",
       "ggtree",
       "ape",
       "adephylo",
       "plotly",
       "ggrepel",
       "apex",
       "adegenet",
       "pegas",
       "mmod",
       "poppr",
       "treeio",
       "data.table",
       "future.apply",
       "scales",
       "quarto",
       "showtext")

# DISEÑO DE PLOTS ####
source(snakemake@params[["design"]])

# datos

demix <- read_csv(snakemake@input[["summary_demixing"]])

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

# PLOT #####

demix_plot <- demix %>%
  mutate(lin_2 = case_when(lineages %in% main_lineages ~ lineages,
                            T ~ "Other")) %>%
  ggplot() + 
  aes(x = factor(sample,date_order), y = as.numeric(abundances), fill = lin_2) + 
  scale_fill_manual(values = lineage_colors) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "bottom") + 
  labs(x = "Sample", y = "Relative abundance", color = "Lineage")


ggsave(filename = snakemake@output[["fig"]], 
        plot = demix_plot, width=159.2, 
        height=119.4, units="mm", 
        dpi=250)