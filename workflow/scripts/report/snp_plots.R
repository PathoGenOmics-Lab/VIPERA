# LIBRERIAS #######

library(tidyverse)
library(stringi)
library(ggrepel)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")


# DISEÑO DE PLOTS ####
source(snakemake@params[["design"]])

# DATOS ####
vcf <- read_delim(snakemake@input[["vcf"]])
data <- read_csv(snakemake@params[["metadata"]]) %>%
  filter(ID %in% vcf$REGION) %>%
  dplyr::select(ID, CollectionDate)

# ANÁLISIS ####

# orden temporal 
date_order <- read_csv(snakemake@params[["metadata"]]) %>%
arrange(CollectionDate) %>%
pull(ID) %>%
unique()

# Modificar datos de variantes
vcf <- vcf %>% 
  mutate(GFF_FEATURE = gsub(":.*","",GFF_FEATURE),
         SNP = case_when(!is.na(REF_AA) ~ paste(GFF_FEATURE,":",REF_AA,POS_AA,ALT_AA, sep = ""),
                         T ~ paste(REF,POS,ALT, sep = ""))) %>%
                         unique() %>%
  dplyr::select(SNP,REGION,ALT_FREQ, POS)

IDs <- pull(vcf,REGION) %>%
  unique()

vcf <- vcf %>%
  pivot_wider(names_from = REGION, values_from = ALT_FREQ, values_fill = 0) %>% # Obtener los 0 en los puntos sin variantes
  pivot_longer(IDs, names_to = "REGION", values_to = "ALT_FREQ") %>%
  rowwise() %>%
  ungroup() 

# Unir datos
vcf <- left_join(vcf,data, by = c("REGION" = "ID"))

# Tiempo en dias entre muestras
vcf <- arrange(vcf,CollectionDate) %>%
  mutate(interval = as.numeric(CollectionDate - min(CollectionDate)))

# FIGURAS ####
## VOLCANO PLOT ####

# Cálculo de la correlación en el tiempo para cada SNP
SNPs <- pull(vcf,SNP) %>%
  unique()

cor.df <- data.frame(snp = "", cor = 0, p.value = 0) %>%
  filter(p.value != 0)

for (SNP_new in SNPs){
  df <- filter(vcf, SNP == SNP_new)
  
 test <- cor.test(df$ALT_FREQ,df$interval)
 pvalue <- p.adjust(test$p.value, method = "BH", n = length(SNPs))
 cor.df <- add_row(cor.df, snp = SNP_new, cor = test$estimate, p.value = pvalue)
}

volcano <- cor.df %>%
  mutate(trans.p = -log10(p.value),
         label = case_when(p.value < 0.05 ~ snp)) %>%
  ggplot() + 
  aes(x = cor, y = trans.p) + 
  geom_point() + 
  geom_label_repel(aes(label = label)) +
  xlim(c(-1,1)) +  
  geom_hline(aes(yintercept = -log10(0.05) ), linetype = 2, color = "red") + 
  labs(x = "Correlation", y = "-log10(p-value)")

ggsave(filename = snakemake@output[["pseudovolcano"]], 
        plot = volcano, width=159.2, 
        height=119.4, units="mm", 
        dpi=250)


## SNP PANEL ####

sign <- filter(cor.df, p.value < 0.05) %>% # SNPs significativamente correlacionados
    pull(snp) %>%
    unique()

dup <- vcf %>% # SNPs que comparten posición
  select(SNP,POS) %>%
  unique() %>%
  group_by(POS) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  pull(SNP) %>% 
  unique()

subset <- c(sign,dup) %>%
  unique()

plot.height = ceiling(length(subset)/4)*42 # Altura del plot para que se vea bién

panel <- vcf %>%
        filter(SNP %in% subset) %>%
        ggplot(aes(x = interval, y = ALT_FREQ, color = SNP)) + 
        scale_color_viridis_d() + 
        geom_point() + 
        geom_line() + 
        theme(legend.position = "bottom",
        legend.text = element_text(size = 9)) + 
        labs(x = "Days since first sample",
             y = "Frequency",
             color = "NV") + 
             guides(color = guide_legend(ncol = 4))

if (length(subset) > 1) {
  panel <- panel + 
    facet_wrap(
      vars(POS),
      nrow = ceiling(length(subset)/4),
      ncol = 4
    )
}

ggsave(filename = snakemake@output[["snp_panel"]], 
        plot = panel, width=159.2, 
        height = max(100, plot.height), units = "mm",
        dpi=250)


print("saving tables")

cor.df %>% 
transmute(
  NV = snp,
  PearsonCor = cor,
  adj_pvalue = p.value
) %>% 
write.csv(snakemake@output[["table_1"]], row.names = F)

vcf %>%
  filter(SNP %in% subset) %>%
  transmute(
    sample = REGION,
    POS = POS,
    NV = SNP,
    ALT_FREQ = ALT_FREQ,
    DaysSinceFirst = interval
    ) %>% 
    write.csv(snakemake@output[["table_2"]], row.names = F)

