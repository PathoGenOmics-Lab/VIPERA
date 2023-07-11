# LIBRERIAS #######

library(tidyverse)
library(stringi)
library(ggrepel)


# DISEÑO DE PLOTS ####
source(snakemake@params[["design"]])

# DATOS ####

# orden temporal 

date_order <- read_csv(snakemake@params[["metadata"]]) %>%
arrange(CollectionDate) %>%
pull(ID) %>%
unique()



vcf <- read_delim(snakemake@input[["vcf"]])

vcf <- vcf %>% 
  mutate(SNP = paste(REF,POS,ALT, sep = "-")) %>%
  dplyr::select(SNP,REGION,ALT_FREQ, GFF_FEATURE, synonimous)

IDs <- pull(vcf,REGION) %>%
  unique()
vcf <- vcf %>%
  pivot_wider(names_from = REGION, values_from = ALT_FREQ, values_fill = 0) %>%
  pivot_longer(IDs, names_to = "REGION", values_to = "ALT_FREQ") %>%
  rowwise() %>%
  mutate(POS = strsplit(SNP,"-")[[1]][2]) %>%
  ungroup() 

data <- read_csv(snakemake@params[["metadata"]]) %>%
  filter(ID %in% vcf$REGION) %>%
  dplyr::select(ID, CollectionDate)

vcf <- left_join(vcf,data, by = c("REGION" = "ID"))

vcf <- arrange(vcf,CollectionDate) %>%
  mutate(interval = as.numeric(CollectionDate - min(CollectionDate)))

## VOLCANO PLOT ####

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


# SNP PANEL ####
sign <- filter(cor.df, p.value < 0.05) %>%
    pull(snp) %>%
    unique()

dup <- vcf %>%
  select(SNP,POS) %>%
  unique() %>%
  group_by(POS) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  pull(SNP) %>% 
  unique()

subset <- c(sign,dup) %>%
  unique()

length = ceiling(length(subset)/4)*40

panel <- vcf %>%
        filter(SNP %in% subset) %>%
        ggplot(aes(x = interval, y = ALT_FREQ, color = SNP)) + 
        scale_color_viridis_d() + 
        geom_point() + 
        geom_line() + 
        facet_wrap(vars(POS), nrow = ceiling(length(subset)/4), ncol = 4) + 
        theme(legend.position = "bottom")

ggsave(filename = snakemake@output[["snp_panel"]], 
        plot = panel, width=159.2, 
        height=length, units="mm", 
        dpi=250)