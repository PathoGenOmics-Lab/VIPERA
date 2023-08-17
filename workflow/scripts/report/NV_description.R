# LIBRERIAS #######
library(tidyverse)
library(stringi)
library(ggpubr)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

# DISEÑO DE PLOTS ####
source(snakemake@params[["design"]])



# Anotacion de sars ####

SCov2_annotation = list(
  "five_prime_UTR"  = c(1:265),
  "orf1ab"          = c(266:21555),
  "Intergenic_1"    = c(21556:21562),
  "S"               = c(21563:25384),
  "Intergenic_2"    = c(25385:25392),
  "ORF3a"           = c(25393:26220),
  "Intergenic_3"    = c(26221:26244),
  "E"               = c(26245:26472),
  "Intergenic_4"    = c(26473:26522),
  "M"               = c(26523:27191),
  "Intergenic_5"    = c(27192:27201),
  "ORF6"            = c(27202:27387),
  "Intergenic_6"    = c(27388:27393),
  "ORF7"            = c(27394:27759),
  "Intergenic_7"    = c(27760:27893),
  "ORF8"            = c(27894:28259),
  "Intergenic_8"    = c(28260:28273),
  "N"               = c(28274:29533),
  "Intergenic_9"    = c(29534:29557),
  "ORF10"           = c(29558:29674),
  "three_prime_UTR" = c(29675:29903))


# DATOS ####
vcf <- read_delim(snakemake@input[["vcf"]])
vcf_snp <- vcf
window <- read_csv(snakemake@input[["window"]])

# ANÁLISIS ####

# orden temporal 
date_order <- read_csv(snakemake@params[["metadata"]]) %>%
arrange(CollectionDate) %>%
pull(ID) %>%
unique()


# Modificar datos de variantes
vcf <- vcf %>% 
  mutate(SNP = paste(REF,POS,ALT, sep = "-")) %>% # SNP
  dplyr::select(SNP,REGION,ALT_FREQ, GFF_FEATURE, synonimous) %>%
  rowwise() %>%
  mutate(POS = strsplit(SNP,"-")[[1]][2]) %>%
  ungroup()


# Df con las longitudes de los genes 
notation = data.frame(gene = "", len = 0) %>%
  filter(len != 0)
for (name in names(SCov2_annotation)){
  notation <- notation %>%
    add_row(gene = name, len = length(SCov2_annotation[[name]]))
}


# Clasificación de las variantes
vcf <- vcf %>%
  mutate(NV_class = case_when(str_detect(SNP,fixed("--")) | str_detect(SNP,fixed("+")) ~ "INDEL", # Diferenciar entre INDELS y SNP
                         T ~ "SNP"),
         Class = case_when(is.na(GFF_FEATURE) ~ "Intergenic", # Clasificación para los SNPs
                           T ~ synonimous),
         POS = as.numeric(POS)) %>%
  rowwise() %>%
  mutate( gene = as.character(window[window$position == POS,"gen"]), # Anotación
          indel_len = case_when(NV_class == "INDEL" & str_detect(SNP,fixed("--")) ~ str_length(strsplit(SNP,"--")[[1]][2]) -1, # Longitud de los INDELS
                                NV_class == "INDEL" & str_detect(SNP,fixed("-+")) ~ str_length(strsplit(SNP,"-+")[[1]][2]) -1),
          indel_class = case_when(gene == "Intergenic" ~ "Intergenic",  # Clasificación de los indels
                                  NV_class == "INDEL" & indel_len %% 3 == 0 ~ "In frame",
                                  NV_class == "INDEL" & indel_len %% 3 > 0 ~ "Frameshift")) %>%
  ungroup() %>%
  mutate(group = case_when(gene == "Intergenic" ~ "Intergenic", # Clasificación general
                           NV_class == "SNP" ~ Class,
                           NV_class == "INDEL" ~ indel_class))


# datos de nsps en ORF1ab
npc <- read_csv(snakemake@params[["nsp"]]) %>% 
            mutate(summary_nsp = case_when(NSP %in% paste("nsp",seq(4,12,1),sep = "") ~ "nsp4-12",
                                          NSP %in% paste("nsp",seq(14,16,1),sep = "") ~ "nsp14-16",
                                          T ~ NSP),
                  summaary_start =case_when(NSP %in% paste("nsp",seq(4,12,1), sep = "") ~ 8555,
                                        NSP %in% paste("nsp",seq(14,16,1),sep = "") ~ 18040,
                                        T ~ POS_i),
                  summaary_end =case_when(NSP %in% paste("nsp",seq(4,12,1), sep = "") ~ 16236,
                                      NSP %in% paste("nsp",seq(14,16,1),sep = "") ~ 21552,
                                      T ~ POS_f)) %>%
  filter(NSP != "nsp1")


# FIGURAS ####
# plot variantes en el tiempo 

variants <- vcf %>%
  filter(ALT_FREQ > 0) %>%
ggplot() + 
  aes(x = POS, y = factor(REGION,date_order), shape = factor(NV_class,c("SNP","INDEL")), color = group, alpha = ALT_FREQ) +
  geom_point(size = 3) + 
  geom_col(data = notation, aes(x = len,y = 0.3, fill = factor(gene,rev(names(SCov2_annotation)))), inherit.aes = F, width = 0.3) +
  scale_fill_manual(values = gene_colors) + 
  xlim(c(0,29903)) + 
  scale_color_manual(labels = NV_names, 
                    values  = NV_colors) + 
  labs(x = "SARS-CoV-2 genome position", y = "Sample", shape = "Variant class", color = "Classification", alpha = "Frequency", fill = "Region") 

# porcentaje por ventanas

window_plot <- ggplot(window) + 
  aes(x = position, y = fractions, color = gen) + 
  geom_point() +
  geom_line(aes(group = 1), colour = "black", alpha = 0.3) +
  scale_y_continuous(label = scales::percent, limits = c(0,max(window$fractions) + 0.005)) + 
  xlim(c(0,29903)) + 
  scale_color_manual(values = gene_colors) +
  labs(y = "Proportion of \n sites with SNV", x = "", color = "Gen")


window_plot_nsp <- window_plot + 
  geom_vline(data = npc, aes(xintercept = summaary_start), color = "red") + 
  geom_vline(data = npc, aes(xintercept = summaary_end), color = "red") + 
  geom_text(data = npc, aes(x = (summaary_start + summaary_end)/2, y = max(window$fractions) + 0.002, label = summary_nsp), inherit.aes = F, size = 5, angle = 60)


figura <- ggarrange(window_plot_nsp,
          variants,nrow = 3,
          align = "v" ,
          legend.grob = get_legend(variants)
          , heights = c(2,6), 
          legend = "right",
          labels = c("A","B"))


ggsave(filename = snakemake@output[["fig"]], 
        plot = figura, width=250, 
        height=240, units="mm", 
        dpi=250)


# Figura SNPs con el tiempo

figur_SNP_time <- vcf_snp %>%
                  filter(ALT_FREQ <= 0.95) %>%
                  left_join(read_csv(snakemake@params[["metadata"]]), by = c("REGION" = "ID")) %>%
                  group_by(REGION) %>%
                  summarise(CollectionDate = min(as.Date(CollectionDate)),
                            n = n()) %>%
                  ungroup() %>%
                  ggplot() + 
                  aes(x = CollectionDate, y = n) +
                  geom_smooth(method = "lm",fill = "gray95", alpha = 0.6) + 
                  geom_point() +
                  stat_cor(geom = "label") + 
                  labs(x = "Date", y = "Nº of polimorphic sites")

ggsave(filename = snakemake@output[["fig_cor"]], 
        plot = figur_SNP_time, width=250, 
        height=119.4, units="mm", 
        dpi=250)

# TABLA ####

n_indels <- filter(vcf, NV_class == "INDEL") %>% length()
n_snv <- length(unique(vcf$SNP)) - n_indels

df <- data.frame(nv = c("SNP","INDEL"),n = c(n_snv,n_indels))
write.csv(df, snakemake@output[["summary_nv"]], row.names = F)

# Zoom in in spike
spike.pos <- 
  window %>% 
  filter(gen == "S") %>% 
  pull(position)

print(spike.pos)

window_plot_spike <- window %>% 
      filter(gen == "S") %>% 
      ggplot() + 
      aes(x = position, y = fractions, color = gen) + 
      geom_point() +
      geom_line(aes(group = 1), colour = "black", alpha = 0.3) +
      scale_y_continuous(label = scales::percent, limits = c(0,max(window$fractions) + 0.005)) + 
      xlim(c(min(spike.pos),max(spike.pos))) + 
      scale_color_manual(values = gene_colors, guide = "none") +
      labs(y = "Proportion of \n sites with SNV", x = "")

variants_spike <- vcf %>%
  filter(ALT_FREQ > 0,
        POS %in% spike.pos) %>%
  ggplot() + 
    aes(x = POS, y = factor(REGION,date_order), shape = factor(NV_class,c("SNP","INDEL")), color = group, alpha = ALT_FREQ) +
    geom_point(size = 3) +  
    xlim(c(min(spike.pos),max(spike.pos))) + 
    scale_color_manual(labels = NV_names, 
                    values  = NV_colors) +
    labs(x = "SARS-CoV-2 genome position", 
        y = "Sample", 
        shape = "Variant class", 
        color = "Classification", 
        alpha = "Frequency", 
        fill = "Region") 


figura_spike <- ggarrange(window_plot_spike,
                variants_spike,nrow = 3,
                align = "v" ,
                legend.grob = get_legend(variants)
                , heights = c(2,6), 
                legend = "right",
                labels = c("A","B"))

ggsave(filename = snakemake@output[["fig_s"]], 
        plot = figura_spike, width=250, 
        height=240, units="mm", 
        dpi=250)

print("Saving tables")

# figure 7b
vcf %>%
select(
  REGION,
  POS,
  SNP,
  ALT_FREQ,
  NV_class,
  group
) %>% 
rename(
  sample = REGION,
  NV = SNP,
  Class = group
) %>% filter(ALT_FREQ >0) %>%
mutate(
  Class = case_when(
    Class == "yes" ~ "synonymous",
    Class == "No" ~ "non_synonymous",
    T ~ Class
    )
  ) %>% 
write.csv(snakemake@output[["table_2"]], row.names = FALSE)

# figure 7a

window %>% 
transmute(
  POS = position,
  feature = gen,
  prop_PolimorphicSites = fractions
) %>% 
write.csv(snakemake@output[["table_1"]], row.names = FALSE)

vcf_snp %>%
  filter(ALT_FREQ <= 0.95) %>%
  left_join(read_csv(snakemake@params[["metadata"]]), by = c("REGION" = "ID")) %>%
  group_by(REGION) %>%
  summarise(CollectionDate = min(as.Date(CollectionDate)),
            n_PolimorphicSites = n()) %>%
  ungroup() %>% 
  rename(sample = REGION) %>% 
  write.csv(snakemake@output[["table_3"]], row.names = FALSE)