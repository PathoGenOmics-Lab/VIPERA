# Jordi Sevilla
# LIBRERIAS ####
source("i2sysbio/Case-study-SARS-CoV-2/config/design_plots.R")


library(pacman)
p_load("tidyverse",
       "jsonlite",
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
       "treeio")
library(data.table)


# DATA####

vcf <- read_delim("i2sysbio/Case-study-SARS-CoV-2/output/case_study.masked.filtered.tsv")


# Eliminar malas combinaciones 

samples <- pull(vcf,REGION) %>% unique()
positions <- pull(vcf,POS) %>% unique()
depth <- data.frame(pos = sort(positions))

for (sample in samples) {
  
  path <- sprintf("i2sysbio/Case-study-SARS-CoV-2/output/demixing/%s/%s_depth.txt",sample,sample)
  data <- read_delim(path, col_names = c("SAMPLE","position","Base","Depth"))
  data <- filter(data, position %in% positions)
  depth[sprintf("%s",sample)] <- data$Depth
}

depth <- pivot_longer(depth, contains("COV"), names_to = "COV", values_to = "depth") %>%
  filter(depth < 30)

bad_combinations <- paste(depth$COV,depth$pos, sep = "-")

#
vcf <- vcf %>% 
  mutate(SNP = paste(REF,POS,ALT, sep = "-")) %>%
  dplyr::select(SNP,REGION,ALT_FREQ, GFF_FEATURE, synonimous)

vcf <- vcf %>%
  pivot_wider(names_from = REGION, values_from = ALT_FREQ, values_fill = 0) %>%
  pivot_longer(contains("COV"), names_to = "REGION", values_to = "ALT_FREQ") %>%
  rowwise() %>%
  mutate(POS = strsplit(SNP,"-")[[1]][2]) %>%
  ungroup() %>%
  filter(!(paste(REGION,POS,sep="-") %in% bad_combinations))

data <- read_csv("i2sysbio/Case-study-SARS-CoV-2/data/report_files/JSF_COV_metadata_master.csv") %>%
  filter(ID %in% vcf$REGION) %>%
  dplyr::select(ID, CollectionDate)

vcf <- left_join(vcf,data, by = c("REGION" = "ID"))

vcf <- arrange(vcf,CollectionDate)


# Orden temporal de las muestras
date_order <- dplyr::select(vcf, REGION, CollectionDate) %>%
  unique() %>%
  pull(REGION)

d <- table(vcf$SNP, vcf$REGION) %>% as.data.frame() %>%
  group_by(Var1) %>%
  arrange(factor(Var2,date_order), .by_group = T) %>%
  mutate(p = case_when(str_detect(paste(Freq, collapse = ""),"1000*1") ~ "Saltatory",
                       T ~ "Continuous")) %>%
  ungroup() %>%
  filter(p == "Saltatory") %>%
  pull(Var1) %>%
  unique()

vcf <- vcf %>%
  mutate(interval = as.numeric(CollectionDate - min(CollectionDate))) %>%
  group_by(SNP) %>%
  mutate(n = sum(ALT_FREQ > 0)) %>%
  ungroup() %>%
  group_by(SNP) %>%
  mutate(p = cor(ALT_FREQ,interval)) %>%
  arrange(factor(REGION,date_order), .by_group = T) %>%
  ungroup()


vcf <- vcf %>%
  mutate(case = case_when(p > 0.48 ~ "Aumenta",
                          p < -0.3 ~ "Disminuye",
                          T ~ "Aleatorio"),
         number = case_when(SNP %in% d ~ "discontinuo",
                            n %in% c(1:2) ~ "[1-2]",
                            n %in% c(3:5) ~ "[3-5]",
                            n %in% c(6:7) ~ "[6-7]",
                            n %in% c(8:9) ~ "[8-9]",
                            n > 9 ~ "> 9"),
         GFF_FEATURE = gsub(":.*","",GFF_FEATURE))


# GRAFICA POSICION 4230 ####

transformer <- function(x){
  y = x*0.11894735*61977
  
  return(round(y))
}

# Datos de depth
c <- read_delim("i2sysbio/Case-study-SARS-CoV-2/output/case_study.masked.filtered.tsv") %>%
  filter(POS == 4230) %>%
  group_by(REGION) %>%
  mutate(Depth = sum(unique(ALT_DP)) + REF_DP + sum(unique(ALT_RV)) + REF_RV) %>%
  ungroup() %>%
  dplyr::select(REGION,Depth) %>%
  unique() %>%
  add_row(REGION = "COV012638", Depth = 2369) %>%
  mutate(pct = (Depth/sum(Depth))/0.12050053)

vcf %>% rowwise() %>%
  mutate(ALT = strsplit(SNP,"-")[[1]][3],
         REF = strsplit(SNP,"-")[[1]][1]) %>%
  ungroup() %>%
  group_by(REGION) %>%
  filter(POS == 4230) %>%
  transmute(ALT = ALT,
            REF = REF,
            REF_FREQ = 1 - sum(unique(ALT_FREQ)),
            ALT_FREQ = ALT_FREQ,
            interval = interval) %>%
  ungroup() %>%
  pivot_longer(contains("FREQ"), names_to = "type", values_to = "FREQ") %>%
  mutate(allel = case_when(type == "REF_FREQ" ~ paste(REF," (REF)"), 
                           type == "ALT_FREQ" ~ ALT),
         interval = factor(interval,sort(unique(interval)))) %>%
  dplyr::select(!c("ALT","REF")) %>%
  unique() %>%
  left_join(c) %>%
  ggplot(aes(x = interval, y = FREQ, fill = allel)) + 
  geom_col() + 
  geom_point(aes(y = pct,)) + 
  scale_y_continuous(sec.axis = dup_axis(labels = transformer, name = "Depth")) + 
  labs(y = "Frecuencia", x = "Tiempo a la primera muestra (dias)", fill = "Alelo" )


# VERSIONES DE TENDENCIAS ####

## POR POSICION ####
vcf %>%
  filter(SNP %in% names(snp_colors)) %>%
  ggplot() + 
  aes(x = interval, y = ALT_FREQ, color = SNP) +  
  geom_point() + 
  scale_color_manual(values = snp_colors) +
  geom_line(aes(group = SNP)) + 
  facet_grid(case ~  number) + 
  scale_x_continuous(breaks = c(0,50,100,150,200), limits = c(0,250)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(vjust = 4,margin = margin(t = 0, r = -20, b = 0, l = 0)),
        legend.position = "bottom",
        plot.margin = margin(1,1,1,1, "cm")) + 
  labs(x = "", y = "Frequency") + 
  facet_wrap(vars(POS))

## POR TENDENCIA ####

vcf %>%
  filter(SNP %in% names(snp_colors)) %>%
  ggplot() + 
  aes(x = interval, y = ALT_FREQ, color = SNP) +  
  geom_point() + 
  scale_color_manual(values = snp_colors) +
  geom_line(aes(group = SNP)) + 
  facet_grid(case ~  number) + 
  scale_x_continuous(breaks = c(0,50,100,150,200), limits = c(0,250)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.y = element_text(vjust = 4)) + 
  labs(x = "", y = "Frequency") + 
  facet_wrap(vars(case))



# PCA ####


pca.1<- vcf %>%
  filter(synonimous == "No", number != "[1-2]") %>%
  select(SNP,REGION,ALT_FREQ, case) %>%
  pivot_wider(names_from = REGION,values_from = ALT_FREQ, values_fill = 0) %>%
  column_to_rownames(var = "SNP")

res <- princomp(pca.1[c(2:13)])

autoplot(res, data = pca.1, color = "case") + 
  scale_color_hue(labels = c("Disminuyen", "Aumentan","Random"))
