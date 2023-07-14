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

# FIGURA ####

# Anotacion de sars

SCov2_annotation = list(
  # "SCov2_genome"    = (1, 29903),
  "five_prime_UTR"  = c(    1:   265),
  "orf1ab"          = c(  266: 21555),
  "Intergenic_1" = c(21556:21562),
  "S"          = c(21563: 25384),
  "Intergenic_2" = c(25385:25392),
  "ORF3a"           = c(25393: 26220),
  "Intergenic_3" = c(26221:26244),
  "E"          = c(26245: 26472),
  "Intergenic_4" = c(26473:26522),
  "M"          = (26523: 27191),
  "Intergenic_5" = c(27192:27201),
  "ORF6"            = c(27202: 27387),
  "Intergenic_6" = c(27388:27393),
  "ORF7a"           = c(27394: 27759),
  "Intergenic_7" = c(27760:27893),
  "ORF8"            = c(27894: 28259),
  "Intergenic_8" = c(28260:28273),
  "N"          = c(28274: 29533),
  "Intergenic_9" = c(29534:29557),
  "ORF10"           = c(29558: 29674),
  "three_prime_UTR" = c(29675: 29903))

dic = data.frame(pos = c(1:29903))

gene = c()
for (pos in c(1:29903)){
  for (name in names(SCov2_annotation)){
    annotated = F
    if (pos %in% SCov2_annotation[[name]] & !annotated){
      gene <- c(gene,name)
      annotated = T
      break
      
    } 
    
  }
  if (!annotated ){
    gene = c(gene,"Intergenic")
  }
}

dic["gene"] <- gene

dic <- mutate(dic, gene = case_when(str_detect(gene,"prime") | str_detect(gene,"Inter") ~ "Intergenic",
                                    T ~ gene))
# Df con las longitudes de los genes 

notation = data.frame(gene = "", len = 0) %>%
  filter(len != 0)
for (name in names(SCov2_annotation)){
  notation <- notation %>%
    add_row(gene = name, len = length(SCov2_annotation[[name]]))
}
# Clasificación de las variantes

vcf <- vcf %>%
  mutate(SNP_class = case_when(str_detect(SNP,fixed("--")) | str_detect(SNP,fixed("+")) ~ "INDEL",
                         T ~ "SNP"),
         Class = case_when(SNP_class == "INDEL" ~ "Intergenic",
                           is.na(GFF_FEATURE) ~ "Intergenic",
                           T ~ synonimous),
         POS = as.numeric(POS)) %>%
  rowwise() %>%
  mutate( gene = as.character(dic[dic$pos == POS,"gene"]),
          indel_len = case_when(SNP_class == "INDEL" & str_detect(SNP,fixed("--")) ~ str_length(strsplit(SNP,"--")[[1]][2]) -1,
                                SNP_class == "INDEL" & str_detect(SNP,fixed("-+")) ~ str_length(strsplit(SNP,"-+")[[1]][2]) -1),
          indel_class = case_when(gene == "Intergenic" ~ "Intergenic",
                                  SNP_class == "INDEL" & indel_len %% 3 == 0 ~ "In frame",
                                  SNP_class == "INDEL" & indel_len %% 3 > 0 ~ "Frameshift")) %>%
  ungroup() %>%
  mutate(group = case_when(gene == "Intergenic" ~ "Intergenic",
                           SNP_class == "SNP" ~ Class,
                           SNP_class == "INDEL" ~ indel_class))
plot <- vcf %>%
  filter(ALT_FREQ > 0) %>%
ggplot() + 
  aes(x = POS, y = factor(REGION,date_order), shape = factor(SNP_class,c("SNP","INDEL")), color = group, alpha = ALT_FREQ) +
  geom_point(size = 3) + 
  geom_col(data = notation, aes(x = len,y = 0.3, fill = factor(gene,rev(names(SCov2_annotation)))), inherit.aes = F, width = 0.3) +
  scale_fill_manual(values = gene_colors) + 
  xlim(c(0,29903)) + 
  scale_color_manual(labels = c("Frameshift","Inframe","Intergenic","Non synonymous","Synonymous"), values = c("#568D63","black","#B27CF9","#AE584A","#0248FD")) + 
  labs(x = "SARS-CoV-2 genome position", y = "Sample", shape = "Variant class", color = "Classification", alpha = "Frequency", fill = "Region") 

plot

# FIGURA CON LAS VENTANAS ####

window <- read.csv("i2sysbio/Case-study-SARS-CoV-2/output/case_study.window.csv")
window["group"] <- 1

library(scales)
npc <- read_csv("i2sysbio/Case-study-SARS-CoV-2/data/report_files/nsp_annotation.csv")
npc <- mutate(npc,summary_nsp = case_when(NSP %in% paste("nsp",seq(4,12,1),sep = "") ~ "nsp4-12",
                                          NSP %in% paste("nsp",seq(14,16,1),sep = "") ~ "nsp14-16",
                                          T ~ NSP),
              summaary_start =case_when(NSP %in% paste("nsp",seq(4,12,1), sep = "") ~ 8555,
                                        NSP %in% paste("nsp",seq(14,16,1),sep = "") ~ 18040,
                                        T ~ POS_i),
              summaary_end =case_when(NSP %in% paste("nsp",seq(4,12,1), sep = "") ~ 16236,
                                      NSP %in% paste("nsp",seq(14,16,1),sep = "") ~ 21552,
                                      T ~ POS_f)) %>%
  filter(NSP != "nsp1")
window_plot <- ggplot(window) + 
  aes(x = position, y = fractions, color = gen) + 
  geom_point() +
  geom_line(aes(group = group), colour = "black", alpha = 0.3) +
  scale_y_continuous(label = scales::percent, limits = c(0,0.012)) + 
  xlim(c(0,29903)) + 
  scale_color_manual(values = gene_colors) +
  labs(y = "Proportion of \n sites with SNV", x = "", color = "Gen")  +
  theme()
window_plot

window_plot_nsp <- window_plot + 
  geom_vline(data = npc, aes(xintercept = summaary_start), color = "red") + 
  geom_vline(data = npc, aes(xintercept = summaary_end), color = "red") + 
  geom_label(data = npc, aes(x = (summaary_start + summaary_end)/2, y = 0.011, label = summary_nsp), inherit.aes = F, size = 5)
window_plot_nsp



ggarrange(window_plot_nsp,
          plot,nrow = 3,
          align = "hv" ,
          legend.grob = get_legend(plot)
          , heights = c(1,3), 
          legend = "right")
