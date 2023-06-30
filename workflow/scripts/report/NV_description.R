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



# Anotacion de sars

SCov2_annotation = list(
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



# DATOS ####

# orden temporal 

date_order <- read_csv(snakemake@params[["metadata"]]) %>%
arrange(CollectionDate) %>%
pull(ID) %>%
unique()



vcf <- read_delim(snakemake@input[["vcf"]])

vcf <- vcf %>% 
  mutate(SNP = paste(REF,POS,ALT, sep = "-")) %>%
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


# datos de calculo de porcentaje por ventanas
window <- read_csv(snakemake@input[["window"]])

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


# PLOTS ####
# plot variantes en el tiempo 

variants <- vcf %>%
  filter(ALT_FREQ > 0) %>%
ggplot() + 
  aes(x = POS, y = factor(REGION,date_order), shape = factor(SNP_class,c("SNP","INDEL")), color = group, alpha = ALT_FREQ) +
  geom_point(size = 3) + 
  geom_col(data = notation, aes(x = len,y = 0.3, fill = factor(gene,rev(names(SCov2_annotation)))), inherit.aes = F, width = 0.3) +
  scale_fill_manual(values = gene_colors) + 
  xlim(c(0,29903)) + 
  scale_color_manual(labels = c("Frameshift","Inframe","Intergenic","Non synonymous","Synonymous"), values = c("#568D63","black","#B27CF9","#AE584A","#0248FD")) + 
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
  geom_label(data = npc, aes(x = (summaary_start + summaary_end)/2, y = max(window$fractions) + 0.002, label = summary_nsp), inherit.aes = F, size = 5)


figura <-   ggarrange(window_plot_nsp,
          variants,nrow = 3,
          align = "v" ,
          legend.grob = get_legend(variants)
          , heights = c(2,6), 
          legend = "right")


ggsave(filename = snakemake@output[["fig"]], 
        plot = figura, width=250, 
        height=240, units="mm", 
        dpi=250)

