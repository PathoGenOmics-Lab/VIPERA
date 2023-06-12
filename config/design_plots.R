# Jordi Sevilla

# Cargar librerias ####
library(pacman)
p_load("tidyverse",
       "showtext")

# Ajustes ####

# Tema

font_add_google("Montserrat", "Montserrat")
font_files()    
font_families()
showtext_auto()

theme_set(theme_minimal())

theme_update(text = element_text(size = 16,  family="Montserrat"), 
             axis.title = element_text(size = 16), 
             axis.line = element_line(size = 0.5, colour = "grey40", linetype=1, arrow = arrow(length = unit(3,"mm"))),
             panel.grid =element_line(size = 0.17, color = "lightgray")
) 




# AJUSTES ####

lineage_colors = c(B.1 = "#5CC7B2",
                   B.1.165 = "#DB3A34",
                   B.1.239 = "#E6913A",
                   B.1.238= "#FBD966",
                   B.1.399 = "#935FCA",
                   B.1.400 = "#D5A6BD",
                   Other = "gray39")

gene_colors = c(M = "#B4D4B4",
                N = "#B7B7B8",
                orf1ab = "#9CC4DC",
                ORF3a = "#ECB4B7",
                ORF8 = "#996D2B",
                S = "#F5CC9E",
                E = "#B2E1EA",
                ORF6 = "#F0D474",
                ORF7a = "#AA88CB",
                ORF10 = "#CACB5D")

snp_colors = c("C-26029-A" = "#e69138",
               "G-21985-T" = "#ffd966",
               "C-29421-T" = "#935fca",
               "C-23043-T" = "#7c1d6f",
               "C-28854-T" = "#ff4639"
               ,"C-5178-T" = "#d5a6bd"
               ,"C-4230-A" = "#1e3f9d"
               ,"C-16375-T" = "#1EA0AE"
               , "A-11923-C" = "#64acee"
               ,"G-12761-A" = "#084c61"
               ,"C-6033-T" = "chartreuse4"
               ,"C-4230-T" = "#93c47d",
               "Other" = "gray40")

