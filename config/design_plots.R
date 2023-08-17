# Jordi Sevilla


library(tidyverse)
library(showtext)

# Ajustes ####
showtext_auto(enable = FALSE)
showtext_opts(dpi = 200)

# Tema
font_add_google("Montserrat", "Montserrat")
showtext_auto()

theme_set(theme_minimal())

theme_update(text = element_text(size = 16,  family = "Montserrat"),
             axis.title = element_text(size = 16),
             axis.line = element_line(size = 0.5, colour = "grey40", linetype=1, arrow = arrow(length = unit(3,"mm"))),
             panel.grid =element_line(size = 0.17, color = "lightgray")
) 




# AJUSTES ####


gene_colors = c(M = "#B4D4B4",
                N = "#B7B7B8",
                orf1ab = "#9CC4DC",
                ORF3a = "#ECB4B7",
                ORF8 = "#996D2B",
                S = "#F5CC9E",
                E = "#B2E1EA",
                ORF6 = "#F0D474",
                ORF7 = "#AA88CB",
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


tree_colors = c(
    tip_label = "blue",
    Bootstrap_pass = "#ff6600",
    alrt_pass = "#ffbf51",
    boot_alrt_pass = "red"
)