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

theme_update(
  text = element_text(size = 16, family = "Montserrat"),
  axis.title = element_text(size = 16),
  axis.line = element_line(
    linewidth = 0.5,
    colour = "grey40",
    linetype = 1,
    arrow = arrow(length = unit(3, "mm"))
  ),
  panel.grid = element_line(linewidth = 0.17, color = "lightgray")
)

# Gene palette
gene_colors <- c(
  M = "#B4D4B4",
  N = "#B7B7B8",
  orf1ab = "#9CC4DC",
  ORF3a = "#ECB4B7",
  ORF8 = "#996D2B",
  S = "#F5CC9E",
  E = "#B2E1EA",
  ORF6 = "#F0D474",
  ORF7 = "#AA88CB",
  ORF10 = "#CACB5D"
)

# Nucleotide diversity
DIVERSITY_PALETTE <- c(
  density_fill = "#fcbf49",
  density_color = "#eae2b7",
  value_color = "#d62828",
  dnorm_color = "#f77f00"
)

# M-L tree colors and labels
tree_colors <- c(
  tip_label = "#D944AA99",
  boot_alrt_pass = "#64ACEEB2"
)

node.size <- c(
  tip_label = 2,
  boot_alrt_pass = 0.8
)

# Nucleotide variants classification colors and labels
NV_colors <- c(
  Frameshift = "#568D63",
  "In frame" = "black",
  Intergenic = "#B27CF9",
  No = "#AE584A",
  Yes = "#0248FD"
)

NV_names <- c(
  Frameshift = "Frameshift",
  "In frame" = "Inframe",
  Intergenic = "Intergenic",
  No = "Non synonymous",
  Yes = "Synonymous"
)

# dn ds colors and labels
dnds.labels <- c(
  dn = "dN",
  ds = "dS"
)

dnds.colors <- c(
  dn = "#E53E47",
  ds = "#2C47F5"
)

dnds.shapes <- c(
  dn = 2,
  ds = 4
)
