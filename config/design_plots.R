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
GENE_PALETTE <- c(
  M = "#B4D4B4",
  N = "#B7B7B8",
  ORF1ab = "#9CC4DC",
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
TREE_PALETTE <- c(
  tip_label = "#D944AA99",
  boot_alrt_pass = "#64ACEEB2"
)

TREE_NODE_SIZE <- c(
  tip_label = 2,
  boot_alrt_pass = 0.8
)

TREE_LEGEND_NAMES <- c(
  tip_label = "Target samples",
  boot_alrt_pass = sprintf(
    "UFBoot ≥ %s%s & SH-aLRT ≥ %s%s ",
    snakemake@params[["boot_th"]],
    "%",
    snakemake@params[["alrt_th"]],
    "%"
  )
)

# Nucleotide variants classification colors and labels
NV_TYPE_PALETTE <- c(
  Frameshift = "#568D63",
  "In frame" = "black",
  Intergenic = "#B27CF9",
  No = "#AE584A",
  Yes = "#0248FD"
)

NV_TYPE_NAMES <- c(
  Frameshift = "Frameshift",
  "In frame" = "Inframe",
  Intergenic = "Intergenic",
  No = "Non synonymous",
  Yes = "Synonymous"
)

# dn ds colors and labels
DNDS_LABELS <- c(
  dn = "dN",
  ds = "dS"
)

DNDS_COLORS <- c(
  dN = "#E53E47",
  dS = "#2C47F5"
)

DNDS_SHAPES <- c(
  dN = 2,
  dS = 4
)

# Allele frequency trajectories panel color
ALL.COLORS <- grDevices::colors()
TRAJECTORY.PANEL.COLORS <- ALL.COLORS[
  !grepl("(gray|grey|white|snow|azure|beige)", ALL.COLORS)
]
