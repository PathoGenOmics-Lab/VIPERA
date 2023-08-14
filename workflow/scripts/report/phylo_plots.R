library(tidyverse)
library(stringi)
library(ape)
library(ggtree)
library(data.table)
library(ggpubr)
library(pegas)

# Write stdout and stderr to log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output")

# DISEÑO DE PLOTS ####
source(snakemake@params[["design"]])

# DATOS ####
matrix <- read_csv(snakemake@input[["dist"]])
metadata <- read_csv(snakemake@params[["metadata"]])
tree_ml <- read.tree(snakemake@input[["ml"]]) %>%
  root(snakemake@params[["ref_name"]], resolve.root = TRUE)

study_names <- read.dna(
  snakemake@input[["study_fasta"]],
  format = "fasta",
  as.matrix = FALSE,
) %>%
  names()

date_order <- read_csv(snakemake@params[["metadata"]]) %>%
  arrange(CollectionDate) %>%
  filter(ID %in% study_names) %>%
  pull(ID) %>%
  unique()
  tree_tiplab <- data.frame(
    ID = date_order,
    order = seq(1, length(date_order), 1)
  ) %>%
  rowwise() %>%
  mutate(tip_label = sprintf("(%s)-%s", order, ID)) %>%
  ungroup() %>%
  add_row(ID = snakemake@params[["ref_name"]],
          order = 0,
          tip_label = snakemake@params[["ref_name"]])


# ANÁLISIS ####

# Con n-j puede darse que hayan ramas con longitud negativa si la matriz de distancias no es congruente
# Esta función arregla estas distancias negativas en el caso de que ocurran
fix_negative_edge_length <- function(nj.tree) {
  edge_infos <- cbind(nj.tree$edge, nj.tree$edge.length) %>% as.data.table
  colnames(edge_infos) <- c("from", "to", "length")
  nega_froms <- edge_infos[length < 0, sort(unique(from))]
  nega_froms
  for (nega_from in nega_froms) {
    minus_length <- edge_infos[from == nega_from, ][order(length)][1, length]
    edge_infos[from == nega_from, length := length - minus_length]
    edge_infos[to == nega_from, length := length + minus_length]
  }
  nj.tree$edge.length <- edge_infos$length
  nj.tree
}

## Construcción de Arboles ####
tree <- matrix %>%
        column_to_rownames(var = "...1") %>%
        as.matrix() %>%
        as.dist() %>%
        nj() %>%
        root(snakemake@params[["ref_name"]], resolve.root = TRUE)

tree <- fix_negative_edge_length(tree)

## Calculo de distancia partistica al ancestro ####
tempest <- adephylo::distRoot(tree, "all", method = "patristic") %>% as.data.frame() %>% 
                    rownames_to_column(var = "ID") %>%
                    filter(ID != snakemake@params[["ref_name"]] ) %>%
                    rename(distance = ".") %>% 
                     mutate(distance = distance) %>%
                    left_join(select(metadata, ID,CollectionDate)) %>%
                     mutate(date_interval = as.numeric(as.Date(CollectionDate) - min(as.Date(CollectionDate))))

### Datos modelo lineal
model <- lm(distance ~ date_interval, data = tempest)

df <- data.frame(
  stat = c("sub_rate", "r2", "pvalue"),
  values = c(
    model$coefficients[[2]],
    summary(model)$r.squared[[1]],
    cor.test(tempest$distance, tempest$date_interval)$p.value
  )
)

df <- column_to_rownames(df, var = "stat")
write.csv(df, snakemake@output[["stats_lm"]], row.names = FALSE)

# FIGURAS ####
# Árbol por distancias ponderadas
tree_plot <- ggtree(tree) %<+% tree_tiplab +
  geom_tiplab(aes(label = tip_label)) +
  geom_treescale() +
  geom_rootedge(0.5) +
  xlim(0, max(tempest$distance) + 3.5)

ggsave(
  filename = snakemake@output[["tree"]],
  plot = tree_plot,
  width = snakemake@params[["plot_width_mm"]],
  height = snakemake@params[["plot_height_mm"]],
  units = "mm",
  dpi = 250
)

# Análisis temp-est
tempest_fig <- ggplot(tempest, aes(x = date_interval, y = distance)) +
  geom_smooth(method = "lm", fill = "gray95", alpha = 0.6, color = "red") +
  geom_point() +
  labs(y = "Root to tip distance", x = "Time since first sample")

ggsave(filename = snakemake@output[["temest"]],
  plot = tempest_fig,
  width = snakemake@params[["plot_width_mm"]],
  height = snakemake@params[["plot_height_mm"]],
  units = "mm",
  dpi = 250
)

# ML tree para el contexto
tip.color <- ifelse(tree_ml$tip.label %in% study_names, "blue", NA)

# Node labels contain SH-aLRT/UFboot values
aLRT.values <- sapply(
  strsplit(tree_ml$node.label, "/"),
  function(x) {
    as.numeric(x[1])
  }
)
bootstrap.values <- sapply(
  strsplit(tree_ml$node.label, "/"),
  function(x) {
    as.numeric(x[2])
  }
)
aLRT.mask <- aLRT.values >= snakemake@params[["alrt_th"]]
boot.mask <- bootstrap.values >= snakemake@params[["boot_th"]]

node.color <- case_when(
  aLRT.mask & boot.mask ~ "red",
  !aLRT.mask & boot.mask ~ "#ff6600",
  aLRT.mask & !boot.mask ~ "#ffbf51"
)

study.mrca <- getMRCA(tree_ml, study_names)
p <- ggtree(tree_ml, layout = "circular") +
          geom_highlight(node = study.mrca, colour = "red", fill = "red", alpha = 0.2) +
          geom_tippoint(color = tip.color, shape = 1) +
          geom_treescale(x = 0.0008) +
          geom_rootedge(0.0005) +
          xlim(-0.0008, NA) +
          geom_nodepoint(color = node.color, shape = 20)

ggsave(
  filename = snakemake@output[["tree_ml"]],
  plot = p,
  width = snakemake@params[["plot_width_mm"]],
  height = snakemake@params[["plot_height_mm"]],
  units = "mm",
  dpi = 250
)
