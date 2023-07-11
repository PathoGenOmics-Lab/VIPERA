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

# DATOS ####
matrix <- read_csv(snakemake@input[["dist"]])
metadata <- read_csv(snakemake@params[["metadata"]])

# nj tree

tree <- matrix %>%
        column_to_rownames(var = "...1") %>%
        as.matrix() %>%
        as.dist() %>%
        nj() %>%
        root("NC_045512.2", resolve.root = TRUE)

fix_negative_edge_length <- function(nj.tree) {
  edge_infos <- cbind(nj.tree$edge, nj.tree$edge.length) %>% as.data.table
  colnames(edge_infos) <- c('from', 'to', 'length')
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

tree <- fix_negative_edge_length(tree)

tempest <- adephylo::distRoot(tree, "all", method = "patristic") %>% as.data.frame() %>% 
                    rownames_to_column(var = "ID") %>%
                    filter(ID != "NC_045512.2" ) %>%
                    rename(distance = ".") %>% 
                     mutate(distance = distance) %>%
                    left_join(select(metadata, ID,CollectionDate)) %>%
                     mutate(date_interval = as.numeric(as.Date(CollectionDate) - min(as.Date(CollectionDate))))



# PLOTS ####
# tree
tree_plot <- ggtree(tree) + 
geom_tiplab() + 
geom_treescale() + 
geom_rootedge(0.5)

ggsave(filename = snakemake@output[["tree"]], 
        plot = tree_plot, width=159.2, 
        height=119.4, units="mm", 
        dpi=250)


# temp-est

tempest_fig <- ggplot(tempest) +
                aes(x = date_interval, y = distance) + 
                geom_smooth(method = "lm",fill = "gray95", alpha = 0.6) +
                stat_cor(geom = "label") +
                geom_point() +
                labs(y = "Root to tip distance", x = "Time since first sample")


ggsave(filename = snakemake@output[["temest"]], 
        plot = tempest_fig, width=159.2, 
        height=119.4, units="mm", 
        dpi=250)



# ML tree para el contexto ####

tree <- read.tree(snakemake@input[["ml"]]) %>%
  root("NC_045512.2", resolve.root = TRUE)

colors <- c()
for(ID in tree$tip.label){
  if(ID %in% metadata$ID){
    colors <- c(colors,"red")
  } else{
    colors <- c(colors,"gray80")
  }
}

plot <- ggtree(tree) + 
  geom_tippoint(color = colors) + 
  geom_treescale() + 
  geom_rootedge(0.00001)


ggsave(filename = snakemake@output[["tree_ml"]], 
        plot = plot, width=159.2, 
        height=119.4, units="mm", 
        dpi=250)