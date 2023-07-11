# LIBRERIAS #######

library(ape)
library(pegas)
library(future.apply)
library(tidyverse)


# DISEÑO DE PLOTS ####

source(snakemake@params[["design"]])

# DATOS ####
gene_ex <- read.dna(snakemake@params[["outgroup_aln"]],format = "fasta", as.matrix = F) # Genomas del contexto
study_aln <- read.dna(snakemake@input[["study_fasta"]],format = "fasta", as.matrix = F) # Genomas del estudio

gene_ex <- gene_ex[names(gene_ex) != "NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome"]

# Análisis ####

## Funciones ####

# Función para calcular pi a partir de una muestra
nucleotide.diversity <- function(dna_object, record.names, sample.size){
  sample <- sample(record.names, sample.size, replace = F)
  dna_subset <- dna_object[record.names %in% sample]
  return(nuc.div(dna_subset))
}

# función para hacer el bootstrapping en paralelo
plan(multisession, workers = 4)
boot.nd.parallel <- function(aln, sample.size = 12, reps = 100) {
  record.names <- names(aln)
  future_sapply(
    1:reps,
    function(x) nucleotide.diversity(aln, record.names, sample.size),
    future.seed = TRUE
  )
}

## Cálculo de PI


# valor de diversidad (PI) para nuestras muestras
study_aln <- study_aln[names(study_aln) != "NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome"]
diversity = nuc.div(study_aln)
write.table(data.frame(div = diversity), snakemake@output[["value"]], row.names = F) # Guardar el valor para el report


# Valores de diversidad para las muestras del contexto

divs <- boot.nd.parallel(gene_ex, length(names(study_aln)) -1, 1000)

# FIGURA ####

plot <- data.frame(pi = divs) %>%
  ggplot() + 
  aes(pi) + 
  geom_density(aes(x = pi) ,fill = "#2177d4", alpha = 0.7, bw = 0.000001,color = "white" ) + 
  geom_vline(aes(xintercept = diversity), color = "red") + 
  labs(x = "π" , y = "Density")


ggsave(filename = snakemake@output[["fig"]], 
        plot = plot, width=159.2, 
        height=119.4, units="mm", 
        dpi=250)
