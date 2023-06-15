#!/usr/bin/env python3
import pandas as pd


# TODO: use GFF file for annotation
# TODO: don't hard-code the genome length (29903)
SCov2_annotation = {
    # "SCov2_genome"    : [1, 29903],
    "five_prime_UTR"  : [    1,   265],
    "orf1ab"          : [  266, 21555],
    "S"               : [21563, 25384],
    "ORF3a"           : [25393, 26220],
    "E"               : [26245, 26472],
    "M"               : [26523, 27191],
    "ORF6"            : [27202, 27387],
    "ORF7a"           : [27394, 27759],
    "ORF8"            : [27894, 28259],
    "N"               : [28274, 29533],
    "ORF10"           : [29558, 29674],
    "three_prime_UTR" : [29675, 29903]
}


def get_polimorphic_sites(df):
    return set(df.POS)


def window_calculation(sites, step, genome_size, coord):
    positions = []
    pct = []
    genes = []
    lim_sup = genome_size + 1
    for position in range(1, lim_sup):
        # Add gene
        gene = coord[position]
        genes.append(gene)
        # Add percent (excluding initial and final positions)
        if position - step not in range(1, lim_sup):
            pct.append(0)
        else:
            # Calculate SNPs
            num_snp = 0
            for x in sites:
                if x in range(position - step, position + 1):
                    num_snp += 1
            pct.append(num_snp / step)
        # Add positions
        positions.append(position)
    return pd.DataFrame({"position": positions, "fractions": pct, "gen": genes})


def main():
    # Diccionario con relación coordenada --> gen
    coord2gene = {i : "intergenic" for i in range(1, 29903 + 1)}
    for gene in SCov2_annotation:
        start, end = SCov2_annotation[gene]
        for i in range(start, end + 1):
            coord2gene[i] = gene
    # Process
    df = pd.read_table(snakemake.input.vcf)    
    sites = get_polimorphic_sites(df)
    frame = window_calculation(sites, snakemake.params.step, 29903, coord2gene)
    frame.to_csv(snakemake.output.window_df)


if __name__ == "__main__":
    main()
