#!/usr/bin/env python3
import pandas as pd

step = 1000
SCov2_annotation = {
# "SCov2_genome"    : [1, 29903],
"five_prime_UTR"  : [    1,   265],
"orf1ab"          : [  266, 21555],
"S"          : [21563, 25384],
"ORF3a"           : [25393, 26220],
"E"          : [26245, 26472],
"M"          : [26523, 27191],
"ORF6"            : [27202, 27387],
"ORF7a"           : [27394, 27759],
"ORF8"            : [27894, 28259],
"N"          : [28274, 29533],
"ORF10"           : [29558, 29674],
"three_prime_UTR" : [29675, 29903]
}


def cov_list(df):
    
    COV = list(df["REGION"].unique())
    
    return COV

def get_polimorphic_sites(df):
    
    lista = set(df.POS)
    
    
    return lista

def window_calculous(lista,step,genome_size, annotation):
    
    positions = []
    pct = []
    genes =[] 
    
    lim_sup = genome_size + 1
    
    for position in range(1,lim_sup):
        if position - step not in range(1,lim_sup):
            positions.append(position)
            pct.append(0)
            
            gene = 0
        
            for key,item in annotation.items():
                if position in range(item[0],item[1] + 1):
                    gene = key             

            if gene != 0:
                genes.append(gene)
            else: 
                genes.append("Intergenic")
            continue
        
        
        
        
        
        num_snp = 0
        
        for x in lista:
            if x in range(position - step,position + 1):
                
                num_snp += 1
                
                
        gene = 0
        
        for key,item in annotation.items():
            if position in range(item[0],item[1] + 1):
                gene = key             
        
        if gene != 0:
            genes.append(gene)
        else: 
            genes.append("Intergenic")
            
        frac = num_snp/1000
        
        positions.append(position)
        pct.append(frac)
    
    df = {"position":positions, "fractions":pct, "gen":genes}
    
    df = pd.DataFrame(df)
    
    return df
            
def main():
    df = pd.read_table(snakemake.input.vcf)    
    COV = cov_list(df)
    sites = get_polimorphic_sites(df)
    frame = window_calculous(sites,step,29903,SCov2_annotation)
    frame.to_csv(snakemake.output.window_df)


if __name__ == "__main__":
    main()