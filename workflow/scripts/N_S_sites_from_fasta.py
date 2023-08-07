#!/usr/bin/env python3

import logging
from Bio import SeqIO
from Bio import Seq
import copy
import pandas as pd

SCov2_annotation = {
# "SCov2_genome"    : [1, 29903],
#"five_prime_UTR"  : [    1,   265],
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
#"three_prime_UTR" : [29675, 29903]
}

def split_into_codons(seq):
    """Split the complete CDS feature in to a list of codons."""

    codons = [seq[i:i + 3] for i in range(0, len(seq), 3) if 'N' not in seq[i:i + 3] ]
    
    return codons

def geneticCode(name):

    """ Dictionary that maps codons to amino acids """ 

    if name == 'standard':
        gc = {  'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', \
                'AGT':'S','ATA':'I','ATC':'I','ATG':'M','ATT':'I','CAA':'Q','CAC':'H','CAG':'Q','CAT':'H','CCA':'P','CCC':'P','CCG':'P', \
                'CCT':'P','CGA':'R','CGC':'R','CGG':'R','CGT':'R','CTA':'L','CTC':'L','CTG':'L','CTT':'L','GAA':'E','GAC':'D','GAG':'E', \
                'GAT':'D','GCA':'A','GCC':'A','GCG':'A','GCT':'A','GGA':'G','GGC':'G','GGG':'G','GGT':'G','GTA':'V','GTC':'V','GTG':'V', \
                'GTT':'V','TAA':'*','TAC':'Y','TAG':'*','TAT':'Y','TCA':'S','TCC':'S','TCG':'S','TCT':'S','TGA':'*','TGC':'C','TGG':'W', \
                'TGT':'C','TTA':'L','TTC':'F','TTG':'L','TTT':'F'  }
    return gc

def potential_changes_dict(nt_to_aa):

    """ Generate a dictionary, with S and N pre-calculated for all 
    possible codons (key:codon, value: (S,N).
    ARGS:
        nt_to_aa, a dict mapping codons (keys), e.g. 'TTA', to 
            amino-acid letter (values), e.g. 'L'
            
            e.g. geneticCode("standard")
    Notes:
        Sources of formulae:
        http://www.megasoftware.net/mega4/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_nei_gojobori_method.htm
    """

    potential_changes = {   'S': {  'AAA':0.0,'AAC':0.0,'AAG':0.0,'AAT':0.0, 'ACA':0.0, 'ACC':0.0, 'ACG':0.0, 'ACT':0.0, 'AGA':0.0, 'AGC':0.0, \
                                    'AGG':0.0, 'AGT':0.0, 'ATA':0.0, 'ATC':0.0, 'ATG':0.0, 'ATT':0.0, 'CAA':0.0, 'CAC':0.0, 'CAG':0.0, 'CAT':0.0, \
                                    'CCA':0.0,'CCC':0.0,'CCG':0.0,'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0, \
                                    'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0,'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0, \
                                    'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0,'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0, \
                                    'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0,'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0},

                            'N': {  'AAA':0.0, 'AAC':0.0, 'AAG':0.0, 'AAT':0.0, 'ACA':0.0, 'ACC':0.0, 'ACG':0.0, 'ACT':0.0, 'AGA':0.0, 'AGC':0.0, 'AGG':0.0, \
                                    'AGT':0.0,'ATA':0.0,'ATC':0.0,'ATG':0.0,'ATT':0.0,'CAA':0.0,'CAC':0.0,'CAG':0.0,'CAT':0.0,'CCA':0.0,'CCC':0.0,'CCG':0.0, \
                                    'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0,'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0, \
                                    'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0,'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0, \
                                    'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0,'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0, \
                                    'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0}}   


    # Mutate (substitutions) all possible codons in the given genetic code, and count proportions of mutations that are synonymous and non-synonmyous
    for codon in nt_to_aa.keys():

        # for each bp position in the codon...
        for codon_p in range(0,2+1):

            nts = ['A','G','T','C']  
            nts.remove(codon[codon_p]) # we do not consider self substitutions, e.g. A->A

            # ...and for each nucleotide that the bp can change 
            
            for nt in nts:

                codon_mutated = list(copy.deepcopy(codon))
                codon_mutated[codon_p] = nt  # mutate the basepair
                codon_mutated = ''.join(codon_mutated)
                
                # ...count how many of them are synonymous.
                if nt_to_aa[codon]==nt_to_aa[codon_mutated]:
                    potential_changes['S'][codon]+=1.0/3.0 
                else:
                    potential_changes['N'][codon]+=1.0/3.0
    return potential_changes

def get_feature_codons(genome,annotation):

    dct = {key:str(genome.seq)[(item[0] -1):item[1]] for key,item in annotation.items()}

    return {key:split_into_codons(item) for key,item in dct.items()}

def get_df(codons):

    keys = []
    N_sites = []
    S_sites = []

    values =  potential_changes_dict(geneticCode(snakemake.params.genetic_code))

    for key,item in codons.items():
        keys.append(key)
        N = sum([values['N'][x] for x in item if x in values['N'].keys()])
        S = sum([values['S'][x] for x in item if x in values['S'].keys()])
        N_sites.append(N)
        S_sites.append(S) 
    
    df =pd.DataFrame({'gene':keys, 'N':N_sites, 'S':S_sites})

    return df

def main():
    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)
    
    sequence = list(SeqIO.parse(snakemake.input.fasta,"fasta"))[0]
    codons_dict =  get_feature_codons(sequence,SCov2_annotation)
    df = get_df(codons_dict)

    df.to_csv(snakemake.output.csv,index= False)



if __name__ == "__main__":
    main()