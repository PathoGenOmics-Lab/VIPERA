#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Script adaptado a nuestro input de https://github.com/PathoGenOmics-Lab/genetic-distances
import pandas as pd
from Bio import SeqIO
import os
import logging
from typing import List
import numpy as np
import multiprocessing
from joblib import Parallel, delayed, effective_n_jobs


def parse_vcf():
    """
    Parse a VCF containing positions for masking. Assumes the VCF file is
    formatted as:
    github.com/W-L/ProblematicSites_SARS-CoV2/blob/master/problematic_sites_sarsCov2.vcf
    with a "mask" or "caution" recommendation in column 7.
    Masked sites are specified with params.
    """
    vcf = pd.read_csv(
        snakemake.config["PROBLEMATIC_VCF"],
        delim_whitespace=True,
        comment="#",
        names=("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
    )
    positions = tuple(vcf.loc[vcf.FILTER.isin(snakemake.params.mask_class), "POS"])
    return positions

def create_freq_dict(input_file) -> dict:
    alt_values = set()  # Utilizamos un conjunto para evitar duplicados

    # Sacamos todos los posibles alelos alternativos
    
    df = pd.read_table(input_file, sep = "\t")
    alt_values.update(df['ALT'].unique())

    # Crea el diccionario freq a partir del conjunto de valores de 'ALT', con un valor inicial de 0
    freq = {alt_value: 0 for alt_value in alt_values}

    return freq


def tsv_from_seq(tsv_reference,reference, reference_name):

    mask = parse_vcf()
    seq = tsv_reference
    reference = reference

    pos = []
    ALT = []
    
    for i in range(0,len(seq)):
        
        if i not in mask:
              
            pos.append(i+1)
            ALT.append(reference[i])
                      
    df = pd.DataFrame({"POS":pos,"ALT":ALT})
    df["ALT_FREQ"] = 1
    df["REGION"] = reference_name
    
    return df

def get_pos_tup(df, pos, reference, freq):
    freq = freq.copy()  # Hace una copia de freq para no modificar el diccionario original
    alt_keys = sorted(list(freq.keys()))  # Crea una lista ordenada de las llaves

    # Si la posicion esta en el dataframe, obtiene las frecuencias alelicas
    if pos + 1 in df["POS"].values:
        df_ = df[df["POS"] == pos+1]
        for base in alt_keys:
            if base in df_["ALT"].values:
                freq[base] = float(df_["ALT_FREQ"][df_["ALT"] == base].iloc[0])

    # Calcula la frecuencia de la base en la secuencia de referencia
    ref = 1 - sum(freq.values())
    freq[reference[pos]] += ref
    # Retorna una tupla de las frecuencias alelicas en el orden de alt_keys
    return tuple(freq[base] for base in alt_keys)

def calc_heterozygosities(df1, df2, pos, reference, freq):
    freqs1 = get_pos_tup(df1, pos, reference, freq)
    freqs2 = get_pos_tup(df2, pos, reference, freq)
    
    hs1 = heterozygosity(freqs1)
    hs2 = heterozygosity(freqs2)
    hs = (hs1 + hs2) / 2

    total_freqs = [(f1 + f2) / 2 for f1, f2 in zip(freqs1, freqs2)]
    ht = heterozygosity(total_freqs)
    return hs, ht

def heterozygosity(freqs):
    return 1 - sum([f ** 2 for f in freqs])

def calc_fst_weir_cockerham(hs, ht):
    return (ht - hs) / ht if ht != 0 else 0

def get_dif_n(df, COV1, COV2, mask, reference, freq):

    df1 = df[df["REGION"] == COV1]
    df2 = df[df["REGION"] == COV2]

    return sum([calc_fst_weir_cockerham(*calc_heterozygosities(df1, df2, i, reference, freq)) 
                for i in range(len(reference)) if i + 1 not in mask])

def get_matrix(df, cov_list, reference, freq):
    distance_matrix = {}
    mask_positions = parse_vcf()

    def calculate_distance(sample):
        print(f"Calculating distances for sample: {sample}")
        distances = [get_dif_n(df, sample, x, mask_positions, reference, freq) for x in cov_list]
        return distances
    num_jobs = snakemake.threads
    #num_jobs = multiprocessing.cpu_count()  # Obten el numero de nucleos de la CPU
    if len(cov_list) <= effective_n_jobs(num_jobs):
        batch_size = 1  # Establecer batch_size en 1 si el tamanyo de cov_list es menor o igual al numero de trabajos paralelos
    else:
        batch_size = int(len(cov_list) / effective_n_jobs(num_jobs))  # Tamanyo de lote optimo

    print(f"Parallelizing the calculation with {num_jobs} jobs...")
    results = Parallel(n_jobs=num_jobs, batch_size=batch_size, verbose=10, timeout=None)(
        delayed(calculate_distance)(sample) for sample in cov_list
    )

    for i, sample in enumerate(cov_list):
        distance_matrix[sample] = results[i]

    # Llenar la parte superior de la matriz reflejando los valores
    for i in range(len(cov_list)):
        for j in range(i+1, len(cov_list)):
            distance_matrix[cov_list[j]][i] = distance_matrix[cov_list[i]][j]

    return pd.DataFrame(distance_matrix, index=cov_list)

def read_and_concatenate_tsvs(input,tsv_reference,reference, reference_name):
    df_list = []
    positions = parse_vcf()  # Get the mask positions

    
    df_1 = pd.read_table(input, sep = "\t")
    df_list.append(df_1)
    df_list.append(tsv_from_seq(tsv_reference,reference,reference_name))
    concatenated_df = pd.concat(df_list, ignore_index=True)
    return concatenated_df


def main():    
    logging.basicConfig(filename=snakemake.log[0], format=snakemake.config["LOG_PY_FMT"], level=logging.INFO)
    input_file = snakemake.input.tsv
    reference = str(next(SeqIO.parse(snakemake.params.tsv_reference, "fasta")).seq)
    outgroup = str(next(SeqIO.parse(snakemake.params.reference, "fasta")).seq)
    outgroup_name = str(next(SeqIO.parse(snakemake.params.reference, "fasta")).id)
    df = read_and_concatenate_tsvs(input_file,reference,outgroup,outgroup_name)
    cov_list = df["REGION"].unique().tolist()
    freq = create_freq_dict(input_file)
    df = get_matrix(df, cov_list, reference, freq)
    df.to_csv(snakemake.output.distances)

if __name__ == "__main__":
    main()