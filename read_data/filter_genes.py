# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

2019/12/18

FILTER GENES
"""

import numpy as np
import pandas as pd

pybiomart_queried_file = '/DATA/s.mourragui/data/2019_04_cell_line_data/pybiomart_gene_status.csv'
mini_cancer_file = '/DATA/s.mourragui/data/2019_04_cancer_mini_genome/hugo_genes.tsv'


def filter_protein_coding(df, gene_caract_file=pybiomart_queried_file):
    """
    Given a DataFrame df with genes in the columns, filter the genes to
    return only protein coding genes.

    -------
    df: pandas.DataFrame
        Genomic data with samples in the rows and genes in the columns.

    gene_caract_file: str, default to DATA location
        Location of the biomart file, i.e. lookup table with a column 'status'
        and a column 'Hugo' for each gene.

    Returned Values
    -------
    pandas.DataFrame with samples in the rows and genes in the columns.
    """
    gene_caract = pd.read_csv(gene_caract_file, sep='\t', index_col=0)
    gene_caract = gene_caract[gene_caract.status == 'protein_coding']
    pc_genes = gene_caract['Hugo'].values
    pc_genes = np.intersect1d(pc_genes, df.columns).astype(str)

    return df[pc_genes]


def filter_mini_cancer(df, gene_caract_file=mini_cancer_file):
    """
    Given a DataFrame df with genes in the columns, filter the genes to
    return only genes part of the mini-cancer genome.
    Can be applied to any other list of genes.

    -------
    df: pandas.DataFrame
        Genomic data with samples in the rows and genes in the columns.

    gene_caract_file: str, default to DATA location
        Location of the list of genes, i.e. file with one line per gene.

    Returned Values
    -------
    pandas.DataFrame with samples in the rows and genes in the columns.
    """
    mini_cancer_file = '/DATA/s.mourragui/data/2019_04_cancer_mini_genome/hugo_genes.tsv'
    mini_genes = pd.read_csv(mini_cancer_file).values.flatten().astype(str)

    mini_genes = np.intersect1d(mini_genes, df.columns)
    return df[mini_genes]

