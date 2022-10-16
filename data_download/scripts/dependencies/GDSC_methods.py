import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree
import re

# Function to assign a probe to a class of genes
def is_probe_in_tree(pos,tree):
    chromosome = re.findall(r'chr([0-9]+)', pos)[0]
    start = int(re.findall(r':([0-9]+)-', pos)[0])
    end = int(re.findall(r'-([0-9]+)$', pos)[0])
    
    return tree[chromosome].overlap(start, end)


def create_tree(df, lookup_df):
    # We create one IntervalTree per chromosome and fit them into dictionary
    unique_chromosomes = np.unique([e.split(':')[0] for e in df.columns])
    unique_chromosomes = [re.findall(r'([0-9]+)', e)[0] for e in unique_chromosomes]
    interval_trees = {c:IntervalTree() for c in unique_chromosomes}

    # Fill the interval trees
    for gene in lookup_df.iterrows():
        chromosome = gene[1]['chromosome_name']
        if chromosome not in unique_chromosomes:
            continue
        gene_start = int(gene[1]['gene_start'])
        gene_end = int(gene[1]['gene_end'])
        interval_trees[str(chromosome)][gene_start:gene_end] = gene[1]['Hugo']

    return interval_trees


# Function to filter protein coding genes
def filter_protein_coding(df, pc_file):
    gene_status_df = pd.read_csv(pc_file, sep='\t', index_col=0)
    protein_coding_genes = gene_status_df[gene_status_df['status'] == 'protein_coding']
    protein_coding_genes = protein_coding_genes['Hugo'].values
    protein_coding_genes = np.intersect1d(protein_coding_genes, df.columns).astype(str)

    return df[protein_coding_genes]


# Function to filter mini cancer genes
def filter_mini_cancer(df, mini_file):
  mini_cancer_df = pd.read_csv(mini_file, index_col=None)
  mini_cancer_df = mini_cancer_df['Hugo'].values
  mini_cancer_genes = np.intersect1d(df.columns, mini_cancer_df)

  return df[mini_cancer_genes]


# Function to read rnaseq data (count and fpkm)
def read_data_rnaseq(file, model_df):
    df = pd.read_csv(file, index_col=[0,1], header=[0,1,2])
    df.columns = df.columns.droplevel(1).droplevel(1)
    df.index = df.index.droplevel(0)
    df = df.T
    df = df.dropna(axis=1)
    df = df.merge(model_df[['model_id', 'model_name', 'tissue']].set_index('model_id'),
                  how='left',
                  left_index=True,
                  right_index=True)
    df = df.reset_index(drop=True).set_index(['model_name', 'tissue'])
    df = df.dropna(axis=1)
    df.index = pd.MultiIndex.from_tuples([(x[0], str(x[1]).upper()) for x in df.index])
    
    return df


# Function to read mutation data
def read_data_mutations(file, model_df, filter_functional=False):
    df = pd.read_csv(file, index_col=None, header=[0,1,2])
    df.columns = df.columns.droplevel(1).droplevel(1)
    if filter_functional:
        df = df[df['protein_mutation'] != '-']
    df = df[['model_id', 'gene_symbol']]
    df['INDIC'] = 1
    df = df.pivot_table(columns=['gene_symbol'], index='model_id', fill_value=0)
    df.columns = df.columns.droplevel(0)
    df = df.merge(model_df[['model_name', 'tissue', 'model_id']].set_index('model_id'),
                  how='left',
                  left_index=True,
                  right_index=True)
    df = df.reset_index(drop=True).set_index(['model_name', 'tissue'])
    df = df.dropna(axis=1)
    df.index = pd.MultiIndex.from_tuples([(x[0], str(x[1]).upper()) for x in df.index])

    return df

def read_data_copy_number(file, model_df):
    cnv_df = pd.read_csv(file, header=[0,1,2])
    cnv_df.columns = cnv_df.columns.droplevel(1).droplevel(1)
    del cnv_df['model_id']
    cnv_df.columns = ['gene_name'] + list(cnv_df.columns)[1:]
    cnv_df.set_index('gene_name', inplace=True)
    cnv_df = cnv_df.T
    cnv_df = cnv_df.merge(model_df[['model_id', 'tissue', 'model_name']].set_index('model_id'),
                          how='left',
                          left_index=True,
                          right_index=True)
    cnv_df = cnv_df.reset_index(drop=True)
    cnv_df = cnv_df.set_index(['model_name', 'tissue'])
    cnv_df = cnv_df.dropna(axis=1)
    cnv_df.index = pd.MultiIndex.from_tuples([(x[0], str(x[1]).upper()) for x in cnv_df.index])

    return cnv_df

def read_crispr(file, model_df):
    df = pd.read_csv(file, sep='\t', index_col=0).T
    df = df.merge(model_df[['model_name', 'model_id', 'tissue']].set_index('model_id'),
                how='left',
                left_index=True,
                right_index=True)
    df.index.names = ['model_id']
    df.reset_index(drop=True).set_index(['model_name', 'tissue'], inplace=True)
    df.index = pd.MultiIndex.from_tuples([(x[0], str(x[1]).upper()) for x in df.index])

    return df

