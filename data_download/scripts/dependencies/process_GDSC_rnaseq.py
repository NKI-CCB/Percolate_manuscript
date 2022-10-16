import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from GDSC_methods import read_data_rnaseq, filter_protein_coding, filter_mini_cancer

raw_folder = '../../data/'
output_folder = '../../data/'

# Read models
model_file = '%smodel_list_20191104.csv'%(raw_folder)
model_df = pd.read_csv(model_file)

# Read data
rnaseq_df = read_data_rnaseq(file='%sraw/rnaseq_read_count_20191101.csv'%(raw_folder), 
                            model_df=model_df)
fpkm_df = read_data_rnaseq(file='%sraw/rnaseq_fpkm_20191101.csv'%(raw_folder), 
                        model_df=model_df)

# Save all
rnaseq_df.to_pickle('%s/all/count/all_count_data.pkl'%(output_folder), compression='gzip')
fpkm_df.to_pickle('%s/all/fpkm/all_count_data.pkl'%(output_folder), compression='gzip')

# Restrict to protein coding genes
gene_status_df = pd.read_csv('%s/pybiomart_gene_status.csv'%(raw_folder), sep='\t', index_col=0)
protein_coding_genes = gene_status_df[gene_status_df['status'] == 'protein_coding']
protein_coding_genes = protein_coding_genes['Hugo'].values
protein_coding_genes = np.intersect1d(protein_coding_genes, rnaseq_df.columns).astype(str)

# Save protein coding
protein_coding_file = '%s/pybiomart_gene_status.csv'%(raw_folder)
rnaseq_pc_df = filter_protein_coding(rnaseq_df, protein_coding_file)
rnaseq_pc_df.to_pickle('%s/protein_coding/count/all_count_data.pkl'%(output_folder), compression='gzip')
fpkm_pc_df = filter_protein_coding(fpkm_df, protein_coding_file)
fpkm_pc_df.to_pickle('%s/protein_coding/fpkm/all_count_data.pkl'%(output_folder), compression='gzip')

# Save mini cancer
mini_cancer_file = './mini_cancer_lookup_genes.csv'
rnaseq_mini_df = filter_mini_cancer(rnaseq_df, mini_cancer_file)
rnaseq_mini_df.to_pickle('%s/mini_cancer/count/all_count_data.pkl'%(output_folder),
                        compression='gzip')
fpkm_mini_df = filter_mini_cancer(fpkm_df, mini_cancer_file)
fpkm_mini_df.to_pickle('%s/mini_cancer/fpkm/all_count_data.pkl'%(output_folder),
                        compression='gzip')