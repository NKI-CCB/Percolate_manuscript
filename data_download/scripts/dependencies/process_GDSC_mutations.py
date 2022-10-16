import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from GDSC_methods import read_data_mutations, filter_protein_coding, filter_mini_cancer

print(os.getcwd())
raw_folder = '../../data/'
output_folder = '../../data/'

# Read models
model_file = '%smodel_list_20191104.csv'%(raw_folder)
model_df = pd.read_csv(model_file)

# Read mutations
mut_file = '%sraw/mutations_20191101.csv'%(raw_folder)
mut_df = read_data_mutations(file=mut_file, model_df=model_df, filter_functional=False)
ns_mut_df = read_data_mutations(file=mut_file, model_df=model_df, filter_functional=True)

# Save all
mut_df.to_pickle('%s/all/mutations/all_mutations_data.pkl'%(output_folder), compression='gzip')
ns_mut_df.to_pickle('%s/all/ns_mutations/all_ns_mutations_data.pkl'%(output_folder), compression='gzip')

# Save protein coding
protein_coding_file = '%s/pybiomart_gene_status.csv'%(raw_folder)
mut_pc_df = filter_protein_coding(df=mut_df, pc_file=protein_coding_file)
mut_pc_df.to_pickle('%s/protein_coding/mutations/all_mutations_data.pkl'%(output_folder), compression='gzip')
ns_mut_pc_df = filter_protein_coding(df=ns_mut_df, pc_file=protein_coding_file)
ns_mut_pc_df.to_pickle('%s/protein_coding/ns_mutations/all_ns_mutations_data.pkl'%(output_folder), compression='gzip')

# Save mini cancer
mini_cancer_file = './mini_cancer_lookup_genes.csv'
mut_mini_df = filter_mini_cancer(df=mut_df, mini_file=mini_cancer_file)
mut_mini_df.to_pickle('%s/mini_cancer/mutations/all_mutations_data.pkl'%(output_folder), compression='gzip')
ns_mut_mini_df = filter_mini_cancer(df=ns_mut_df, mini_file=mini_cancer_file)
ns_mut_mini_df.to_pickle('%s/mini_cancer/ns_mutations/all_ns_mutations_data.pkl'%(output_folder), compression='gzip')