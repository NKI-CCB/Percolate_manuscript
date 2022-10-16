import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from GDSC_methods import read_data_copy_number, filter_protein_coding, filter_mini_cancer

print(os.getcwd())
raw_folder = '../../data/'
output_folder = '../../data/'

# Read models
model_file = '%smodel_list_20191104.csv'%(raw_folder)
model_df = pd.read_csv(model_file)

# Load data
copy_number_file_gistic = '%sraw/cnv_gistic_20191101.csv'%(raw_folder)
copy_number_file_picnic = '%sraw/cnv_abs_copy_number_picnic_20191101.csv'%(raw_folder)
cnv_df_gistic = read_data_copy_number(file=copy_number_file_gistic, model_df=model_df)
cnv_df_picnic = read_data_copy_number(file=copy_number_file_picnic, model_df=model_df)

# Save all genes
cnv_df_gistic.to_pickle('%s/all/copy_number_binary/all_copy_number_binary_data.pkl'%(output_folder), compression='gzip')
cnv_df_picnic.to_pickle('%s/all/copy_number/all_copy_number_data.pkl'%(output_folder), compression='gzip')

# Protein coding
protein_coding_file = '%s/pybiomart_gene_status.csv'%(raw_folder)
cnv_pc_df_gistic = filter_protein_coding(cnv_df_gistic, pc_file=protein_coding_file)
cnv_pc_df_picnic = filter_protein_coding(cnv_df_picnic, pc_file=protein_coding_file)
cnv_pc_df_gistic.to_pickle('%s/protein_coding/copy_number_binary/all_copy_number_binary_data.pkl'%(output_folder),
                           compression='gzip')
cnv_pc_df_picnic.to_pickle('%s/protein_coding/copy_number/all_copy_number_data.pkl'%(output_folder),
                           compression='gzip')

# Mini cancer
# Mini cancer
mini_cancer_file = './mini_cancer_lookup_genes.csv'
cnv_mini_df_gistic = filter_mini_cancer(cnv_df_gistic, mini_file=mini_cancer_file)
cnv_mini_df_picnic = filter_mini_cancer(cnv_df_picnic, mini_file=mini_cancer_file)
cnv_mini_df_gistic.to_pickle('%s/mini_cancer/copy_number_binary/all_copy_number_binary_data.pkl'%(output_folder),
                             compression='gzip')
cnv_mini_df_picnic.to_pickle('%s/mini_cancer/copy_number/all_copy_number_data.pkl'%(output_folder),
                             compression='gzip')