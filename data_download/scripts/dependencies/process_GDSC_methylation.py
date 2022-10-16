import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from functools import reduce
from intervaltree import Interval, IntervalTree
from GDSC_methods import filter_protein_coding, filter_mini_cancer, is_probe_in_tree, create_tree

print(os.getcwd())
raw_folder = '../../data/'
output_folder = '../../data/'

# Read models
model_file = '%sraw/methSampleId_2_cosmicIds.xlsx'%(raw_folder)
model_df = pd.read_excel(model_file)
model_df = model_df[model_df['Cell_Line'] == 'YES']
model_df['METH_ID'] = ['%s_%s'%(e,f)
                       for e,f in zip(model_df['Sentrix_ID'], model_df['Sentrix_Position'])]
model_df = model_df[['METH_ID', 'Tissue', 'Sample_Name']]

# Read methylation
methylation_file = '%sraw/F2_METH_CELL_DATA.txt'%(raw_folder)
meth_df = pd.read_csv(methylation_file, sep='\t', index_col=0)
meth_df = meth_df.T.merge(model_df.set_index('METH_ID'),
                          left_index=True,
                          right_index=True,
                          how='left')
meth_df = meth_df.groupby(['Sample_Name', 'Tissue']).agg('mean')
meth_df.index = meth_df.index.rename(['model_name', 'tissue'])
meth_df.index = meth_df.index.set_levels([meth_df.index.levels[0].astype(str), 
                                        meth_df.index.levels[-1]])
meth_df.index = pd.MultiIndex.from_tuples([(str(x[0]), str(x[1]).upper()) for x in meth_df.index])
meth_df.to_pickle('%s/all/methylation/all_methylation_data.pkl'%(output_folder),
                  compression='gzip')

# Protein coding
genes_lookup_file = '%s/pybiomart_gene_status.csv'%(raw_folder)
genes_lookup_df = pd.read_csv(genes_lookup_file, sep='\t', index_col=0)
pc_interval_trees = create_tree(meth_df, genes_lookup_df[genes_lookup_df.status == 'protein_coding'])
protein_coding_probes = []
for pos in meth_df.columns:
    if len(is_probe_in_tree(pos, pc_interval_trees)) >= 1:
        protein_coding_probes.append(pos)
meth_pc_df = meth_df[protein_coding_probes]
meth_pc_df.to_pickle('%s/protein_coding/methylation/all_methylation_data.pkl'%(output_folder), 
                     compression='gzip')

# Mini cancer
mini_cancer_file = './mini_cancer_lookup_genes.csv'
mini_cancer_df = pd.read_csv(mini_cancer_file, index_col=None)
mini_cancer_genes = mini_cancer_df['Hugo'].values
mini_interval_trees = create_tree(meth_df, genes_lookup_df[genes_lookup_df['Hugo'].isin(mini_cancer_genes)])
mini_probes = []
for pos in meth_df.columns:
    if len(is_probe_in_tree(pos, mini_interval_trees)) >= 1:
        mini_probes.append(pos)        
meth_mini_df = meth_df[mini_probes]
meth_mini_df.to_pickle('%s/mini_cancer/methylation/all_methylation_data.pkl'%(output_folder), 
                       compression='gzip')