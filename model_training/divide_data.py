"""
DIVIDE DATA PRIOR TO AIC AND COMPUTE DISPERSION

2021/12/22
"""

import os, sys, torch, getopt
import numpy as np
import pandas as pd
from pickle import dump
from sklearn.model_selection import KFold
from tqdm import tqdm
from percolate import GLMPCA
from percolate.negative_binomial_routines import compute_dispersion

sys.path.insert(0, './read_data/')
from read_GDSC_data import read_GDSC_data
from read_GDSC_response import read_GDSC_response
import library_size_normalization

from params import naming, data_folder, n_cv, family, n_jobs, sample_subdiv, n_subdiv

opts, args = getopt.getopt(sys.argv[1:],'o:m:t:d:',['output=', 'maxeval=', 'tmpfolder=', 'datatype='])
for opt, arg in opts:
    if opt in ("-o", "--output"):
        output_folder = str(arg)
    if opt in ("-d", "--datatype"):
        data_types = [str(arg)]

###
# IMPORT DATA
###


data_df = read_GDSC_data(
    data_types=[naming[e] for e in data_types],
    tissues=['all'],
    projects=None,
    filtering=None,
    folder=data_folder
)

data_key = list(data_df.keys())[0]
data_df = data_df[data_key]

#Retrieve names of features and samples
feature_names = data_df.columns
sample_names = data_df.index

# Library size of gene expression
if 'count' in data_key:
    data_df = library_size_normalization.TMM_normalization(
        data_df.values.astype(float)
    )
    data_df = pd.DataFrame(np.array(data_df))

# Transform data into Tensors
data_df = torch.Tensor(data_df.values)

###
# SAVE DATA
###

if 'data' not in os.listdir(output_folder):
    os.mkdir('%s/data'%(output_folder))

if 'cv_folds' not in os.listdir(output_folder):
    os.mkdir('%s/cv_folds'%(output_folder))
cv_folder = '%s/cv_folds'%(output_folder)

# Set up division
split_folds = KFold(n_splits=n_cv, shuffle=True)
for cv_idx, (train, test) in enumerate(split_folds.split(data_df)):
    if 'fold_%s_train_idxs.npy'%(cv_idx) in os.listdir(cv_folder):
        continue
    np.save('%s/fold_%s_train_idxs.npy'%(cv_folder, cv_idx), np.array(train))
    np.save(
        '%s/fold_%s_train_samples.npy'%(cv_folder, cv_idx),
        np.array(sample_names.get_level_values(0)[train]).astype(str)
    )
    np.save('%s/fold_%s_test_idxs.npy'%(cv_folder, cv_idx), np.array(test))
    np.save(
        '%s/fold_%s_test_samples.npy'%(cv_folder, cv_idx),
        np.array(sample_names.get_level_values(0)[test]).astype(str)
    )
    
# Save data
torch.save(data_df, '%s/data.pt'%(output_folder))

###
# COMPUTE DISPERSION
###

exp_family_params = {}
if 'nb' in family[data_types[0]]:
    if sample_subdiv < 0:
        pass
    else:
        data_df = data_df[np.random.choice(np.arange(data_df.shape[0]), sample_subdiv)]

    if n_subdiv > 1:
        sub_idx = np.linspace(0, data_df.shape[1],n_subdiv).astype(int)
        dispersion_parameters = [
            compute_dispersion(pd.DataFrame(data_df[:,sub_idx[idx]:sub_idx[idx+1]])).values
            for idx in tqdm(range(n_subdiv-1))
        ]
        dispersion_parameters = np.concatenate(dispersion_parameters)
        r_coef = 1. / torch.Tensor(dispersion_parameters).flatten()
    else:
        r_coef = 1. / torch.Tensor(compute_dispersion(pd.DataFrame(data_df.detach().numpy())).values).flatten()
    
    gene_filter = torch.where((r_coef > 0.01) & (~torch.isnan(r_coef)) & (r_coef < 10**4))[0]

    exp_family_params['r'] = torch.Tensor(r_coef)
    exp_family_params['gene_filter'] = torch.where((r_coef > 0.01) & (~torch.isnan(r_coef)) & (r_coef < 10**4))[0]
else:
    glmpca_clf_sat = GLMPCA(n_pc=-1,family=family[data_types[0]], n_jobs=n_jobs)
    sat_params = glmpca_clf_sat.compute_saturated_params(
        data_df, with_intercept=False, exp_family_params=None, save_family_params=True
    )
    exp_family_params = glmpca_clf_sat.exp_family_params
    
# Save data
dump(exp_family_params, open('%s/dispersion.pkl'%(output_folder), 'wb'))
