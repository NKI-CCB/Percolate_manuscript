# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

2019/12/18

READ GDSC DATA

NOTES:
- 2020/01/29: Add read_one_GDSC_modality for incorporating
and concatenating more tissues.
- 2020/10/13: Modification to new da
"""

import os
import pandas as pd
import numpy as np
from functools import reduce


def read_GDSC_data(data_types,
                   tissues,
                   projects=None,
                   filtering=None,
                   folder='/DATA/s.mourragui/data/2020_10_GDSC_data'):
    
    # If folder is None, set to default
    folder = folder or '/DATA/s.mourragui/data/2020_10_GDSC_data'

    # Change data_folder
    if filtering in ['mini', 'mini_cancer']:
        data_folder = folder + '_mini_cancer/'
    elif filtering in ['protein_coding', 'pc']:
        data_folder = folder + '_protein_coding/'
    elif filtering in ['cancer_genes', 'cancer']:
        data_folder = folder + '_cancer_genes/'
    else:
        if filtering is not None:
            print('WARNING: %s not recognized as filtering'%(filtering))
        data_folder = folder + '/'
    
    if tissues is None:
        tissues = 'all'
    if type(tissues) == str:
        tissues = [tissues]
    if type(data_types) == str:
        data_types = [data_types]
        
    data_df = {'%s_%s'%('-'.join(tissues),d): read_one_GDSC_modality(tissues, projects, d, data_folder)
               for d in data_types}
    
    # Reduce to same samples
    common_idxs = reduce(lambda x,y: x.intersection(y),
                         [data_df[d].index.get_level_values(0) for d in data_df])
    common_idxs = np.unique(common_idxs)

    data_df = {d:data_df[d].loc[common_idxs] for d in data_df}
    for e in data_df:
        for d in data_df:
            pd.testing.assert_index_equal(data_df[e].index.get_level_values(0),
                                          data_df[d].index.get_level_values(0))
    
    return data_df


def read_one_GDSC_modality(tissues, projects, data_type, folder):
    if type(tissues) != list:
        tissues = [tissues]
    if type(projects) != list:
        projects = [projects]
    
    folders = ['%s/%s/%s%s_%s_data.pkl'%(folder,
                                        data_type,
                                        t.replace(' ', '_').lower(),
                                        '_%s'%(p) if p else '',
                                        data_type)
              for t in tissues for p in projects]
    folders = np.array(folders)
    
    folders = folders[[os.path.exists(f) for f in folders]]
    if folders.shape[0] == 0:
        return pd.DataFrame()

    data_df = [pd.read_pickle(f, compression='gzip') for f in folders]

    # Take same features
    unique_features = reduce(np.intersect1d, [df.columns for df in data_df])
    data_df = pd.concat([df[unique_features] for df in data_df])
    return data_df


def _useful_pivot_columns(return_tissue):
    if return_tissue:
        useful_columns = ['gene_symbol', 'model_name', 'tissue']
        renamed_columns = ['gene', 'sample', 'tissue']
    else:
        useful_columns = ['gene_symbol', 'model_name']
        renamed_columns = ['gene', 'sample']

    return useful_columns, renamed_columns