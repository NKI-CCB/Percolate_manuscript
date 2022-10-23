"""
COMPUTE TRAIN LIKELIHOOD FOR VARIOUS PARAMETERS

Computes the AIC of the data for various parameters.

2021/12/22
"""


import sys, os, getopt, uuid, torch, itertools, datetime, pickle, logging
import numpy as np
import pandas as pd
from sklearn.model_selection import ParameterGrid
from hyperopt import STATUS_FAIL, STATUS_OK, Trials, fmin, hp, tpe
from percolate import GLMPCA
from percolate.exponential_family import *
from params import *

opts, args = getopt.getopt(sys.argv[1:],'o:m:t:d:',['output=', 'maxeval=', 'tmpfolder=', 'datatype='])
for opt, arg in opts:
    if opt in ("-o", "--output"):
        output_folder = str(arg)
    if opt in ("-t", "--tmpfile"):
        tmp_folder = str(arg)
    if opt in ("-d", "--datatype"):
        datatype = str(arg)

if 'gridsearch' not in os.listdir(output_folder):
    os.mkdir('%s/gridsearch/'%(output_folder))

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

family = family[datatype]
grid_search_space = grid_search_space[datatype]

###
# IMPORT AND PROCESS DATA
###

X_data = torch.load(open('%s/data.pt'%(output_folder), 'rb'))
exp_family_params = pickle.load(open('%s/dispersion.pkl'%(output_folder), 'rb'))

###
# COMPUTE SATURATED PARAMETERS
###

glmpca_clf = GLMPCA(n_pc=0, family=family, n_jobs=40)
saturated_parameters = glmpca_clf.compute_saturated_params(
    X_data, with_intercept=False, exp_family_params=exp_family_params, save_family_params=True 
)

###
# GRID SEARCH
###


def _objective_function(params):
    """
    Compute likelihood.
    """
    n_pc = params['n_pc']
    n_glmpca_init = params['n_glmpca_init']
    batch_size = params['batch_size']
    learning_rate = params['learning_rate']
    max_param = params['max_param']
    maxiter = params['maxiter']
    step_size = params['step_size']
    gamma = params['gamma']
    
    glmpca_clf = GLMPCA(
        n_pc=n_pc,
        family=family,
        maxiter=maxiter,
        max_param=max_param,
        learning_rate=learning_rate,
        step_size=step_size,
        gamma=gamma
    )

    try:
        glmpca_clf.compute_saturated_loadings(
            X_data, 
            batch_size=batch_size, 
            n_init=n_glmpca_init, 
            exp_family_params=exp_family_params,
            saturated_params=saturated_parameters
        )
    except ValueError:
        results_dict = {'status': STATUS_FAIL, 'loss': np.iinfo(np.uint64).max}
        print('\n ERROR: \n %s \n'%(params))
        results_dict.update(params)
        pickle.dump(results_dict, open('%s/iter_%s.pkl'%(tmp_folder, str(uuid.uuid4()).split('-')[1]), 'wb'))
        return results_dict

    # Compute performance scores
    params_copy = {
        k: glmpca_clf.exp_family_params[k].to(device) if type(glmpca_clf.exp_family_params[k]) is torch.Tensor else glmpca_clf.exp_family_params[k] 
        for k in glmpca_clf.exp_family_params
    }
    if 'r' in params_copy:
        params_copy['r'] = params_copy['r'][params_copy['gene_filter']].to(device)

    likelihood_computation = make_saturated_loading_cost(
        family, max_value=max_param, params=params_copy, train=False
    )
    del params_copy
    
    data_likelihood = likelihood_computation(
        X=glmpca_clf.saturated_loadings_.to(device),
        data=(X_data[:,glmpca_clf.exp_family_params['gene_filter']] if 'gene_filter' in glmpca_clf.exp_family_params else X_data).to(device),
        parameters=saturated_parameters,
        intercept=glmpca_clf.saturated_intercept_.to(device),
    ).cpu().detach().numpy()

    # AIC correction
    results_dict = {
        'loss': data_likelihood,
        'loss_choice': 'CV_likelihood',
        'AIC_supp': n_pc * X_data.shape[1] - n_pc*(n_pc+1)/2,
        'status': STATUS_OK,
    }
    
    # Remove cache
    del glmpca_clf
    torch.cuda.empty_cache()

    results_dict.update(params)
    pickle.dump(results_dict, open('%s/iter_%s_GS_AIC.pkl'%(tmp_folder, str(uuid.uuid4()).split('-')[1]), 'wb'))
    return results_dict

grid_params = ParameterGrid(grid_search_space)
gridsearch_results_df = []
for idx, params in enumerate(grid_params):
    print('START ITER %s %s/%s'%(datatype, idx, len(grid_params)), flush=True)
    gridsearch_results_df.append(_objective_function(params))
gridsearch_results_df = pd.DataFrame(gridsearch_results_df)

# Save
save_id = '{:%B_%d_%Y_}'.format(datetime.datetime.now()) + '_' + str(uuid.uuid4()).split('-')[1]
gridsearch_results_df.to_csv('%s/gridsearch/results_AIC_%s.csv'%(output_folder, save_id))




