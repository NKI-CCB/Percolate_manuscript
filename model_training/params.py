from hyperopt import hp
import numpy as np

data_folder = './data/mini_cancer/'
naming = {
    'MUT': 'mutations',
    'GE': 'count',
    'CNA': 'copy_number',
    'METH': 'methylation'
}
n_cv = 3
n_jobs = 10
sample_subdiv = 200 # Number of samples to take for computing the dispersion (-1 takes all)
n_subdiv = 50

family = {
    'MUT': 'bernoulli',
    'GE': 'nb_rep',
    'METH': 'beta_rep',
    'CNA': 'gamma'#'gamma'
}

grid_search_space = {}
if sample_subdiv < 250:
    grid_search_space['MUT'] = {
        'n_pc': np.linspace(2,100,5).astype(int),
        'n_glmpca_init': [2],
        'learning_rate': [0.1,0.5,1.],
        'batch_size': [8],
        'max_param': [25],
        'maxiter': [10],
        'gamma': [0.5],
        'step_size': [20]
    }

    grid_search_space['GE'] = {
        'n_pc': np.linspace(2,100,5).astype(int),
        'n_glmpca_init': [2],
        'learning_rate': [0.5, 1., 10.],
        'batch_size': [16],
        'max_param': [25],
        'maxiter': [10],
        'gamma': [0.5],
        'step_size': [20]
    }

    grid_search_space['METH'] = {
        'n_pc': np.linspace(2,150,5).astype(int),
        'n_glmpca_init': [5],
        'learning_rate': [0.05,0.1,0.5],
        'batch_size': [16],
        'max_param': [25],
        'maxiter': [10],
        'gamma': [0.5],
        'step_size': [10]
    }
    grid_search_space['CNA'] = {
        'n_pc': np.linspace(2,100,5).astype(int),
        'n_glmpca_init': [2],
        'learning_rate': [0.5, 1., 10.],
        'batch_size': [16],
        'max_param': [np.inf],
        'maxiter': [10],
        'gamma': [0.5],
        'step_size': [20]
    }
else:
    grid_search_space['MUT'] = {
        'n_pc': np.linspace(2,100,50).astype(int),
        'n_glmpca_init': [2],
        'learning_rate': [0.1,0.5,1.],
        'batch_size': [8],
        'max_param': [25],
        'maxiter': [150],
        'gamma': [0.5],
        'step_size': [20]
    }

    grid_search_space['GE'] = {
        'n_pc': np.linspace(2,100,50).astype(int),
        'n_glmpca_init': [2],
        'learning_rate': [0.5, 1., 10.],
        'batch_size': [16],
        'max_param': [25],
        'maxiter': [150],
        'gamma': [0.5],
        'step_size': [20]
    }

    grid_search_space['METH'] = {
        'n_pc': np.linspace(2,150,50).astype(int),
        'n_glmpca_init': [5],
        'learning_rate': [0.05,0.1,0.5],
        'batch_size': [16],
        'max_param': [25],
        'maxiter': [150],
        'gamma': [0.5],
        'step_size': [10]
    }
    grid_search_space['CNA'] = {
        'n_pc': np.linspace(2,100,50).astype(int),
        'n_glmpca_init': [2],
        'learning_rate': [0.5, 1., 10.],
        'batch_size': [16],
        'max_param': [np.inf],
        'maxiter': [150],
        'gamma': [0.5],
        'step_size': [20]
    }
