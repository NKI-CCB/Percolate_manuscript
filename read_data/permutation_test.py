"""
PERMUTATION TEST

Given a number of components, assess the number of joint components by permutation
testing.
"""

import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from joblib import Parallel, delayed
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from copy import deepcopy
import torch
import umap

sys.path.insert(0, '../../src/genomic_data_reader/')
from read_GDSC_data import read_GDSC_data
from read_GDSC_response import read_GDSC_response
import library_size_normalization

sys.path.insert(0, '../../src/')
# from trickle_down.pred_perf_routine import make_grid_cv
from GLM_JIVE import GLMJIVE
from GLM_JIVE.exponential_family import *
from GLM_JIVE.generalized_SVD import reconstruct_generalized_SVD