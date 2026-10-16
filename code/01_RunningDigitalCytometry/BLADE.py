'''
Python script to run BLADE using Apptainer container on HPC
To include BLADE in evaluation manuscript (during revision)
Leaving input parameters with default values to run deconv.
With harmonized input data as the rest of tools in the study.

Elia Alonso
28-09-2026
'''

print('... setting up')

# import modules for analysis
import sys, os
from Deconvolution.BLADE import Framework
import numpy as np
from numpy import transpose as t
import itertools
import pickle
from scipy.optimize import nnls
from sklearn.svm import SVR
from sklearn.svm import NuSVR
from sklearn.metrics import mean_squared_error as mse
import pandas as pd

# modules for visualization
# import qgrid
from matplotlib import pyplot as plt
import seaborn as sns

# define working dir
WD = "/home/idc/elia.alonso/projects/invitro_mixtures/" # HPC

# define output file to export results
out_path = WD + 'data/output_dig_cyt/BLADE/'

print("... loading input data")

# define path to load input data
inp_path = WD + "data/input_dig_cyt/BLADE/"

# load reference (X, stdX)
X = pd.read_csv(inp_path + "X_SCC040_Fib_TC_NK.csv", index_col=0)
stdX = pd.read_csv(inp_path + "stdX_SCC040_Fib_TC_NK.csv", index_col=0)

# load mixtures (Y)
mix_id = "I_a"
Y = pd.read_csv(inp_path + "Y_" + mix_id + ".csv", index_col=0)

# ----
cells = X.columns
genes = X.index
samples = Y.columns
# ----

# define input arguments

# hyperpars = {
#     'Alpha': [1, 10],                     # leave default values
#     'Alpha0': [0.1, 0.5, 1, 5, 10],       # leave default values
#     'Kappa0': [1, 0.5, 0.1],              # leave default values
#     'SY': np.sqrt([0.1, 0.5, 1, 1.5, 2])  # leave default values
# }

# Nrep = 3        # default values
# Nrepfinal = 10  # default values
# Njob = 50       # default values

# run BLADE
print('... running BLADE')

Y = np.asarray(Y, dtype = float)
X = np.asarray(X, dtype = float)
stdX = np.asarray(stdX, dtype = float)

final_obj, best_obj, best_set, outs = Framework(
    X, stdX, Y
    # Alphas=hyperpars['Alpha'], Alpha0s=hyperpars['Alpha0'], # default settings
    # Kappa0s=hyperpars['Kappa0'], SYs=hyperpars['SY'],       # default settings
    # Nrep=Nrep, Njob=Njob, Nrepfinal=Nrepfinal               # default settings
    )

print("... BLADE execution finished")

# export results
print('... exporting deconvolution results to' + out_path)

## all results from BLADE execution
outfile = out_path + 'BLADE_output_' + mix_id + '.pickle'
pickle.dump(
    {
        'final_obj': final_obj,
        'best_obj': best_obj,
        'best_set': best_set,
        'outs' : outs
    },
        open(outfile, 'wb')
)

## matrices with relevant output for evalaution
obj = final_obj

fraction = t(obj.ExpF(obj.Beta))  # estimated cell fractions
fraction = pd.DataFrame(fraction, index = cells, columns = samples)

GEPs = np.mean(obj.Nu, 0)         # group mode purification
GEPs = pd.DataFrame(GEPs, index = genes, columns = cells)

fraction.to_csv(out_path + 'BLADE_fractions_' + mix_id + '.csv', index=True)
GEPs.to_csv(out_path + 'BLADE_GEPs_' + mix_id + '.csv', index=True)

print("Done.")
