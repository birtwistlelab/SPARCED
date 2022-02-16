import pandas as pd
import numpy as np
# import re
# import libsbml
import os
import sys
import importlib
import amici

# from scipy.signal import find_peaks

# import matplotlib as mpl

# from mpi4py import MPI




# mpl.rcParams['figure.dpi'] = 300


wd = str(os.getcwd()).replace("jupyter_notebooks","")




sbml_file = "SPARCED_au565.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.join(wd,model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()

ic_check = pd.Series(index=model.getStateIds(),data=model.getInitialStates())
params_check = pd.Series(index=model.getFixedParameterIds(),data=model.getFixedParameters())

ic_check.to_csv(os.path.join(wd,'ic_check.txt'),sep='\t',index=True,header=False)
params_check.to_csv(os.path.join(wd,'params_check.txt'),sep='\t',index=True,header=False)

#%%

ic_check_palmetto = pd.read_csv(os.path.join(wd,'ic_check_palmetto.txt'),sep='\t',index_col=0,header=None,squeeze=True)
params_check_palmetto = pd.read_csv(os.path.join(wd,'params_check_palmetto.txt'),sep='\t',index_col=0,header=None,squeeze=True)

#%%

ic_check_all = np.zeros(len(ic_check))

for i in range(len(ic_check_all)):
    if abs(ic_check[ic_check.index[i]] - ic_check_palmetto[ic_check.index[i]]) > 1e-6:
        ic_check_all[i] = 1
        
ic_check_all_idx = np.where(ic_check_all)[0]

params_check_all = np.zeros(len(params_check))

for i in range(len(params_check_all)):
    if abs(params_check[params_check.index[i]] - params_check_palmetto[params_check.index[i]]) > 1e-6:
        params_check_all[i] = 1
        
params_check_all_idx = np.where(params_check_all)[0]