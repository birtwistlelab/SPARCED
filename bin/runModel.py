#!/usr/bin/env python3

import libsbml
import importlib
import amici
import amici.plotting
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import scipy.stats
import sys

from modules.SGEmodule import SGEmodule
from modules.RunPrep import RunPrep
from modules.RunSPARCED import RunSPARCED


# input_data_folder = 'input_data/'
# if len(sys.argv) > 1:
#     input_data_folder = sys.argv[1]

# SBML model we want to import
sbml_file = 'SPARCEDv6.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4] # 'BigModel_byparts_v1'
# Directory to which the generated model code is written
model_output_dir = model_name

Vn = 1.75E-12
Vc = 5.25E-12



flagD = 0
th = 2
ts = 30
numStocCells = 1
STIMligs = [100,100.0,100.0,100.0,100.0,100.0,1721.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
flagWr = 1
nmxlsfile = 'GrowthStim_stoc_'

sys.path.insert(0, os.path.abspath(model_output_dir))


species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

species_initializations = []
for row in species_sheet[1:]:
    species_initializations.append(float(row[2]))
species_initializations = np.array(species_initializations)

for nn in range(numStocCells):
    species_initializations[155:162] = STIMligs

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model)

    if flagWr==1:
        columnsS=[ele for ele in model.getStateIds()]
        columnsG = columnsS[773:914]
        resa = [sub.replace('m_', 'ag_') for sub in columnsG]
        resi = [sub.replace('m_', 'ig_') for sub in columnsG]
        columnsG2 = np.concatenate((resa, resi), axis=None)
        condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
        condsSDF.to_excel(nmxlsfile+'S_'+str(nn)+'.xlsx')
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
        condsGDF.to_excel(nmxlsfile+'G_'+str(nn)+'.xlsx')
        condsGDF = None


flagD = 1
th = 2
ts = 30
numStocCells = 1
STIMligs = [100,100.0,100.0,100.0,100.0,100.0,1721.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
flagWr = 1
nmxlsfile = 'GrowthStim_det_'

sys.path.insert(0, os.path.abspath(model_output_dir))

species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0

for nn in range(numStocCells):
    species_initializations[155:162] = STIMligs

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model)

    if flagWr==1:
        columnsS=[ele for ele in model.getStateIds()]
        columnsG = columnsS[773:914]
        resa = [sub.replace('m_', 'ag_') for sub in columnsG]
        resi = [sub.replace('m_', 'ig_') for sub in columnsG]
        columnsG2 = np.concatenate((resa, resi), axis=None)
        condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
        condsSDF.to_excel(nmxlsfile+'S_'+str(nn)+'.xlsx')
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
        condsGDF.to_excel(nmxlsfile+'G_'+str(nn)+'.xlsx')
        condsGDF = None
