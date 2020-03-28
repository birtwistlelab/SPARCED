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

from SGEmodule import SGEmodule
from RunPrep import RunPrep
from RunSPARCED import RunSPARCED


# SBML model we want to import
sbml_file = 'SPARCEDv6.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4] # 'BigModel_byparts_v1'
# Directory to which the generated model code is written
model_output_dir = model_name

Vn = 1.75E-12
Vc = 5.25E-12

currentDirectory = os.getcwd()


flagD = 0
th = 2
ts = 30
numStocCells = 1
STIMligs = [100,100.0,100.0,100.0,100.0,100.0,1721.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
flagWr = 1
nmxlsfile = 'GrowthStim_stoc_'

sys.path.insert(0, os.path.abspath(model_output_dir))

spdata0 = pd.read_excel('Species_v6.xlsx',header=0,index_col=0)
spdata = np.double(spdata0.values[:,1])
spdata[np.argwhere(spdata <= 1e-6)] = 0.0

startTime = datetime.now()
print(startTime)
for nn in range(numStocCells):
    spdata[155:162] = STIMligs

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,spdata,[],Vn,Vc,model)

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
    print(datetime.now() - startTime)
print(datetime.now())


flagD = 1
th = 2
ts = 30
numStocCells = 1
STIMligs = [100,100.0,100.0,100.0,100.0,100.0,1721.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
flagWr = 1
nmxlsfile = 'GrowthStim_det_'

sys.path.insert(0, os.path.abspath(model_output_dir))

spdata0 = pd.read_excel('Species_v6.xlsx',header=0,index_col=0)
spdata = np.double(spdata0.values[:,1])
spdata[np.argwhere(spdata <= 1e-6)] = 0.0

startTime = datetime.now()
print(startTime)
for nn in range(numStocCells):
    spdata[155:162] = STIMligs

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,spdata,[],Vn,Vc,model)

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
    print(datetime.now() - startTime)
print(datetime.now())
