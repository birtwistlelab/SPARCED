import sys
import os
sys.path.append(os.getcwd()[0:os.getcwd().rfind('/')]+'/bin')
#
import libsbml
import importlib
import amici
import amici.plotting
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import scipy.stats
import argparse

from modules.RunSPARCED import RunSPARCED

# SBML model we want to import
sbml_file = 'SPARCED_U87.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name
# used for filename to save
cellNumber = 0  

# deterministic=1, stochastic=0
flagD = 1

# deterministic='GrowthStim_det_', stochastic='GrowthStim_stoc_'
nmxlsfile = 'U87SPARCED_1nmGMDet_'

ts = 30
th = 72
Vn = 1.75E-12
Vc = 5.25E-12
STIMligs = [1.0,1.0,1.0,1.0,1.0,1.0,1.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS

#%%

sys.path.insert(0, os.path.abspath(model_output_dir))
species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

species_initializations = []
for row in species_sheet[1:]:
    species_initializations.append(float(row[2]))
species_initializations = np.array(species_initializations)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
species_initializations[155:162] = STIMligs

#%%
model_module = importlib.import_module(model_name)
model = model_module.getModel()
solver = model.getSolver()          # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model)

# saves output
columnsS=[ele for ele in model.getStateIds()]
columnsG = columnsS[773:914]
resa = [sub.replace('m_', 'ag_') for sub in columnsG]
resi = [sub.replace('m_', 'ig_') for sub in columnsG]
columnsG2 = np.concatenate((resa, resi), axis=None)
condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
condsSDF.to_csv(nmxlsfile+'S_'+str(cellNumber)+'.txt',sep="\t")  
condsSDF = None
condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
condsGDF.to_csv(nmxlsfile+'G_'+str(cellNumber)+'.txt',sep="\t") 
condsGDF = None

#%%

# Import new starting point (U87_UoE cells / stst with GM, timepoint=0) and run 24hr Starve, 48hr GM, 72hr MEKi+GM
flagD = 0
ts = 30
Vn = 1.75E-12
Vc = 5.25E-12
Ve = 5.00E-05
th1 = 24
th2 = 48
th3 = 72
NSteps1 = th1*3600/ts
NSteps1 = int(NSteps1)
NSteps2 = th2*3600/ts
NSteps2 = int(NSteps2)
NSteps3 = th3*3600/ts
NSteps3 = int(NSteps3)

numofCells = 100
nmxlsfile = 'MEKi037uM72hrSims_1nmGM_'

STIMligs = [1.0,1.0,1.0,1.0,1.0,1.0,1.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
MEKiIdx = 767
k1813val = 2.5e-3
k497val = 1e-8 # k497 (kA77) values
k507val = 3.5e-8  # k507 (kA87) values
k1815val = 2.75e-8 # k1815 value
k1816val = 9e-1  # k1816 value
MEKidoses = 0.037037 # uM
MEKidose = np.multiply(MEKidoses,1000.0*Ve/Vc) # nM

sys.path.insert(0, os.path.abspath(model_output_dir))
species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

spdata = []
for row in species_sheet[1:]:
    spdata.append(float(row[2]))
spdata = np.array(spdata)
spdata[np.argwhere(spdata <= 1e-6)] = 0.0

startTime = datetime.now()
print(startTime)

for nbr in range(0,numofCells):
    spIn = spdata
    spIn[155:162] = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver()          
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts))    
    model.setInitialStates(spIn)
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th1,spIn,genedata,Vn,Vc,model)
    lenSim, lenSps = xoutS_all.shape
    print([nbr, 1, lenSim])
    
    if lenSim>NSteps1:
        spIn = xoutS_all[-1,:]
        spIn[155:162] = STIMligs
        model_module2 = importlib.import_module(model_name)
        model2 = model_module2.getModel()
        solver2 = model2.getSolver()          
        solver2.setMaxSteps = 1e10
        model2.setTimepoints(np.linspace(0,ts))        
        model2.setInitialStates(spIn)
        xoutS_all2, xoutG_all2, tout_all2 = RunSPARCED(flagD,th2,spIn,xoutG_all[-1,:],Vn,Vc,model2)
        lenSim2, lenSps2 = xoutS_all2.shape
    
        if lenSim2>NSteps2:
            spIn = xoutS_all2[-1,:]
            spIn[155:162] = STIMligs
            spIn[MEKiIdx] = np.float(MEKidose)
            model_module3 = importlib.import_module(model_name)
            model3 = model_module3.getModel()
            solver3 = model3.getSolver()          
            solver3.setMaxSteps = 1e10
            model3.setTimepoints(np.linspace(0,ts))        
            model3.setInitialStates(spIn)
            xoutS_all3, xoutG_all3, tout_all3 = RunSPARCED(flagD,th3,spIn,xoutG_all2[-1,:],Vn,Vc,model3)
        
    xoutS = np.concatenate((xoutS_all, xoutS_all2, xoutS_all3), axis=0)
    xoutG = np.concatenate((xoutG_all, xoutG_all2, xoutG_all3), axis=0)
    lenSim3, lenSps3  = xoutS.shape
    print([nbr,2,lenSim3])
    
    # saves output
    columnsS=[ele for ele in model.getStateIds()]
    columnsG = columnsS[773:914]
    resa = [sub.replace('m_', 'ag_') for sub in columnsG]
    resi = [sub.replace('m_', 'ig_') for sub in columnsG]
    columnsG2 = np.concatenate((resa, resi), axis=None)
    condsSDF = pd.DataFrame(data=xoutS,columns=columnsS)
    condsSDF.to_csv(nmxlsfile+'S_'+str(nbr)+'.txt',sep="\t")  
    condsSDF = None
    condsGDF = pd.DataFrame(data=xoutG,columns=columnsG2)
    condsGDF.to_csv(nmxlsfile+'G_'+str(nbr)+'.txt',sep="\t") 
    condsGDF = None
    print(datetime.now())
print(datetime.now())

#%%

# cell pop test

# Deterministic Run

flagD = 1
solver = model.getSolver()  # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints

#%%