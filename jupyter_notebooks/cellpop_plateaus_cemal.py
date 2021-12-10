#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 15:12:09 2021

@author: arnab
"""
import libsbml
import importlib
import amici
import amici.plotting
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import pandas as pd
from datetime import datetime
import scipy.stats
from scipy.signal import find_peaks
import argparse
import sys
import os
mpl.rcParams['figure.dpi'] = 300


wd = str(os.getcwd()).replace("jupyter_notebooks","")
sys.path.append(wd+'/bin')
from modules.RunSPARCED import RunSPARCED

#%%

def find_plateaus(F, min_length=200, tolerance = 0.75, smoothing=25):
#     from Rune Højlund on Oct3,2021 (https://stackoverflow.com/questions/53492508/find-plateau-in-numpy-array)
#     Finds plateaus of signal using second derivative of F.
#     Parameters
#     F : Signal.
#     min_length: Minimum length of plateau.
#     tolerance: Number between 0 and 1 indicating how tolerant
#         the requirement of constant slope of the plateau is.
#     smoothing: Size of uniform filter 1D applied to F and its derivatives.
#     Returns
#     plateaus: array of plateau left and right edges pairs
#     dF: (smoothed) derivative of F
#     d2F: (smoothed) Second Derivative of F

    import numpy as np
    from scipy.ndimage.filters import uniform_filter1d
    
    # calculate smooth gradients
    smoothF = uniform_filter1d(F, size = smoothing)
    dF = uniform_filter1d(np.gradient(smoothF),size = smoothing)
    d2F = uniform_filter1d(np.gradient(dF),size = smoothing)
    
    def zero_runs(x):
        # Helper function for finding sequences of 0s in a signal
        # https://stackoverflow.com/questions/24885092/finding-the-consecutive-zeros-in-a-numpy-array/24892274#24892274
        iszero = np.concatenate(([0], np.equal(x, 0).view(np.int8), [0]))
        absdiff = np.abs(np.diff(iszero))
        ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
        return ranges
    
    # Find ranges where second derivative is zero
    # Values under eps are assumed to be zero.
    eps = np.quantile(abs(d2F),tolerance) 
    smalld2F = (abs(d2F) <= eps)
    
    # Find repititions in the mask "smalld2F" (i.e. ranges where d2F is constantly zero)
    p = zero_runs(np.diff(smalld2F))
    
    # np.diff(p) gives the length of each range found.
    # only accept plateaus of min_length
    plateaus = p[(np.diff(p) > min_length).flatten()]
    
    return plateaus #, dF, d2F


#%%
# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name

#%%

sys.path.insert(0, os.path.abspath(model_output_dir))

species_input = os.path.join(wd,'input_files','Species.txt')

species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(species_input, encoding='latin-1')])

species_initializations = []
for row in species_sheet[1:]:
    species_initializations.append(float(row[2]))
species_initializations = np.array(species_initializations)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
species_names = [species_sheet[k][0] for k in range(1,len(species_sheet))]
# for k in range(len(STIMligs)):
#     species_initializations[species_names.index(STIMligs_id[k])] = STIMligs[k]

ts = 30
model_module = importlib.import_module(model_name)
model = model_module.getModel()
solver = model.getSolver()          # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0,ts,2)) # np.linspace(0, 30) # set timepoints

# model.get
# xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],sbml_file,model)

columnsS=[ele for ele in model.getStateIds()]
columnsG = columnsS[773:914]
resa = [sub.replace('m_', 'ag_') for sub in columnsG]
resi = [sub.replace('m_', 'ig_') for sub in columnsG]
columnsG2 = np.concatenate((resa, resi), axis=None)    

#%%

# Create initial conditions for simulations - 48 hr serum+GF starved 

################### User defined #####################
flagD = 0 # deterministic=1, stochastic=0
numofCells = 30
th1 = 24
STIMligs_id = ['E', 'H', 'HGF', 'P', 'F', 'I', 'INS']
STIMligs = [10.0,10.0,10.0,10.0,10.0,10.0,10.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
nmxlsfile = '/app/jupyter_notebooks/SimsOut/MEKi1uM_24st_48gf48inh_10nmLig_'
######################################################

################### Context fitted ###################
k1815val = 2.75e-8 # k1815 value
k1816val = 9e-1  # k1816 value
k1813val = 2.5e-3
k497val = 1e-8 # k497 (kA77) values
k507val = 3.5e-8  # k507 (kA87) values
k337_1val = 7.083333 # Cd degradation (x3 of default)
######################################################

Vn = 1.75E-12
Vc = 5.25E-12
Ve = 5.00E-05
NSteps1 = th1*3600/ts
NSteps1 = int(NSteps1)

speciesInits = np.tile(species_initializations, (1, 1))
counterC = 0

startTime = datetime.now()
print(startTime)

#%% test - runSPARCED
spIn = speciesInits[0]

model.setFixedParameterById('k1815',np.float(k1815val)) 
model.setFixedParameterById('k1816',np.float(k1816val)) 
model.setFixedParameterById('k1813',np.float(k1813val))
model.setFixedParameterById('k497',np.float(k497val))
model.setFixedParameterById('k507',np.float(k507val))
model.setFixedParameterById('k337_1',np.float(k337_1val))

xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th1,spIn,Vn,Vc,model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')



#%%


while counterC < numofCells:
    spIn = speciesInits[0]
    model.setInitialStates(spIn) 
    model.setFixedParameterById('k1815',np.float(k1815val)) 
    model.setFixedParameterById('k1816',np.float(k1816val)) 
    model.setFixedParameterById('k1813',np.float(k1813val))
    model.setFixedParameterById('k497',np.float(k497val))
    model.setFixedParameterById('k507',np.float(k507val))
    model.setFixedParameterById('k337_1',np.float(k337_1val))
    
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th1,spIn,[],sbml_file,model)
    if (len(xoutS_all[:,0])==(NSteps1+1)): 
        speciesInits = np.vstack([speciesInits,xoutS_all[-1,:]])
        
        condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
        condsSDF.to_csv(nmxlsfile+'S_'+str(counterC)+'.txt',sep="\t")  
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
        condsGDF.to_csv(nmxlsfile+'G_'+str(counterC)+'.txt',sep="\t") 
        condsGDF = None
        
        counterC += 1 
    else:
        print('deadcell')
    print(counterC)
    print(datetime.now())
print(datetime.now())

speciesInits = speciesInits[1::]
condsinits = pd.DataFrame(data=speciesInits)
condsinits.to_csv(nmxlsfile+'StarvedInits.txt',sep="\t") 

#%%

th2 = 48

nmxlsfile = '/app/jupyter_notebooks/SimsOut/MEKi1uM_24st48gf_48inh_10nmLig_'

NSteps2 = th2*3600/ts
NSteps2 = int(NSteps2)

cell_lineage = (np.zeros([numofCells,3]))
th2run = np.tile(th2, (numofCells, 1))
speciesInitsInhn = np.tile(species_initializations, (1, 1))
counterC = 0
cellsTot = numofCells

startTime = datetime.now()
print(startTime)

while counterC < cellsTot:
    spIn = speciesInits[counterC]
    for k in range(len(STIMligs)):
        spIn[species_names.index(STIMligs_id[k])] = STIMligs[k]
    
    model_module2 = importlib.import_module(model_name)
    model2 = model_module2.getModel()
    solver2 = model2.getSolver()          
    solver2.setMaxSteps = 1e10
    model2.setTimepoints(np.linspace(0,ts))        
    model2.setFixedParameterById('k1815',np.float(k1815val)) 
    model2.setFixedParameterById('k1816',np.float(k1816val)) 
    model2.setFixedParameterById('k1813',np.float(k1813val))
    model2.setFixedParameterById('k497',np.float(k497val))
    model2.setFixedParameterById('k507',np.float(k507val))
    model2.setFixedParameterById('k337_1',np.float(k337_1val))
    model2.setInitialStates(spIn)
    
    th = th2run[counterC][0]
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,spIn,[],sbml_file,model2)
    cb_peaks, propps = find_peaks(xoutS_all[:, list(species_names).index('Mb')],height=30)
    tt = 0
    if (len(cb_peaks)>0):
        plateaus = find_plateaus(xoutS_all[:, list(species_names).index('Mb')], min_length=120, tolerance = 0.95, smoothing=15)
        plateausT = plateaus[:,0]*30/3600
        plat_idxs = []
        for idx,val in enumerate(cb_peaks):
            if (idx+1==len(cb_peaks)) and (plateausT[-1]<val):
                continue
            else:
                qq = np.nonzero((np.where((val*30.0/3600.0)<plateausT,1,0)))
                plat_idxs.append(qq[0][0])
        for ii in range(len(plat_idxs)):
            tt = th-plateausT[plat_idxs[ii]]
            th2run = np.vstack([th2run,tt])
            ss = xoutS_all[plateaus[plat_idxs[ii],0],:]
            speciesInits = np.vstack([speciesInits,ss])
            cell_lineage = np.vstack([cell_lineage, [counterC+1, 0, 0]])
            cellsTot += 1
            
    condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
    condsSDF.to_csv(nmxlsfile+'S_'+str(counterC)+'.txt',sep="\t")  
    condsSDF = None
    condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
    condsGDF.to_csv(nmxlsfile+'G_'+str(counterC)+'.txt',sep="\t") 
    condsGDF = None
    
    cell_lineage[counterC,1] = len(cb_peaks) # store # of division events for each cell sim.

    if (len(xoutS_all[:,0])==(NSteps2+1)): # store if the cell dies
        cell_lineage[counterC,2] = 0 # did not die
        speciesInitsInhn = np.vstack([speciesInitsInhn,xoutS_all[-1,:]])
    else:
        cell_lineage[counterC,2] = 1 # dead cell
    
    counterC += 1 
    print(counterC)
    print(cellsTot)
    print(datetime.now())
print(datetime.now())

speciesInitsInhn = speciesInitsInhn[1::]
condsinits = pd.DataFrame(data=speciesInitsInhn)
condsinits.to_csv(nmxlsfile+'speciesInitsInhn.txt',sep="\t") 

condsths = pd.DataFrame(data=th2run)
condsths.to_csv(nmxlsfile+'th2run.txt',sep="\t") 

condslineages = pd.DataFrame(data=cell_lineage)
condslineages.to_csv(nmxlsfile+'lineage.txt',sep="\t") 

#%%

th3 = 48
nmxlsfile = '/app/jupyter_notebooks/SimsOut/MEKi1uM_24st48gf48inh_10nmLig_'

NSteps3 = th3*3600/ts
NSteps3 = int(NSteps3)

MEKiIdx = 767
MEKidoses = 1.0 # uM
MEKidose = np.multiply(MEKidoses,1000.0*Ve/Vc) # nM

cellsTot2 = len(speciesInitsInhn)
cell_lineage2 = (np.zeros([cellsTot2,3]))
th2run2 = np.tile(th3, (cellsTot2, 1))
speciesEnd = np.tile(species_initializations, (1, 1))
counterC = 0

startTime = datetime.now()
print(startTime)

while counterC < cellsTot2:
    spIn = speciesInitsInhn[counterC]
    for k in range(len(STIMligs)):
        spIn[species_names.index(STIMligs_id[k])] = STIMligs[k]
    spIn[MEKiIdx] = np.float(MEKidose)
    
    model_module3 = importlib.import_module(model_name)
    model3 = model_module3.getModel()
    solver3 = model3.getSolver()          
    solver3.setMaxSteps = 1e10
    model3.setTimepoints(np.linspace(0,ts))        
    model3.setFixedParameterById('k1815',np.float(k1815val)) 
    model3.setFixedParameterById('k1816',np.float(k1816val)) 
    model3.setFixedParameterById('k1813',np.float(k1813val))
    model3.setFixedParameterById('k497',np.float(k497val))
    model3.setFixedParameterById('k507',np.float(k507val))
    model3.setFixedParameterById('k337_1',np.float(k337_1val))
    model3.setInitialStates(spIn)
    
    th = th2run2[counterC][0]
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,spIn,[],sbml_file,model3)
    cb_peaks, propps = find_peaks(xoutS_all[:, list(species_names).index('Mb')],height=30)
    tt = 0
    if (len(cb_peaks)>0):
        plateaus = find_plateaus(xoutS_all[:, list(species_names).index('Mb')], min_length=120, tolerance = 0.95, smoothing=15)
        plateausT = plateaus[:,0]*30/3600
        plat_idxs = []
        for idx,val in enumerate(cb_peaks):
            if (idx+1==len(cb_peaks)) and (plateausT[-1]<val):
                continue
            else:
                qq = np.nonzero((np.where((val*30.0/3600.0)<plateausT,1,0)))
                plat_idxs.append(qq[0][0])
        for ii in range(len(plat_idxs)):
            tt = th-plateausT[plat_idxs[ii]]
            th2run2 = np.vstack([th2run2,tt])
            ss = xoutS_all[plateaus[plat_idxs[ii],0],:]
            speciesInitsInhn = np.vstack([speciesInitsInhn,ss])
            cell_lineage2 = np.vstack([cell_lineage2, [counterC+1, 0, 0]])
            cellsTot2 += 1
            
    condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
    condsSDF.to_csv(nmxlsfile+'S_'+str(counterC)+'.txt',sep="\t")  
    condsSDF = None
    condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
    condsGDF.to_csv(nmxlsfile+'G_'+str(counterC)+'.txt',sep="\t") 
    condsGDF = None
    
    cell_lineage2[counterC,1] = len(cb_peaks) # store # of division events for each cell sim.

    if (len(xoutS_all[:,0])==(NSteps3+1)): # store if the cell dies
        cell_lineage2[counterC,2] = 0 # did not die
        speciesEnd = np.vstack([speciesInitsInhn,xoutS_all[-1,:]])
    else:
        cell_lineage2[counterC,2] = 1 # dead cell
    
    counterC += 1 
    print(counterC)
    print(cellsTot2)
    print(datetime.now())
print(datetime.now())

speciesEnd = speciesEnd[1::]
condsinits = pd.DataFrame(data=speciesEnd)
condsinits.to_csv(nmxlsfile+'speciesEnd.txt',sep="\t") 

condsths = pd.DataFrame(data=th2run2)
condsths.to_csv(nmxlsfile+'th2run.txt',sep="\t") 

condslineages = pd.DataFrame(data=cell_lineage2)
condslineages.to_csv(nmxlsfile+'lineage.txt',sep="\t") 