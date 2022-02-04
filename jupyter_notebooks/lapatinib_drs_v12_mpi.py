#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 22:18:22 2021

@author: arnab
"""
##

import pandas as pd
import numpy as np
# import re
import libsbml
import os
import sys
import importlib
import amici
import argparse
from scipy.signal import find_peaks
import itertools
import time
import multiprocessing
import concurrent.futures
import matplotlib.pyplot as plt
import matplotlib as mpl

from mpi4py import MPI

MY_RANK = MPI.COMM_WORLD.Get_rank()


mpl.rcParams['figure.dpi'] = 300

parser = argparse.ArgumentParser(description='Input doses in uM')
parser.add_argument('--dose', metavar='dose', help='input laptinib dose in uM', default = 0.0)
parser.add_argument('--egf', metavar='egf', help='input EGF conc in nM', default = 100.0)
parser.add_argument('--cellpop', metavar='cellpop', help='starting cell population', default = 5)
parser.add_argument('--td',metavar='td', help='cell line doubling time (hrs) ', default = 48)
parser.add_argument('--sim_name',metavar='sim_name', help='insert exp name', default = 'lapatinib_drs_test_mpi')
args = parser.parse_args()

wd = str(os.getcwd()).replace("jupyter_notebooks","")


sim_name = str(args.sim_name)

output_path = os.path.join(wd,'output',sim_name)

if MY_RANK==0:
    if not os.path.exists(output_path):
        os.mkdir(output_path)


sbml_file = "SPARCED_au565.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.join(wd,model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()

#%%
sys.path.append(wd+'/bin')
from modules.RunSPARCED import RunSPARCED


omics_input = 'OmicsData_extended_au565.txt'
genereg_input = 'GeneReg_au565.txt'

flagD = 0

# deterministic='GrowthStim_det_', stochastic='GrowthStim_stoc_'
# nmxlsfile = 'U87SPARCED_1nmGMDet_'

ts = 30
# th = 24
Vn = 1.75E-12
Vc = 5.25E-12
# STIMligs = [0.0,0,0,0,0,0,0.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS




#%%

species_all = list(model.getStateIds())

solver = model.getSolver()          # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

cell_pop = int(args.cellpop)


#%% preincubate

STIMligs = [100.0,0,0,0,0,0,0.0]


dose = float(args.dose)*10e2

sp_input = pd.read_csv(os.path.join(wd,'initializer','species_au565.txt'),sep='\t',header=None,index_col=0,squeeze=True)
# sp_input['lapatinib'] = dose


species_initializations = np.array(sp_input)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
species_initializations[155:162] = STIMligs


output_dose = os.path.join(output_path,'lapatinib_'+str(float(args.dose)))

if MY_RANK==0:
    if not os.path.exists(output_dose):
        os.mkdir(output_dose)

# np.linspace(0, 30) # set timepoint




th = 48

output_dir = output_dose

#%%

# def pre_incubate(cell_n,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input):
#     xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    
#     np.savetxt(os.path.join(output_dose,'c'+str(cell_n+1)+'_preinc.txt'),xoutS_all[-1],delimiter='\t')
    
    
# start = time.perf_counter()




#%% test - pool executor

# from mpi4py.futures import MPIPoolExecutor

# def pre_incubate(cell_n):
#     xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    
#     np.savetxt(os.path.join(output_dose,'c'+str(cell_n+1)+'_preinc.txt'),xoutS_all[-1],delimiter='\t')

# if __name__ == "__main__":
    
#     # flagD = 1


    
#     # processes = []
    
#     # for c in range(cell_pop):
#     #     p = multiprocessing.Process(target=pre_incubate, args = [c,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])
#     #     p.start()
#     #     processes.append(p)
        
#     # for process in processes:
#     #     process.join()
    
#     with MPIPoolExecutor() as executor:
#         res = executor.map(pre_incubate, range(cell_pop), unordered=True)
#         num_done = len(list(res))
#         print("Finished", num_done, "pre_incubate tasks.")
    

#%% test - comm executor

from mpi4py.futures import MPICommExecutor



def pre_incubate(cell_n):
    print("Running preincubate(%d) on rank %d" %(cell_n+1, MY_RANK))
    # arr = ["This is rank",str(MY_RANK)]
    # np.savetxt(os.path.join(output_dose,str("cell"+str(cell_n)+"_preinc.txt")),np.array(arr),delimiter='\t',fmt='%s')
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    np.savetxt(os.path.join(output_dose,str("g0c"+str(cell_n+1)+"_xoutS.txt")),xoutS_all,delimiter='\t')
    np.savetxt(os.path.join(output_dose,str("g0c"+str(cell_n+1)+"_tout.txt")),tout_all,delimiter='\t')
    
with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
    if executor is not None:
        print("Submitting pre_incubate stage from rank", MY_RANK)
        res = executor.map(pre_incubate, range(cell_pop), unordered=True)
        num_done = len(list(res))
        print("Finished", num_done,"pre_incubate tasks.")
        
        
x = None

if MY_RANK == 0:
    x = 123
    
print("Before broadcast: Rank %d has: x = %s" %(MY_RANK, x))
x = MPI.COMM_WORLD.bcast(x, root = 0)
print("After broadcast: Rank %d has: x = %s" %(MY_RANK, x))