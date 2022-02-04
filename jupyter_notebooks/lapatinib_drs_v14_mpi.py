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

STIMligs = [0.0,0,0,0,0,0,0.0]


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



th = 48

output_dir = output_dose

    

#%% test - comm executor

from mpi4py.futures import MPICommExecutor



def pre_incubate(cell_n):
    print("Running preincubate(%d) on rank %d" %(cell_n+1, MY_RANK))
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    np.savetxt(os.path.join(output_dose,'c'+str(cell_n+1)+'_preinc.txt'),xoutS_all[-1],delimiter='\t')
    
with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
    if executor is not None:
        preinc_begins = time.time()
        print("Submitting pre_incubate stage from rank", MY_RANK)
        res = executor.map(pre_incubate, range(cell_pop), unordered=True)
        num_done = len(list(res))
        preinc_ends = time.time()
        preinc_runtime = float(preinc_ends - preinc_begins)
        print("Finished", num_done,"pre_incubate tasks in",preinc_runtime,"seconds.")
        
#%% broadcast example        
# x = None

# if MY_RANK == 0:
#     x = 123
    
# print("Before broadcast: Rank %d has: x = %s" %(MY_RANK, x))
# x = MPI.COMM_WORLD.bcast(x, root = 0)
# print("After broadcast: Rank %d has: x = %s" %(MY_RANK, x))

#%% initiate gen0


STIMligs = [0.0,0,0,0,0,0,0.0]

dose_egf = float(args.egf)

STIMligs[0] = dose_egf

th = 48

cellpop_g1 = cell_pop

output_dir = output_dose

def cell_g0(cell_n):
    
    cell_name = 'g0_c'+str(cell_n+1)
    
    s_preinc_i = np.loadtxt(os.path.join(output_dir,'c'+str(cell_n+1)+'_preinc.txt'),delimiter='\t')
    sp_input = s_preinc_i
    species_initializations = np.array(sp_input)
    species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
    species_initializations[155:162] = STIMligs
    
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)

    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutS.txt'),xoutS_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutG.txt'),xoutG_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_tout.txt'),tout_all,delimiter='\t')
    
    
with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
    if executor is not None:
        gen0_begins = time.time()
        print("Submitting gen0 stage from rank", MY_RANK)
        res = executor.map(cell_g0, range(cellpop_g1), unordered=True)
        num_done = len(list(res))
        gen0_ends = time.time()
        gen0_runtime = float(gen0_ends - gen0_begins)
        print("Finished", num_done,"gen0 tasks in",gen0_runtime,"seconds.")
        
#%% aux functions
def find_dp(xoutS,tout,species_all=species_all):
    data = xoutS[:,list(species_all).index('Mb')]
    p,_ = find_peaks(data,height=30)
    b = (np.diff(np.sign(np.diff(data))) > 0).nonzero()[0] + 1
    
    if sum(b>p[0]) > 0:
    
        dp = int(b[b>p[0]][0])
        
    else:
        dp = np.nan
    
    return(dp)

def find_dp_all(xoutS,species_all=species_all):
    data = xoutS[:,list(species_all).index('Mb')]
    p,_ = find_peaks(data,height=30)
    b = (np.diff(np.sign(np.diff(data))) > 0).nonzero()[0] + 1
    
    dp_all = []
    for i in range(len(p)):
        b2 = np.where(b>p[i])[0]
        if len(b2)!=0:
            dp_all.append(b[b2[0]])
    
    return(dp_all)

#%% initiate gen1

th = 72


cellpop_g1 = cell_pop

output_dir = output_dose

def cell_g1(cell_n,flagD,th,Vn,Vc,model,wd,omics_input,genereg_input):
    
    cell_name = 'g1_c'+str(cell_n+1)
    
    # s_preinc_i = np.loadtxt(os.path.join(output_dose,'c'+str(cell_n+1)+'_preinc.txt'),delimiter='\t')
    x_s_g0 = np.loadtxt(os.path.join(output_dose,'g0_c'+str(cell_n+1)+'_xoutS.txt'),delimiter='\t')
    np.random.seed()
    tp_g0 = np.random.randint(0,np.shape(x_s_g0)[0])    
    sp_input = x_s_g0[tp_g0,:]
    species_initializations = np.array(sp_input)
    species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
    # species_initializations[155:162] = STIMligs
    
    species_initializations[list(model.getStateIds()).index('lapatinib')] = dose
    
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    
    np.savetxt(os.path.join(output_dir,'g0_c'+str(int(cell_n+1))+'_tp.txt'),np.array([tp_g0]),delimiter='\t')

    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutS.txt'),xoutS_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutG.txt'),xoutG_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_tout.txt'),tout_all,delimiter='\t')
