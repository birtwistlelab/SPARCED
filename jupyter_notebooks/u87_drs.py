#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 22:18:22 2021

@author: arnab
"""
##

import pandas as pd
import numpy as np
import re
import libsbml
import os
import sys
import importlib
import amici
import amici.plotting
import argparse
import multiprocessing as mpr

parser = argparse.ArgumentParser(description='Input doses in uM')
parser.add_argument('--dose', metavar='dose', help='input dose in uM', default = 0.0025)
args = parser.parse_args()

wd = str(os.getcwd()).replace("jupyter_notebooks","")


sim_name = 'u87_drs_test'

output_path = os.path.join(wd,sim_name)

if not os.path.exists(output_path):
    os.mkdir(output_path)


sbml_file = "SPARCED_u87i.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.join(wd,model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()

#%%
sys.path.append(wd+'/bin')
from modules.RunSPARCED import RunSPARCED


omics_input = 'OmicsData_extended_u87.txt'
genereg_input = 'GeneReg.txt'

flagD = 1

# deterministic='GrowthStim_det_', stochastic='GrowthStim_stoc_'
# nmxlsfile = 'U87SPARCED_1nmGMDet_'

ts = 30
th = 72
Vn = 1.75E-12
Vc = 5.25E-12
STIMligs = [10.0,0,0,0,0,0,10.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS

#%% dose

## ccle doese (um): 0.0025, 0.008, 0.025, 0.08, 0.25, 0.80, 2.53, 8

# dose = float(args.dose) * 10e2

dose = 0.0

sp_input = pd.read_csv(os.path.join(wd,'initializer','species_u87i.txt'),sep='\t',header=None,index_col=0,squeeze=True)
sp_input['lapatinib'] = dose

output_dose = os.path.join(output_path,'lapatinib_'+str(float(args.dose)))

if not os.path.exists(output_dose):
    os.mkdir(output_dose)



#%%

species_initializations = np.array(sp_input)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
species_initializations[155:162] = STIMligs


solver = model.getSolver()          # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

cell_pop = 100

#%%


# for c in range(cell_pop):
    
#     xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
#     np.savetxt(os.path.join(output_dose,'c'+str(c)+'_xoutS_all.txt',xoutS_all,delimiter='\t'))
#     np.savetxt(os.path.join(output_dose,'c'+str(c)+'_xoutG_all.txt',xoutG_all,delimiter='\t'))
#     np.savetxt(os.path.join(output_dose,'c'+str(c)+'_tout_all.txt',tout_all,delimiter='\t'))

#%% test - pre-incubation

import time

flagD = 0

th = 24

def pre_incubate(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input):
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    return (xoutS_all, xoutG_all, tout_all)

# args_pool = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input]
start = time.perf_counter()

k, l, m = pre_incubate(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)

finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} seconds(s)')



#%% test - queue

STIMligs = [0.0,0,0,0,0,0,0.0]


dose = 0.0

sp_input = pd.read_csv(os.path.join(wd,'initializer','species_u87i.txt'),sep='\t',header=None,index_col=0,squeeze=True)
sp_input['lapatinib'] = dose

output_dose = os.path.join(output_path,'lapatinib_'+str(float(dose)))

if not os.path.exists(output_dose):
    os.mkdir(output_dose)





species_initializations = np.array(sp_input)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
species_initializations[155:162] = STIMligs


solver = model.getSolver()          # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

cell_pop = 5


import time
import multiprocessing

th = 24



def pre_incubate(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input,s_q,g_q):
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    s_q.put(xoutS_all[-1])
    g_q.put(xoutG_all[-1])
    # return (xoutS_all, xoutG_all, tout_all)

flagD=0
# args_pool = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input]
start = time.perf_counter()

# k, l, m = pre_incubate(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)

# p1 = multiprocessing.Process(target=pre_incubate, args = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])
# p2 = multiprocessing.Process(target=pre_incubate, args = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])

# p1.start()
# p2.start()

# p1.join()
# p2.join()
if __name__ == "__main__":

    s_q = multiprocessing.Queue()
    g_q = multiprocessing.Queue()
    
    s_all = []
    g_all = []
    
    processes = []
    
    for _ in range(cell_pop):
        p = multiprocessing.Process(target=pre_incubate, args = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input,s_q,g_q])
        p.start()
        processes.append(p)
        
    for process in processes:
        process.join()
    
    while not s_q.empty():
        s_all.append(s_q.get())
        
    while not g_q.empty():
        g_all.append(g_q.get())

finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} seconds(s)')

#%% save preincubation result
s_preinc = s_all
np.savetxt(os.path.join(output_dose,'pre_inc_'+str(cell_pop)+'.txt'), s_preinc, delimiter = '\t')

#%%

s_preinc = np.loadtxt(os.path.join(output_dose,'pre_inc_5.txt'), delimiter = '\t')
#%%
import itertools

STIMligs = [10.0,0,0,0,0,0,10.0]
sp_input = pd.read_csv(os.path.join(wd,'initializer','species_u87i.txt'),sep='\t',header=None,index_col=0,squeeze=True)

th = 6


species_initializations = np.array(sp_input)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
species_initializations[155:162] = STIMligs


solver = model.getSolver()          # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints



def single_cell(s_preinc_i,STIMligs,cell_n,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input,sq):
    
    sp_input = s_preinc_i
    
    species_initializations = np.array(sp_input)
    species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
    species_initializations[155:162] = STIMligs
    
    
    
    
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    
    xoutS_lite = list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),100))
    # xoutG_lite = itertools.islice(xoutG_all,0,(len(xoutG_all)-1),100)
    # tout_lite = itertools.islice(tout_all,0,(len(tout_all)-1),100)
    
    sq.put(xoutS_lite)
    # gq.put(xoutG_lite)
    # tq.put(tout_lite)
    
    
    # np.savetxt(os.path.join(output_dose,'c'+str(cell_n)+'_xoutS.txt'),xoutS_lite,delimiter='\t')
    # np.savetxt(os.path.join(output_dose,'c'+str(cell_n)+'_xoutG.txt'),xoutG_lite,delimiter='\t')
    # np.savetxt(os.path.join(output_dose,'c'+str(cell_n)+'_tout.txt'),tout_lite,delimiter='\t')


#%% u87 dose 72 hrs - q
cell_pop = 4

sq = multiprocessing.Queue()
# gq = multiprocessing.Queue()
# tq = multiprocessing.Queue()


s_all = []
# g_all = []

processes = []

for c in range(cell_pop):
    p = multiprocessing.Process(target=single_cell, args = [s_preinc[c],STIMligs,c,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input,sq])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()

while not s_q.empty():
    s_all.append(s_q.get())
    
# while not g_q.empty():
#     g_all.append(g_q.get())


#%%
# np.savetxt('xoutS_text.txt',xoutS_all,delimiter='\t')
# xoutS_read = np.loadtxt('xoutS_text.txt',delimiter='\t')

def single_cell(cell_n,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input):
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    np.savetxt(os.path.join(output_dose,'c'+str(cell_n)+'_xoutS_all.txt'),xoutS_all,delimiter='\t')
    np.savetxt(os.path.join(output_dose,'c'+str(cell_n)+'_xoutG_all.txt'),xoutG_all,delimiter='\t')
    np.savetxt(os.path.join(output_dose,'c'+str(cell_n)+'_tout_all.txt'),tout_all,delimiter='\t')
    
# p1 = mpr.Process(target=single_cell)

# p1.start()

# p1.join()

processes = []

for c in range(cell_pop):
    cell_n = c
    p = mpr.Process(target=single_cell,args=[cell_n,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])
    p.start()
    processes.append(p)
    
    
for process in processes:
    process.join()
    
