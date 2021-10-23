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


th = 24

def pre_incubate(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input):
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    return (xoutS_all, xoutG_all, tout_all)

# args_pool = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input]
start = time.perf_counter()

k, l, m = pre_incubate(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)

finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} seconds(s)')


#%% test - pre-incubate/multi process

import time
import multiprocessing

th = 24

def pre_incubate(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input):
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    return (xoutS_all, xoutG_all, tout_all)

# args_pool = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input]
start = time.perf_counter()

# k, l, m = pre_incubate(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)

# p1 = multiprocessing.Process(target=pre_incubate, args = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])
# p2 = multiprocessing.Process(target=pre_incubate, args = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])

# p1.start()
# p2.start()

# p1.join()
# p2.join()


processes = []

for _ in range(5):
    p = multiprocessing.Process(target=pre_incubate, args = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()

finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} seconds(s)')

#%% pre-incubate / process pool

import concurrent.futures
import time

def pre_incubate(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input):
    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    return (xoutS_all, xoutG_all, tout_all)

start = time.perf_counter()

with concurrent.futures.ProcessPoolExecutor() as executor:
    f1 = executor.submit(pre_incubate, [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])
    # f2 = executor.submit(pre_incubate, [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])

k,l,m = f1.results()
finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} seconds(s)')


#%%
import concurrent.futures

cell_n = 5

args_pool = [flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input]

with concurrent.futures.ProcessPoolExecutor() as executor:
    # f1 = executor.submit(do_something, 1)
    # f2 = executor.submit(do_something, 1)
    # print(f1.result())
    # print(f2.result())
    results = [executor.submit(pre_incubate,args_pool) for _ in range(2)]
    
    # for f in concurrent.futures.as_completed(results):
    #     print(f.result())

# run process pool and view results


#%% test - multi processing, concurrent futures
import multiprocessing
import time

start = time.perf_counter()

def do_something(seconds):
    print(f'Sleeping {seconds} second(s)...')
    time.sleep(seconds)
    print('Done Sleeping...')
    
# p1 = multiprocessing.Process(target=do_something)
# p2 = multiprocessing.Process(target=do_something)

# p1.start()
# p2.start()

# p1.join()
# p2.join()

processes = []

for _ in range(1):
    p = multiprocessing.Process(target=do_something, args = [1])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()


finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} seconds(s)')

#%% test - process pool executor
import concurrent.futures
import multiprocessing
import time

start = time.perf_counter()

def do_something(seconds):
    print(f'Sleeping {seconds} second(s)...')
    time.sleep(seconds)
    return 'Done Sleeping...'
    
# p1 = multiprocessing.Process(target=do_something)
# p2 = multiprocessing.Process(target=do_something)

# p1.start()
# p2.start()

# p1.join()
# p2.join()

with concurrent.futures.ProcessPoolExecutor() as executor:
    # f1 = executor.submit(do_something, 1)
    # f2 = executor.submit(do_something, 1)
    # print(f1.result())
    # print(f2.result())
    results = [executor.submit(do_something,1) for _ in range(10)]
    
    for f in concurrent.futures.as_completed(results):
        print(f.result())

# processes = []

# for _ in range(10):
#     p = multiprocessing.Process(target=do_something, args = [1.5])
#     p.start()
#     processes.append(p)
    
# for process in processes:
#     process.join()


finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} seconds(s)')



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