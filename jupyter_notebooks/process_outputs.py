#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 11:56:26 2021

@author: arnab
"""

import pandas as pd
import numpy as np
import re
import libsbml
import os
import sys
import importlib
import amici
# import amici.plotting
# import argparse
# import multiprocessing as mpr

#%%
wd = str(os.getcwd()).replace("jupyter_notebooks","")


sim_name = 'lapatinib_dose_response'

output_path = os.path.join('/media/arnab/bigchungus/projects/ccle_laptinib_drs/sparced',sim_name)

# if not os.path.exists(output_path):
#     os.mkdir(output_path)


sbml_file = "SPARCED_au565.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.join(wd,model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()

#%%

species_all = list(model.getStateIds())

n_cells = 100

dose = 0.0025

#%%

# xoutS_read = np.loadtxt(os.path.join(output_path,'lapatinib_0.0025','c2_xoutS_all.txt'),delimiter='\t')

parp_idx = species_all.index('PARP')
cparp_idx = species_all.index('cPARP')

# xoutS_read[-1,cparp_idx]

cells_living = 0
cells_dead = 0

for c in range(n_cells):
    xoutS_c = np.loadtxt(os.path.join(output_path,'lapatinib_'+str(dose),'c'+str(c)+'_xoutS_all.txt'),delimiter='\t')
    parp = float(xoutS_c[-1,parp_idx])
    cparp = float(xoutS_c[-1,cparp_idx])
    
    if parp > cparp:
        cells_living = cells_living + 1
    elif cparp > parp:
        cells_dead = cells_dead + 1
        
#%%

doses = [0.0025, 0.008, 0.025, 0.08, 0.25, 0.8, 2.53, 8.0]
response = pd.DataFrame(data=doses)

living = []
dead = []

for dose in doses:
    cells_living = 0
    cells_dead = 0
    
    for c in range(n_cells):
        xoutS_c = np.loadtxt(os.path.join(output_path,'lapatinib_'+str(dose),'c'+str(c)+'_xoutS_all.txt'),delimiter='\t')
        parp = float(xoutS_c[-1,parp_idx])
        cparp = float(xoutS_c[-1,cparp_idx])
        
        if parp > cparp:
            cells_living = cells_living + 1
        elif cparp > parp:
            cells_dead = cells_dead + 1
        
    living.append(cells_living)
    dead.append(cells_dead)

#%%
response['living'] = np.array(living)
response['dead'] = np.array(dead)

response['viability'] = response['living']/n_cells

response.to_csv('lapatinib_response.txt',sep='\t')

#%%

doses = [0.0025, 0.008, 0.025, 0.08, 0.25, 0.8, 2.53, 8.0]

output_path = os.path.join('/scratch1/amutsud/modeling/ccle/sparced/',sim_name)

def single_cell(c):
    xoutS_c = np.loadtxt(os.path.join(output_path,'lapatinib_'+str(dose),'c'+str(c)+'_xoutS_all.txt'),delimiter='\t')
    parp = float(xoutS_c[-1,parp_idx])
    cparp = float(xoutS_c[-1,cparp_idx])
    
    if parp > cparp:
        cells_living = 1
        # cells_dead = 0
    elif cparp > parp:
        cells_living = 0
        # cells_dead = 1
    return (cells_living)

#%%

import concurrent.futures

living = []

for dose in doses:

    with concurrent.futures.ProcessPoolExecutor() as executor:
        cells = range(n_cells)
        results = executor.map(single_cell,cells)
        
        living_cells = []
        for result in results():
            living_cells.append(int(result))
            
    living_cells_n = sum(living_cells)
    living.append(living_cells_n)
    
response = pd.DataFrame(data=doses)
response['living'] = living
response['viability'] = response['living']/n_cells

response.to_csv('lapatinib_response.txt',sep='\t')
    
    

