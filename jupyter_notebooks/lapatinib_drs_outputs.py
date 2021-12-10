#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 13:39:25 2021

@author: arnab
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import importlib
import random

mpl.rcParams['figure.dpi'] = 300

wd = str(os.getcwd()).replace("jupyter_notebooks","")

wd_output = '/media/arnab/bigchungus/projects/ccle_laptinib_drs/sparced'


sim_name = 'lapatinib_dose_response'

#%% import model
sbml_file = "SPARCED_au565.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.join(wd,model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()

#%%

dose_output_all = os.listdir(os.path.join(wd_output,sim_name))

dose_output_all = dose_output_all.sort()

cell_pop = 100

dose_resp = []

for dose in dose_output_all:
    
    
    flagA_dose = []
    
    for c in range(cell_pop):
        flagA_c = float(np.loadtxt(os.path.join(wd_output,sim_name,dose,'c'+str(c)+'_flagA.tx')))
        flagA_dose.append(flagA_c)
        
    flagA_dose = np.array(flagA_dose)
    
    dose_resp.append(sum(flagA_dose))
    
#%%

xoutS_cellpop = []
for c in range(cell_pop):
    xoutS_cellpop.append(np.loadtxt(os.path.join(wd_output,sim_name,dose,'c'+str(c)+'_xoutS.txt')))

plt.plot(xoutS_cellpop[2][:,list(model.getStateIds()).index('ppERK')])


#%%
tout_lite = np.loadtxt(os.path.join(wd_output,sim_name,dose,'c'+str(c)+'_tout.txt'))

xoutS_lite = np.loadtxt(os.path.join(wd_output,sim_name,dose,'c'+str(c)+'_xoutS.txt'))


#%%
def timecourse(species,tout,dose,y_ax_lim='default'): 
    
    x_t_pop = []
    
    for r in range(10):
        
        c = random.randint(0,99)
    
        timeh = tout/3600
        species_ind = list(model.getStateIds()).index(species)
        xoutS = np.loadtxt(os.path.join(wd_output,sim_name,dose,'c'+str(c)+'_xoutS.txt'))
        
        # xoutS_pop.append(xoutS)
    
    
         
        x_t = xoutS[:,species_ind]
        x_t_pop.append(x_t)
        plt.plot(timeh,x_t)
        
    x_t_pop = np.array(x_t_pop)
    
    # for k in range(np.shape(xoutS_pop)[0]):
        # x_t = xoutS[k,species_ind]
        # plt.plot(timeh,x_t)
        
    
    # c = random.randint(0,99)
    
    # timeh = tout/3600
    # species_ind = list(model.getStateIds()).index(species)
    # xoutS = np.loadtxt(os.path.join(wd,sim_name,dose,'c'+str(c)+'_xoutS.txt'))


     
    # x_t = xoutS[:,species_ind]
    # plt.plot(timeh,x_t)
    plt.ylabel(str(species))
    if type(y_ax_lim) == 'str': 
        plt.ylim(0,1.2*max(x_t_pop.max(axis=1)))
    elif type(y_ax_lim) == int or type(y_ax_lim) == float:
        plt.ylim(0,y_ax_lim)
    plt.xlabel('time(h)')
    plt.title(str(dose)+str('_um'))
    plt.show
    # return(x_t_pop)

#%%