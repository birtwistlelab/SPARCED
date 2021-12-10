#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 14:59:16 2021

@author: arnab
"""
import re
import matplotlib as mpl

import argparse
import scipy.stats
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import libsbml
import importlib
import amici
import amici.plotting
import os
import sys

mpl.rcParams['figure.dpi'] = 300


wd = str(os.getcwd()).replace("jupyter_notebooks","")
sys.path.append(wd+'/bin')

from modules.RunSPARCED import RunSPARCED



# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name
sys.path.insert(0, os.path.join(wd,model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()

#%%

# th = 48

Vn = float(1.7500E-12)
Vc = float(5.2500E-12)
STIMligs = [100.0,0.0,0.0,0.0,0.0,0.0,100.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
STIMligs_id = ['E', 'H', 'HGF', 'P', 'F', 'I', 'INS']


#%%

species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(
    os.path.join(wd,'input_files','Species.txt'), encoding='latin-1')])


species_all = [species_sheet[k][0] for k in range(1,len(species_sheet))]

species_initializations = []
for row in species_sheet[1:]:
    species_initializations.append(float(row[2]))
species_initializations = np.array(species_initializations)



for k in range(len(STIMligs)):
    species_initializations[species_all.index(STIMligs_id[k])] = STIMligs[k]

#%% stochastic run

ts = 30
flagD = 0
solver = model.getSolver()  # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints

# th = 36

# xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD, th, species_initializations, Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')

#%%
from scipy.signal import find_peaks
import itertools


cellpop_g1 = 5

th = 72

flagD = 0
th_g2 = []
sp_g2 = []

cellpop_g2 = 0

for cell in range(cellpop_g1):
    cell_name = 'g1_c'+str(cell+1)
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD, th, species_initializations, Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')
    
    xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),20)))
    xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),20)))
    tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),20)))
    
    np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_xoutS.txt'),xoutS_lite,delimiter='\t')
    np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_xoutG.txt'),xoutG_lite,delimiter='\t')
    np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_tout.txt'),tout_lite,delimiter='\t')
    
    cb_peaks, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
    
    if len(cb_peaks)>0:
    
        th_g2_cell = tout_all[cb_peaks[0]]/3600
        sp_g2_cell = xoutS_all[cb_peaks[0]]
    
        th_g2.append(th_g2_cell)
        th_g2.append(th_g2_cell)
        sp_g2.append(sp_g2_cell)
        sp_g2.append(sp_g2_cell)
        cellpop_g2 = cellpop_g2 + 2
    
# cellpop_total = cellpop_g1 + cellpop_g2

#%%


# mpl.rcParams['figure.dpi'] = 300

def timecourse(species, x_s=xoutS_all, tout_all=tout_all):

    x_t = x_s[:, list(species_all).index(species)]
    plt.plot(tout_all/3600, x_t)
    plt.ylabel('species: '+str(species))
    plt.xlabel('time(h)')
    plt.ylim(0, max(x_t)*1.25)

    plt.show
    
#%%

tout = np.loadtxt(os.path.join(wd,'output','cellpop_test','g1_c1_tout.txt'),delimiter='\t')
xoutS1 = np.loadtxt(os.path.join(wd,'output','cellpop_test','g1_c1_xoutS.txt'),delimiter='\t')
xoutS2 = np.loadtxt(os.path.join(wd,'output','cellpop_test','g1_c2_xoutS.txt'),delimiter='\t')
xoutS3 = np.loadtxt(os.path.join(wd,'output','cellpop_test','g1_c3_xoutS.txt'),delimiter='\t')


timecourse('Mb',xoutS1,tout)
timecourse('Mb',xoutS2,tout)
timecourse('Mb',xoutS3,tout)

#%%

output_dir = os.path.join(wd,'output','cellpop_test')

# sp_gn0 = np.array([item for sublist in sp_g2 for item in sublist])
sp_gn0 = np.array(sp_g2)

# th_gn0 = th - np.array([item for sublist in th_g2 for item in sublist])
th_gn0 = th - np.array(th_g2)
   
cellpop_gn0 = cellpop_g2

g = 2

while cellpop_gn0 > 0:
    cellpop_gn = 0
    th_gn = []
    sp_gn = []

    for cell_n in range(cellpop_gn0):
        cell_name = 'g'+str(g)+'_c'+str(cell_n+1)
        xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th_gn0[cell_n],sp_gn0[cell_n],Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')
        
        tout_all = tout_all + (th-th_gn0[cell_n])*3600
        
        xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),20)))
        xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),20)))
        tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),20)))
        
        
        np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutS.txt')), xoutS_lite, delimiter='\t')
        np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutG.txt')), xoutG_lite, delimiter='\t')
        # tout_all = tout_all+th_gn[cell_n]*3600
        np.savetxt(os.path.join(output_dir,str(cell_name+'_tout.txt')), tout_lite, delimiter='\t')
        
        cb_peaks_n, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
        
        if len(cb_peaks_n > 0):
            
            th_g2_n = tout_all[cb_peaks_n[0]]/3600
            
            
            sp_g2_n = xoutS_all[cb_peaks_n[0]]
            
            
            th_gn.append(th-th_g2_n)
            th_gn.append(th-th_g2_n)
            # for p in range(len(th_g2_n)):
            #     th_gn.append(th_g2_n[p])
            
            
            sp_gn.append(sp_g2_n)
            sp_gn.append(sp_g2_n)
            
            cellpop_gn = cellpop_gn + 2
    
    if len(sp_gn) == 1:
        sp_gn0 = sp_gn[0]
        th_gn0 = np.array(th_gn).flatten()
    elif len(sp_gn) > 1:    
        # for i in range(1,len(sp_gn)):
        #      sp_gn0 = np.concatenate((sp_gn0,sp_gn[i]),axis=0)
        # sp_gn0 = np.array([item for sublist in sp_gn for item in sublist])
        sp_gn0 = sp_gn
        
        th_gn0 = np.array(th_gn).flatten()
    
    cellpop_gn0 = cellpop_gn
    if cellpop_gn0 > 0:
        g += 1
    # cellpop_total = cellpop_total + cellpop_gn
    
#%% observe output


def timecourse(species, gx_cx):

    x_s = np.loadtxt(os.path.join(wd,'output','cellpop_test',str(gx_cx+'_xoutS.txt')),delimiter='\t')
    tout_all = np.loadtxt(os.path.join(wd,'output','cellpop_test',str(gx_cx+'_tout.txt')),delimiter='\t')
    
    x_t = x_s[:, list(species_all).index(species)]
    plt.plot(tout_all/3600, x_t)
    plt.ylabel('species: '+str(species))
    plt.xlabel('time(h)')
    plt.ylim(0, max(x_t)*1.25)

    plt.show

#%%

# tout_g3c1 = np.loadtxt(os.path.join(wd,'output','cellpop_test','g3_c1_tout.txt'),delimiter='\t')
# xoutS_g3c1 = np.loadtxt(os.path.join(wd,'output','cellpop_test','g3_c1_xoutS.txt'),delimiter='\t')
# xoutS2 = np.loadtxt(os.path.join(wd,'output','cellpop_test','g1_c2_xoutS.txt'),delimiter='\t')
# xoutS3 = np.loadtxt(os.path.join(wd,'output','cellpop_test','g1_c3_xoutS.txt'),delimiter='\t')


timecourse('Mb','g2_c6')
# timecourse('Mb',xoutS2,tout)
# timecourse('Mb',xoutS3,tout)
    

#%% test - find plateaus


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
    
    return plateaus 

#%%
from scipy.signal import find_peaks

th = 72
xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD, th, species_initializations, Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')


plateaus = find_plateaus(xoutS_all[:, list(species_all).index('Mb')], min_length=120, tolerance = 0.95, smoothing=15)
cb_peaks, propps = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
# xoutS_all_g1c2 = np.loadtxt(os.path.join(wd,'output','cellpop_test','g1_c2_xoutS.txt'),delimiter='\t')
# tout_all_g1c2 = np.loadtxt(os.path.join(wd,'output','cellpop_test','g1_c2_tout.txt'),delimiter='\t')

# plateaus_g1c2 = find_plateaus(xoutS_all_g1c2[:, list(species_all).index('Mb')], min_length=120, tolerance = 0.95, smoothing=15)
#%%

plateausT = plateaus[:,0]*30/3600
plat_idxs = []
for idx,val in enumerate(cb_peaks):
    if (idx+1==len(cb_peaks)) and (plateausT[-1]<val):
        continue
    else:
        qq = np.nonzero((np.where((val*30.0/3600.0)<plateausT,1,0)))
        plat_idxs.append(qq[0][0])
        
#%%

for ii in range(len(plat_idxs)):
    tt = th-plateausT[plat_idxs[ii]]
    # th2run = np.vstack([th2run,tt])
    ss = xoutS_all[plateaus[plat_idxs[ii],0],:]
    # speciesInits = np.vstack([speciesInits,ss])
    # cell_lineage = np.vstack([cell_lineage, [counterC+1, 0, 0]])
    # cellsTot += 1