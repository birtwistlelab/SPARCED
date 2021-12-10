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
import time
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

#%%

def find_dp(xoutS,tout,species_all=species_all):
    data = xoutS[:,list(species_all).index('Mb')]
    p,_ = find_peaks(data,height=30)
    b = (np.diff(np.sign(np.diff(data))) > 0).nonzero()[0] + 1
    
    if sum(b>p[0]) > 0:
    
        dp = int(b[b>p[0]][0])
        
    else:
        dp = np.nan
    
    return(dp)


#%% stochastic run

ts = 30
flagD = 0
solver = model.getSolver()  # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints

# th = 36

# xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD, th, species_initializations, Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')

#
from scipy.signal import find_peaks
import itertools
#%%

cellpop_g1 = 5

th = 72

flagD = 0
th_g2 = []
sp_g2 = []
lin_g2 = []

cellpop_g2 = 0

for cell in range(cellpop_g1):
    cell_name = 'g1_c'+str(cell+1)
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD, th, species_initializations, Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')
    
    
    
    
    
    xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),2)))
    xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),2)))
    tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),2)))
    
    np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_xoutS.txt'),xoutS_lite,delimiter='\t')
    np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_xoutG.txt'),xoutG_lite,delimiter='\t')
    np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_tout.txt'),tout_lite,delimiter='\t')
    
    cb_peaks, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
    
    if len(cb_peaks)>0:
        

        dp = find_dp(xoutS_all,tout_all)
    
        if ~np.isnan(dp):
            th_g2_cell = tout_all[dp]/3600
            
            sp_g2_cell = xoutS_all[dp]
        
            th_g2.append(th_g2_cell)
            th_g2.append(th_g2_cell)
            sp_g2.append(sp_g2_cell)
            sp_g2.append(sp_g2_cell)
            lin_g2.append('c'+str(cell+1))
            lin_g2.append('c'+str(cell+1))
            
            xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,dp,2)))
            xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,dp,2)))
            tout_lite = np.array(list(itertools.islice(tout_all,0,dp,2)))
            
            np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_xoutS.txt'),xoutS_lite,delimiter='\t')
            np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_xoutG.txt'),xoutG_lite,delimiter='\t')
            np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_tout.txt'),tout_lite,delimiter='\t')
            
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

tout = np.loadtxt(os.path.join(wd,'output','cellpop_test_0','g1_c1_tout.txt'),delimiter='\t')
xoutS1 = np.loadtxt(os.path.join(wd,'output','cellpop_test_0','g1_c1_xoutS.txt'),delimiter='\t')
xoutS2 = np.loadtxt(os.path.join(wd,'output','cellpop_test_0','g1_c2_xoutS.txt'),delimiter='\t')
xoutS3 = np.loadtxt(os.path.join(wd,'output','cellpop_test_0','g1_c3_xoutS.txt'),delimiter='\t')


timecourse('Mb',xoutS1,tout)
timecourse('Mb',xoutS2,tout)
timecourse('Mb',xoutS3,tout)


#%%

output_dir = os.path.join(wd,'output','cellpop_test')

sp_gn0 = np.array(sp_g2)

lin_gn0 = lin_g2

th_gn0 = th - np.array(th_g2)
   
cellpop_gn0 = cellpop_g2

g = 2

while cellpop_gn0 > 0:
    cellpop_gn = 0
    th_gn = []
    sp_gn = []
    lin_gn = []

    for cell_n in range(cellpop_gn0):
        cell_name = 'g'+str(g)+'_c'+str(cell_n+1)+'_lin_'+str(lin_gn0[cell_n])
        xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th_gn0[cell_n],sp_gn0[cell_n],Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')
        
        tout_all = tout_all + (th-th_gn0[cell_n])*3600
        
        xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),20)))
        xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),20)))
        tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),20)))
        
        
        np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutS.txt')), xoutS_lite, delimiter='\t')
        np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutG.txt')), xoutG_lite, delimiter='\t')

        np.savetxt(os.path.join(output_dir,str(cell_name+'_tout.txt')), tout_lite, delimiter='\t')
        
        cb_peaks_n, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
        
        if len(cb_peaks_n > 0):
            

            dp = find_dp(xoutS_all,tout_all)
            
            
            if ~np.isnan(dp):
            
                th_g2_n = tout_all[dp]/3600
                
                sp_g2_n = xoutS_all[dp]
                
                lin_g2_n = str(lin_gn0[cell_n])+'c'+str(cell_n+1)
                
                th_gn.append(th-th_g2_n)
                th_gn.append(th-th_g2_n)

                
                
                sp_gn.append(sp_g2_n)
                sp_gn.append(sp_g2_n)
                
                lin_gn.append(lin_g2_n)
                lin_gn.append(lin_g2_n)
                
                xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,dp,20)))
                xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,dp,20)))
                tout_lite = np.array(list(itertools.islice(tout_all,0,dp,20)))
                
                np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_xoutS.txt'),xoutS_lite,delimiter='\t')
                np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_xoutG.txt'),xoutG_lite,delimiter='\t')
                np.savetxt(os.path.join(wd,'output','cellpop_test',cell_name+'_tout.txt'),tout_lite,delimiter='\t')
                    
                cellpop_gn = cellpop_gn + 2
    
    if len(sp_gn) == 1:
        sp_gn0 = sp_gn[0]
        th_gn0 = np.array(th_gn).flatten()
    elif len(sp_gn) > 1:    

        sp_gn0 = sp_gn
        
        th_gn0 = np.array(th_gn).flatten()
    
    lin_gn0 = lin_gn
    
    cellpop_gn0 = cellpop_gn
    if cellpop_gn0 > 0:
        g += 1

    
#%% observe output

output_dir = "cellpop_test_parallel"
def timecourse_gc(output_dir,species, gx_cx,tx='default'):
    
    output_ls = os.listdir(os.path.join(wd,'output',output_dir))
    
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx),output_ls))

    # x_s = np.loadtxt(os.path.join(wd,'output','cellpop_test',str(gx_cx+'_xoutS.txt')),delimiter='\t')
    # tout_all = np.loadtxt(os.path.join(wd,'output','cellpop_test',str(gx_cx+'_tout.txt')),delimiter='\t')
    
    x_s = np.loadtxt(os.path.join(wd,'output',output_dir,str(outputs_gc[2])),delimiter='\t')
    tout_all = np.loadtxt(os.path.join(wd,'output',output_dir,str(outputs_gc[0])),delimiter='\t')
    
    x_t = x_s[:, list(species_all).index(species)]
    plt.plot(tout_all/3600, x_t)
    plt.ylabel('species: '+str(species))
    plt.xlabel('time(h)')
    plt.ylim(0, max(x_t)*1.25)
    if type(tx)==str:
        plt.xlim(0, round(max(tout_all)/3600))
    elif type(tx)==int or type(tx)==float:
        plt.xlim(0,tx)
        
    plt.title(gx_cx)

    plt.show

#%%

timecourse_gc(output_dir,'Mb','g1_c1')


    
#%% cell lineage tracking

def find_dp_gc(x_s,species_all=species_all):
    
    # output_ls = os.listdir(os.path.join(wd,'output',output_dir))
    
    # outputs_gc = list(filter(lambda x:x.startswith(gx_cx),output_ls))
    
    # x_s = np.loadtxt(os.path.join(wd,'output',output_dir,str(outputs_gc[2])),delimiter='\t')
    
    
    
    data = x_s[:,list(species_all).index('Mb')]
    p,_ = find_peaks(data,height=30)
    # x = tout_g1c3/3600
    if len(p)!=0:
    
        b = (np.diff(np.sign(np.diff(data))) > 0).nonzero()[0] + 1
        
        if len(b)!=0:
            dp = int(b[b>p[0]][0])
        else:
            dp = np.nan
    else:
        dp = np.nan
    # th_dp = x[b[b>p[0]][0]]
    
    return(dp)


def timecourse_dp(outout_dir,species,gx_cx,species_all=species_all):
    output_ls = os.listdir(os.path.join(wd,'output',output_dir))
    
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),output_ls))
    
    xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
    tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    # xoutG_file = list(filter(lambda x:x.endswith('xoutG.txt'),outputs_gc))[0]
    
    
    x_s = np.loadtxt(os.path.join(wd,'output',output_dir,xoutS_file),delimiter='\t')
    th = int(max(np.loadtxt(os.path.join(wd,'output',output_dir,'g1_c1_tout.txt'),delimiter='\t')/3600))
    

    
    dp = find_dp_gc(x_s,species_all=species_all)

    # x_s = np.loadtxt(os.path.join(wd,'output','cellpop_test',str(gx_cx+'_xoutS.txt')),delimiter='\t')
    # tout_all = np.loadtxt(os.path.join(wd,'output','cellpop_test',str(gx_cx+'_tout.txt')),delimiter='\t')
    
    # x_s = np.loadtxt(os.path.join(wd,'output',output_dir,str(outputs_gc[2])),delimiter='\t')[:(dp+1)]
    tout_all = np.loadtxt(os.path.join(wd,'output',output_dir,tout_file),delimiter='\t')
    
    # th_dp0 = th - max(tout_all)/3600
    
    if ~np.isnan(dp):
        x_s = x_s[:(dp+1)]
        
        tout_all = tout_all[:(dp+1)]
    
    # tout_all = tout_all + th_dp0*3600
    
    x_t = x_s[:, list(species_all).index(species)]
    plt.plot(tout_all/3600, x_t)
    plt.ylabel('species: '+str(species))
    plt.xlabel('time(h)')
    plt.ylim(0, max(x_t)*1.25)
    plt.xlim(0, th)
    plt.title(gx_cx)

    plt.show
    



# cell = 'g1_c1'

# timecourse_gc('Mb',cell)

#%% clean


cellpop_g1 = 5

th = 72

flagD = 0

output_dir = os.path.join(wd,'output','cellpop_test')

th_g2 = []
sp_g2 = []
lin_g2 = []

cellpop_g2 = 0

for cell in range(cellpop_g1):
    cell_name = 'g1_c'+str(cell+1)
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD, th, species_initializations, Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')
    
    
    
    
    
    xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),20)))
    xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),20)))
    tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),20)))
    
    np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutS.txt')),xoutS_lite,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutG.txt')),xoutG_lite,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name+'_tout.txt')),tout_lite,delimiter='\t')
    
    cb_peaks, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
    
    if len(cb_peaks)>0:
        

        dp = find_dp(xoutS_all,tout_all)
    
        if ~np.isnan(dp):
            th_g2_cell = tout_all[dp]/3600
            
            sp_g2_cell = xoutS_all[dp]
        
            th_g2.append(th_g2_cell)
            th_g2.append(th_g2_cell)
            sp_g2.append(sp_g2_cell)
            sp_g2.append(sp_g2_cell)
            lin_g2.append('c'+str(cell+1))
            lin_g2.append('c'+str(cell+1))
            
            xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,dp,20)))
            xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,dp,20)))
            tout_lite = np.array(list(itertools.islice(tout_all,0,dp,20)))
            
            np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutS.txt')),xoutS_lite,delimiter='\t')
            np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutG.txt')),xoutG_lite,delimiter='\t')
            np.savetxt(os.path.join(output_dir,str(cell_name+'_tout.txt')),tout_lite,delimiter='\t')
            
            cellpop_g2 = cellpop_g2 + 2





sp_gn0 = np.array(sp_g2)

lin_gn0 = lin_g2

th_gn0 = th - np.array(th_g2)
   
cellpop_gn0 = cellpop_g2

g = 2

while cellpop_gn0 > 0:
    cellpop_gn = 0
    th_gn = []
    sp_gn = []
    lin_gn = []

    for cell_n in range(cellpop_gn0):
        cell_name = 'g'+str(g)+'_c'+str(cell_n+1)+'_lin_'+str(lin_gn0[cell_n])
        xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th_gn0[cell_n],sp_gn0[cell_n],Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')
        
        tout_all = tout_all + (th-th_gn0[cell_n])*3600
        
        xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),20)))
        xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),20)))
        tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),20)))
        
        
        np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutS.txt')), xoutS_lite, delimiter='\t')
        np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutG.txt')), xoutG_lite, delimiter='\t')

        np.savetxt(os.path.join(output_dir,str(cell_name+'_tout.txt')), tout_lite, delimiter='\t')
        
        cb_peaks_n, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
        
        if len(cb_peaks_n > 0):
            

            dp = find_dp(xoutS_all,tout_all)
            
            
            if ~np.isnan(dp):
            
                th_g2_n = tout_all[dp]/3600
                
                sp_g2_n = xoutS_all[dp]
                
                lin_g2_n = str(lin_gn0[cell_n])+'c'+str(cell_n+1)
                
                th_gn.append(th-th_g2_n)
                th_gn.append(th-th_g2_n)

                
                
                sp_gn.append(sp_g2_n)
                sp_gn.append(sp_g2_n)
                
                lin_gn.append(lin_g2_n)
                lin_gn.append(lin_g2_n)
                
                xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,dp,20)))
                xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,dp,20)))
                tout_lite = np.array(list(itertools.islice(tout_all,0,dp,20)))
                
                np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutS.txt')),xoutS_lite,delimiter='\t')
                np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutG.txt')),xoutG_lite,delimiter='\t')
                np.savetxt(os.path.join(output_dir,str(cell_name+'_tout.txt')),tout_lite,delimiter='\t')
                    
                cellpop_gn = cellpop_gn + 2
    
    if len(sp_gn) == 1:
        sp_gn0 = sp_gn[0]
        th_gn0 = np.array(th_gn).flatten()
    elif len(sp_gn) > 1:    

        sp_gn0 = sp_gn
        
        th_gn0 = np.array(th_gn).flatten()
    
    lin_gn0 = lin_gn
    
    cellpop_gn0 = cellpop_gn
    if cellpop_gn0 > 0:
        g += 1




#%% parallel processes
import multiprocessing

omics_input = 'OmicsData.txt'
genereg_input = 'GeneReg.txt'

cellpop_g1 = 5

th = 84

flagD = 0

output_dir = os.path.join(wd,'output','cellpop_test_parallel')
int_dir = os.path.join(output_dir,'intermediates')

th_g2 = []
sp_g2 = []
lin_g2 = []

cellpop_g2 = 0


def cell_g1(cell_n,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input):
    
    cell_name = 'g1_c'+str(cell_n+1)   
    
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    
    # xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),20)))
    # xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),20)))
    # tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),20)))
    
    # q_flagA.put(flagA)
    # gq.put(xoutG_lite)
    # tq.put(tout_lite)
    
    
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutS.txt'),xoutS_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutG.txt'),xoutG_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_tout.txt'),tout_all,delimiter='\t')

start_p = time.perf_counter()

processes = []

for c in range(cellpop_g1):
    p = multiprocessing.Process(target=cell_g1, args = [c,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()
    
#% test - threading
output_dir0 = os.path.join(wd,'output','cellpop_test_parallel')
def read_cell_g1(output_dir,g,cell_n):
    
    gx_cx = 'g'+str(g)+'_c'+str(int(cell_n+1))    
    
    outputs_ls = os.listdir(os.path.join(wd,'output',output_dir))
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),outputs_ls))
    xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
    tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    xoutG_file = list(filter(lambda x:x.endswith('xoutG.txt'),outputs_gc))[0]
    
    xoutS_all = np.loadtxt(os.path.join(output_dir,xoutS_file),delimiter='\t')
    xoutG_all = np.loadtxt(os.path.join(output_dir,xoutG_file),delimiter='\t')
    tout_all = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    
    

    
    cb_peaks, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
    
    results = {}
    
    if len(cb_peaks)>0:
        

        dp = find_dp(xoutS_all,tout_all)

    
        if ~np.isnan(dp):
            tdp_g2_cell = tout_all[dp]/3600
            
            sp_g2_cell = xoutS_all[dp]
            results['cell'] = int(cell_n+1)
            results['dp'] = dp
            results['th_g2'] = th- tdp_g2_cell    
            results['sp_g2'] = sp_g2_cell
            
            
    return results

# import threading
import concurrent.futures

results_all = []


with concurrent.futures.ThreadPoolExecutor() as executor:
    results_tpe = [executor.submit(read_cell_g1, output_dir0,1,cell_n) for cell_n in range(cellpop_g1)]
    
    for f in concurrent.futures.as_completed(results_tpe):
        results_all.append(f)

#%
results_f = [bool(results_all[i].result()) for i in range(len(results_all))]

results_notempty = np.array(results_all)[np.where(results_f)[0]]

results_actual = [r.result() for r in results_notempty]

th_g2 = [r['th_g2'] for r in results_actual]

th_g2 = th_g2 + th_g2

lin_g2 = [r['cell'] for r in results_actual]
lin_g2 = lin_g2 + lin_g2

sp_g2 = [r['sp_g2'] for r in results_actual]
sp_g2 = sp_g2 + sp_g2

#%

sp_gn0 = np.array(sp_g2)

lin_gn0 = lin_g2

th_gn0 = th_g2
   
cellpop_gn0 = len(th_g2)

g = 2

#%

def cell_gn(cell_n,lin_gn0,flagD,th,th_gc,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input):
    
    cell_name = 'g'+str(g)+'_c'+str(cell_n+1)+'_lin_'+str(lin_gn0[cell_n])  
    
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th_gc[cell_n],species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    
    tout_all = tout_all + (th-th_gc[cell_n])*3600

   
    
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutS.txt'),xoutS_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutG.txt'),xoutG_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_tout.txt'),tout_all,delimiter='\t')


def read_cell_gn(output_dir,g,cell_n):
    
    gx_cx = 'g'+str(g)+'_c'+str(int(cell_n+1))    
    
    outputs_ls = os.listdir(os.path.join(wd,'output',output_dir))
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),outputs_ls))
    xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
    tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    xoutG_file = list(filter(lambda x:x.endswith('xoutG.txt'),outputs_gc))[0]
    
    xoutS_all = np.loadtxt(os.path.join(output_dir,xoutS_file),delimiter='\t')
    xoutG_all = np.loadtxt(os.path.join(output_dir,xoutG_file),delimiter='\t')
    tout_all = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    
    

    
    cb_peaks, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
    
    results = {}
    
    if len(cb_peaks)>0:
        

        dp = find_dp(xoutS_all,tout_all)

    
        if ~np.isnan(dp):
            tdp_gn_cell = tout_all[dp]/3600
            
            sp_gn_cell = xoutS_all[dp]
            
            lin_gn_cell = str(lin_gn0[cell_n])+'c'+str(cell_n+1)
            
            results['cell'] = int(cell_n+1)
            results['dp'] = dp
            results['th_gn'] = th- tdp_gn_cell    
            results['sp_gn'] = sp_gn_cell
            results['lin'] = lin_gn_cell
            
            
    return results

#%
while cellpop_gn0 > 0:
    # cellpop_gn = 0
    # th_gn = []
    # sp_gn = []
    # lin_gn = []
    
    
        
    processes = []
    
    for c in range(cellpop_gn0):
        p = multiprocessing.Process(target=cell_gn, args = [c,lin_gn0,flagD,th,th_gn0,sp_gn0[c],Vn,Vc,model,wd,omics_input,genereg_input])
        p.start()
        processes.append(p)
        
    for process in processes:
        process.join()
    
    
    results_all = []

    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results_tpe = [executor.submit(read_cell_gn, output_dir0,g,cell_n) for cell_n in range(cellpop_gn0)]
        
        for f in concurrent.futures.as_completed(results_tpe):
            results_all.append(f)
            
    results_f = [bool(results_all[i].result()) for i in range(len(results_all))]

    results_notempty = np.array(results_all)[np.where(results_f)[0]]
    
    results_actual = [r.result() for r in results_notempty]
    
    th_gn = [r['th_gn'] for r in results_actual]
    
    th_gn = th_gn + th_gn
    
    lin_gn = [r['lin'] for r in results_actual]
    lin_gn = lin_gn + lin_gn
    

    
    
    sp_gn = [r['sp_gn'] for r in results_actual]
    sp_gn = sp_gn + sp_gn
    
    cellpop_gn = len(th_gn)
    
    cellpop_gn0 = cellpop_gn
    
    if cellpop_gn0 > 0:
        g += 1
        
        lin_gn0 = lin_gn
        sp_gn0 = sp_gn
        th_gn0 = th_gn


finish_p = time.perf_counter()

#%%

time_min = (finish_p - start_p)/60




#%%
# cell_n = 4
# x1 = np.loadtxt(os.path.join(output_dir,'c'+str(cell_n)+'_xoutS.txt'),delimiter='\t')
# # xoutG_all = np.loadtxt(os.path.join(output_dir,'c'+str(cell_n)+'_xoutG.txt'),delimiter='\t')
# t1 = np.loadtxt(os.path.join(output_dir,'c'+str(cell_n)+'_tout.txt'),delimiter='\t')

# timecourse('cPARP',x1,t1)


#%% test multi-threaded reading on previous outputs:
output_dir0 = os.path.join(wd,'output','cellpop_test_parallel')

# output_dir = os.path.join(wd,'output','cellpop_test_0')

def read_cell_gn(output_dir,g,cell_n):
    
    gx_cx = 'g'+str(g)+'_c'+str(int(cell_n+1))
    
    outputs_ls = os.listdir(output_dir)
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx),outputs_ls))
    xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
    tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    xoutG_file = list(filter(lambda x:x.endswith('xoutG.txt'),outputs_gc))[0]
    
    xoutS_all = np.loadtxt(os.path.join(output_dir,xoutS_file),delimiter='\t')
    xoutG_all = np.loadtxt(os.path.join(output_dir,xoutG_file),delimiter='\t')
    tout_all = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
        
    
    cb_peaks, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
    
    
    
    # results = []
    
    # results.append
    
    results = {}
    
    if len(cb_peaks)>0:
        

        dp = find_dp(xoutS_all,tout_all)
        
        # results.append(dp)
    
        if ~np.isnan(dp):
            tdp_g2_cell = tout_all[dp]/3600
            
            sp_g2_cell = xoutS_all[dp]
            
            # lin_g2_n = str(lin_gn0[cell_n])+'c'+str(cell_n+1)
    
            
            results['cell'] = int(cell_n+1)
            results['dp'] = dp
            results['th_g2'] = th - tdp_g2_cell    
            results['sp_g2'] = sp_g2_cell
            
            
    return results

results_all = []

g = 1


cellpop_g1 = 8

with concurrent.futures.ThreadPoolExecutor() as executor:
    results_tpe = [executor.submit(read_cell_gn, output_dir0, g, cell_n) for cell_n in range(cellpop_g1)]
    
    for f in concurrent.futures.as_completed(results_tpe):
        results_all.append(f)

#%%
results_f = [bool(results_all[i].result()) for i in range(len(results_all))]

results_notempty = np.array(results_all)[np.where(results_f)[0]]

results_actual = [r.result() for r in results_notempty]

th_g2 = [r['th_g2'] for r in results_actual]
    
#%% temp - not ready



for cell in range(cellpop_g1):
    cell_name = 'g1_c'+str(cell+1)
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD, th, species_initializations, Vn, Vc, model,wd,omics_input='OmicsData.txt',genereg_input='GeneReg.txt')
    
    
    
    
    
    xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),20)))
    xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),20)))
    tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),20)))
    
    np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutS.txt')),xoutS_lite,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutG.txt')),xoutG_lite,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name+'_tout.txt')),tout_lite,delimiter='\t')
    
    cb_peaks, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
    
    if len(cb_peaks)>0:
        

        dp = find_dp(xoutS_all,tout_all)
    
        if ~np.isnan(dp):
            th_g2_cell = tout_all[dp]/3600
            
            sp_g2_cell = xoutS_all[dp]
        
            th_g2.append(th_g2_cell)
            th_g2.append(th_g2_cell)
            sp_g2.append(sp_g2_cell)
            sp_g2.append(sp_g2_cell)
            lin_g2.append('c'+str(cell+1))
            lin_g2.append('c'+str(cell+1))
            
            xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,dp,20)))
            xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,dp,20)))
            tout_lite = np.array(list(itertools.islice(tout_all,0,dp,20)))
            
            np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutS.txt')),xoutS_lite,delimiter='\t')
            np.savetxt(os.path.join(output_dir,str(cell_name+'_xoutG.txt')),xoutG_lite,delimiter='\t')
            np.savetxt(os.path.join(output_dir,str(cell_name+'_tout.txt')),tout_lite,delimiter='\t')
            
            cellpop_g2 = cellpop_g2 + 2





sp_gn0 = np.array(sp_g2)

lin_gn0 = lin_g2

th_gn0 = th - np.array(th_g2)
   
cellpop_gn0 = cellpop_g2