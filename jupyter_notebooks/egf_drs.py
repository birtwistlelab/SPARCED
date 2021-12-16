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
from scipy.signal import find_peaks
import itertools

parser = argparse.ArgumentParser(description='Input doses in uM')
parser.add_argument('--dose', metavar='dose', help='input dose in nM', default = 100.0)
parser.add_argument('--cellpop', metavar='cellpop', help='starting cell population', default = 5)
parser.add_argument('--td',metavar='td', help='cell line doubling time (hrs) ', default = 48)
args = parser.parse_args()

wd = str(os.getcwd()).replace("jupyter_notebooks","")


sim_name = 'egf_dose_response'

output_path = os.path.join(wd,'output',sim_name)

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


dose = float(args.dose)

sp_input = pd.read_csv(os.path.join(wd,'initializer','species_au565.txt'),sep='\t',header=None,index_col=0,squeeze=True)
# sp_input['lapatinib'] = dose


species_initializations = np.array(sp_input)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
species_initializations[155:162] = STIMligs


output_dose = os.path.join(output_path,'egf_'+str(float(args.dose)))

if not os.path.exists(output_dose):
    os.mkdir(output_dose)

# np.linspace(0, 30) # set timepoint


import time
import multiprocessing

th = 48



def pre_incubate(cell_n,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input):
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)
    
    np.savetxt(os.path.join(output_dose,'c'+str(cell_n+1)+'_preinc.txt'),xoutS_all[-1],delimiter='\t')
    
    
# start = time.perf_counter()


if __name__ == "__main__":


    
    processes = []
    
    for c in range(cell_pop):
        p = multiprocessing.Process(target=pre_incubate, args = [c,flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input])
        p.start()
        processes.append(p)
        
    for process in processes:
        process.join()
    


# finish = time.perf_counter()

# print(f'Finished in {round(finish-start,2)} seconds(s)')




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

#%% drs no queue - 72 hr

import itertools

STIMligs = [0.0,0,0,0,0,0,0.0]

doubling_time = float(args.td)

STIMligs[0] = dose

th = doubling_time + 10.0


cellpop_g1 = cell_pop

output_dir = output_dose

def cell_g1(cell_n,flagD,th,Vn,Vc,model,wd,omics_input,genereg_input):
    
    cell_name = 'g1_c'+str(cell_n+1)
    
    s_preinc_i = np.loadtxt(os.path.join(output_dose,'c'+str(cell_n+1)+'_preinc.txt'),delimiter='\t')
    sp_input = s_preinc_i
    species_initializations = np.array(sp_input)
    species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
    species_initializations[155:162] = STIMligs
    
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th,species_initializations,Vn,Vc,model,wd,omics_input,genereg_input)

    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutS.txt'),xoutS_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutG.txt'),xoutG_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_tout.txt'),tout_all,delimiter='\t')


start_p = time.perf_counter()

processes = []

for c in range(cellpop_g1):
    p = multiprocessing.Process(target=cell_g1, args = [c,flagD,th,Vn,Vc,model,wd,omics_input,genereg_input])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()

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
    
    xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),20)))
    xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),20)))
    tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),20)))
    
    if len(cb_peaks)>0:
        

        dp = find_dp(xoutS_all,tout_all)

    
        if ~np.isnan(dp):
            
            parp_dp = float(xoutS_all[dp,list(species_all).index('PARP')])
            cparp_dp = float(xoutS_all[dp,list(species_all).index('cPARP')])
            
            if parp_dp > cparp_dp:
            
                
                tdp_g2_cell = tout_all[dp]/3600
                
                sp_g2_cell = xoutS_all[dp]
                
                lin_g2_cell = 'c'+str(int(cell_n+1))
                
                results['cell'] = int(cell_n+1)
                results['dp'] = dp
                results['th_g2'] = th- tdp_g2_cell    
                results['sp_g2'] = sp_g2_cell
                results['lin'] = lin_g2_cell
                
                xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(dp+1),20)))
                xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(dp+1),20)))
                tout_lite = np.array(list(itertools.islice(tout_all,0,(dp+1),20)))


    
    np.savetxt(os.path.join(output_dir,xoutS_file),xoutS_lite,delimiter='\t')
    np.savetxt(os.path.join(output_dir,xoutG_file),xoutG_lite,delimiter='\t')
    np.savetxt(os.path.join(output_dir,tout_file),tout_lite,delimiter='\t')
        
            
            
    return results


import concurrent.futures

results_all = []


with concurrent.futures.ThreadPoolExecutor() as executor:
    results_tpe = [executor.submit(read_cell_g1, output_dir,1,cell_n) for cell_n in range(cellpop_g1)]
    
    for f in concurrent.futures.as_completed(results_tpe):
        results_all.append(f)

#%
results_f = [bool(results_all[i].result()) for i in range(len(results_all))]

results_notempty = np.array(results_all)[np.where(results_f)[0]]

results_actual = [r.result() for r in results_notempty]

th_g2 = [r['th_g2'] for r in results_actual]

th_g2 = th_g2 + th_g2

lin_g2 = [r['lin'] for r in results_actual]
lin_g2 = lin_g2 + lin_g2

sp_g2 = [r['sp_g2'] for r in results_actual]
sp_g2 = sp_g2 + sp_g2

#%

sp_gn0 = np.array(sp_g2)

lin_gn0 = lin_g2

th_gn0 = th_g2
   
cellpop_gn0 = len(th_g2)

g = 2


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
    
    xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(len(xoutS_all)-1),20)))
    xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(len(xoutG_all)-1),20)))
    tout_lite = np.array(list(itertools.islice(tout_all,0,(len(tout_all)-1),20)))   

    
    cb_peaks, _ = find_peaks(xoutS_all[:, list(species_all).index('Mb')],height=30)
    
    results = {}
    
    if len(cb_peaks)>0:
        

        dp = find_dp(xoutS_all,tout_all)

    
        if ~np.isnan(dp):
            
            parp_dp = float(xoutS_all[dp,list(species_all).index('PARP')])
            cparp_dp = float(xoutS_all[dp,list(species_all).index('cPARP')])
            
            if parp_dp > cparp_dp:
                
            
                tdp_gn_cell = tout_all[dp]/3600
                
                sp_gn_cell = xoutS_all[dp]
                
                lin_gn_cell = str(lin_gn0[cell_n])+'c'+str(cell_n+1)
                
                results['cell'] = int(cell_n+1)
                results['dp'] = dp
                results['th_gn'] = th- tdp_gn_cell    
                results['sp_gn'] = sp_gn_cell
                results['lin'] = lin_gn_cell

            xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(dp+1),20)))
            xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(dp+1),20)))
            tout_lite = np.array(list(itertools.islice(tout_all,0,(dp+1),20)))   

                

    
    np.savetxt(os.path.join(output_dir,xoutS_file),xoutS_lite,delimiter='\t')
    np.savetxt(os.path.join(output_dir,xoutG_file),xoutG_lite,delimiter='\t')
    np.savetxt(os.path.join(output_dir,tout_file),tout_lite,delimiter='\t')
            
            
    return results



#%
while cellpop_gn0 > 0:
 
       
    processes = []
    
    for c in range(cellpop_gn0):
        p = multiprocessing.Process(target=cell_gn, args = [c,lin_gn0,flagD,th,th_gn0,sp_gn0[c],Vn,Vc,model,wd,omics_input,genereg_input])
        p.start()
        processes.append(p)
        
    for process in processes:
        process.join()
    
    
    results_all = []

    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results_tpe = [executor.submit(read_cell_gn, output_dir,g,cell_n) for cell_n in range(cellpop_gn0)]
        
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

#%% debug
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300

def timecourse_gc(output_dir,species, gx_cx,tx='default'):
    
    output_ls = os.listdir(os.path.join(wd,output_dir))
    
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

def find_dp_gc(x_s,species_all=species_all):
    

    data = x_s[:,list(species_all).index('Mb')]
    p,_ = find_peaks(data,height=30)

    if len(p)!=0:
    
        b = (np.diff(np.sign(np.diff(data))) > 0).nonzero()[0] + 1
        
        if len(b)!=0:
            dp = int(b[b>p[0]][0])
        else:
            dp = np.nan
    else:
        dp = np.nan

    
    return(dp)

    
def timecourse_dp(outout_dir,species,gx_cx,species_all=species_all):
    output_ls = os.listdir(os.path.join(wd,'output',output_dir))
    
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),output_ls))
    
    xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
    tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    
    
    x_s = np.loadtxt(os.path.join(wd,'output',output_dir,xoutS_file),delimiter='\t')
    

    
    dp = find_dp_gc(x_s,species_all=species_all)


    tout_all = np.loadtxt(os.path.join(wd,'output',output_dir,tout_file),delimiter='\t')

    
    if ~np.isnan(dp):
        x_s = x_s[:(dp+1)]
        
        tout_all = tout_all[:(dp+1)]
    

    
    x_t = x_s[:, list(species_all).index(species)]
    plt.plot(tout_all/3600, x_t)
    plt.ylabel('species: '+str(species))
    plt.xlabel('time(h)')
    plt.ylim(0, max(x_t)*1.25)
    plt.xlim(0, th)
    plt.title(gx_cx)

    plt.show
    
#%% find number of cells at doubling time



