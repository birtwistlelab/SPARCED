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

mpl.rcParams['figure.dpi'] = 300

parser = argparse.ArgumentParser(description='Input doses in uM')
parser.add_argument('--dose', metavar='dose', help='input dose in nM', default = 100.0)
parser.add_argument('--cellpop', metavar='cellpop', help='starting cell population', default = 5)
parser.add_argument('--td',metavar='td', help='cell line doubling time (hrs) ', default = 48)
parser.add_argument('--sim_name',metavar='sim_name', help='insert exp name', default = 'egf_dose_response')
args = parser.parse_args()

wd = str(os.getcwd()).replace("jupyter_notebooks","")


sim_name = str(args.sim_name)

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




th = 48

output_dir = output_dose

#%%

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
#%%
# results_tmax = []
# with concurrent.futures.ThreadPoolExecutor() as executor:
#     results_tout = [executor.submit(read_tout, output_dir, tout_file) for tout_file in tout_files]

# for f in concurrent.futures.as_completed(results_tout):
#     results_tmax.append(f.result())
output_palmetto = os.path.join('/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/egf_drs_test/egf_100.0')

output_test1 = os.path.join('output/egf_dose_response_0/egf_100.0')

#%% debug


def read_tout(output_dir,tout_file):
    tout = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    
    if np.shape(tout):
        tmax = max(tout)
    else:
        tmax = float(tout)

    return (tmax)



def get_tmax(output_dir,read_tout=read_tout):
    
    outputs_ls = os.listdir(output_dir)
    tout_files = list(filter(lambda x:x.endswith('tout.txt'),outputs_ls))
                      
    results_tmax = []                 
                  
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results_tout = [executor.submit(read_tout, output_dir, tout_file) for tout_file in tout_files]

    for f in concurrent.futures.as_completed(results_tout):
        results_tmax.append(f.result())
        
    results_tmax = np.array(results_tmax)
    
    tmax = max(results_tmax)
    
    return(tmax)



#%%
output_dir = output_palmetto
timecourse_tmax = get_tmax(output_dir)  
# timecourse_tmax = get_tmax(output_palmetto)


#%%
def timecourse_gc(output_dir,species, gx_cx,get_tmax=get_tmax,read_tout=read_tout,tx='default'):
    
    outputs_ls = os.listdir(os.path.join(wd,output_dir))
    
    outputs_gc = list(filter(lambda x:x.startswith(str(gx_cx+'_')),outputs_ls))
    
    xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
    tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    
    if 'timecourse_tmax' not in globals():
        tmax = get_tmax(output_dir)
    else:
        tmax = timecourse_tmax



    # x_s = np.loadtxt(os.path.join(wd,'output','cellpop_test',str(gx_cx+'_xoutS.txt')),delimiter='\t')
    # tout_all = np.loadtxt(os.path.join(wd,'output','cellpop_test',str(gx_cx+'_tout.txt')),delimiter='\t')
    
    x_s = np.loadtxt(os.path.join(output_dir,xoutS_file),delimiter='\t')
    tout_all = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    
    x_t = x_s[:, list(species_all).index(species)]
    plt.plot(tout_all/3600, x_t)
    plt.ylabel('species: '+str(species))
    plt.xlabel('time(h)')
    plt.ylim(0, max(x_t)*1.25)
    if type(tx)==str:
        # plt.xlim(0, round(max(tout_all)/3600))
        plt.xlim(0, round(tmax/3600))
        
    elif type(tx)==int or type(tx)==float:
        plt.xlim(0,tx)
        
    plt.title(gx_cx)

    plt.show

#%%

timecourse_gc(output_palmetto,'Mb','g1_c1')

#%% test dataset
outputs_ls_test = os.listdir(os.path.join(wd,output_test1))
outputs_gc_test = list(filter(lambda x:x.startswith(str('g1_c3_')),outputs_ls_test))
xoutS_file_test = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc_test))[0]
tout_file_test = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc_test))[0]

x_s_test = np.loadtxt(os.path.join(output_test1,xoutS_file_test),delimiter='\t')
tout_test = np.loadtxt(os.path.join(output_test1,tout_file_test),delimiter='\t')

data_test = x_s_test[:,list(species_all).index('Mb')]
p,_ = find_peaks(data_test,height=30)


if len(p)!=0:
    
        b = (np.diff(np.sign(np.diff(data_test))) > 0).nonzero()[0] + 1
        
        if len(b)!=0:
            dp = int(b[b>p[0]][0])
        else:
            dp = np.nan
else:
        dp = np.nan

x_s_test2 = x_s_test[181:,:]
tout_test2 = tout_test [181:]

data_test2 = x_s_test2[:,list(species_all).index('Mb')]
p2,_ = find_peaks(data_test,height=30)


if len(p2)!=0:
    
        b2 = (np.diff(np.sign(np.diff(data_test2))) > 0).nonzero()[0] + 1
        
        if len(b2)!=0:
            dp2 = int(b2[b2>p2[0]][0])
        else:
            dp2 = np.nan
else:
        dp2 = np.nan
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
    
#%% find number of alive cells at given time



output_ls = os.listdir(output_dir)

xoutS_files = list(filter(lambda x:x.endswith('xoutS.txt'),output_ls))

tout_files = list(filter(lambda x:x.endswith('tout.txt'),output_ls))

tout_test = np.loadtxt(os.path.join(wd,'output',output_dir,tout_files[0]),delimiter='\t')

ti = 46.31

cells_all = np.zeros(len(tout_files))

for c in range(len(cells_all)):
    tout_c = np.loadtxt(os.path.join(output_dir,tout_files[c]),delimiter='\t')/3600
    tc_max = max(tout_c)
    tc_min = min(tout_c)
    if tc_max > ti and tc_min < ti:
        idx_ti = np.abs(tout_c - ti).argmin()
        
        gx_cx = str(tout_files[c].split('_')[:2][0]) + '_' +  str(tout_files[c].split('_')[:2][1])
        
        xoutS_file = list(filter(lambda x:x.startswith(gx_cx),xoutS_files))[0]
        
        xout_c = np.loadtxt(os.path.join(output_dir,xoutS_file),delimiter='\t')[idx_ti,:]
        
        cPARP_ti = xout_c[species_all.index('cPARP')]
        PARP_ti = xout_c[species_all.index('PARP')]
        
        if PARP_ti > cPARP_ti:
        
            cells_all[c] = 1.0

cells_alive_ti = int(sum(cells_all))

#%%

output_test3 = os.path.join('/media','arnab','bigchungus','projects','ccle_egf_drs','SPARCED','output','egf_drs_r1','egf_0.1')

output_ls3 = os.listdir(output_test3)

xoutS_files3 = list(filter(lambda x:x.endswith('xoutS.txt'),output_ls3))
tout_files3 = list(filter(lambda x:x.endswith('tout.txt'),output_ls3))

cells_all = np.zeros(len(tout_files3))

ti = 46.31

for c in range(len(cells_all)):
    tout_c = np.loadtxt(os.path.join(output_test3,tout_files3[c]),delimiter='\t')/3600
    
    if np.shape(tout_c):
    
        tc_max = max(tout_c)
        tc_min = min(tout_c)
    else:
        tc_max = tout_c
        tc_min = tout_c
       
        
    if tc_max > ti and tc_min < ti:
        idx_ti = np.abs(tout_c - ti).argmin()
        
        gx_cx = str(tout_files3[c].split('_')[:2][0]) + '_' +  str(tout_files3[c].split('_')[:2][1])
        
        xoutS_file = list(filter(lambda x:x.startswith(gx_cx),xoutS_files3))[0]
        
        xout_c = np.loadtxt(os.path.join(output_test3,xoutS_file),delimiter='\t')[idx_ti,:]
        
        cPARP_ti = xout_c[species_all.index('cPARP')]
        PARP_ti = xout_c[species_all.index('PARP')]
        
        if PARP_ti > cPARP_ti:
        
            cells_all[c] = 1.0

cells_alive_ti = int(sum(cells_all))

#%%

timecourse_gc(output_test3,'Mb','g2_c15')


def cellcount_dir(output_dir,ti):
    output_ls = os.listdir(output_dir)
    xoutS_files = list(filter(lambda x:x.endswith('xoutS.txt'),output_ls))
    tout_files = list(filter(lambda x:x.endswith('tout.txt'),output_ls))
    cells_all = np.zeros(len(tout_files))
    
    for c in range(len(cells_all)):
        tout_c = np.loadtxt(os.path.join(output_dir,tout_files[c]),delimiter='\t')/3600
        
        if np.shape(tout_c):
        
            tc_max = max(tout_c)
            tc_min = min(tout_c)
        else:
            tc_max = tout_c
            tc_min = tout_c
           
            
        if tc_max > ti and tc_min < ti:
            idx_ti = np.abs(tout_c - ti).argmin()
            
            gx_cx = str(tout_files[c].split('_')[:2][0]) + '_' +  str(tout_files[c].split('_')[:2][1])
            
            xoutS_file = list(filter(lambda x:x.startswith(str(gx_cx+'_')),xoutS_files))[0]
            
            xout_c = np.loadtxt(os.path.join(output_dir,xoutS_file),delimiter='\t')[idx_ti,:]
            
            cPARP_ti = xout_c[species_all.index('cPARP')]
            PARP_ti = xout_c[species_all.index('PARP')]
            
            if PARP_ti > cPARP_ti:
            
                cells_all[c] = 1.0
        
        cells_alive_ti = int(sum(cells_all))
        
    return cells_alive_ti 

#%%

output_palmetto_main = os.path.join('/media','arnab','bigchungus','projects','ccle_egf_drs','SPARCED','output','egf_drs')

dir_doses = os.listdir(os.path.join(output_palmetto_main,os.listdir(output_palmetto_main)[0]))

egf_doses = [float(str(x).split('_')[1]) for x in dir_doses]

ti = 46.31

results = pd.DataFrame()

cell_count_doses = np.zeros(len(dir_doses))

for d in range(len(dir_doses)):
    
    cell_count = cellcount_dir(os.path.join(output_palmetto_main,'egf_drs_r1',dir_doses[d]),ti)
    cell_count_doses[d] = cell_count
    
results['doses'] = egf_doses

results['r1'] = cell_count_doses

#%%

ti = 46.31

output_palmetto_main = os.path.join('/media','arnab','bigchungus','projects','ccle_egf_drs','SPARCED','output','egf_drs')

dir_doses = os.listdir(os.path.join(output_palmetto_main,os.listdir(output_palmetto_main)[0]))

egf_doses = [float(str(x).split('_')[1]) for x in dir_doses]

results = pd.DataFrame()

for r in range(3):
    
    cell_count_doses = np.zeros(len(dir_doses))
    
    for d in range(len(dir_doses)):
        
        cell_count = cellcount_dir(os.path.join(output_palmetto_main,'egf_drs_r'+str(r+1),dir_doses[d]),ti)
        cell_count_doses[d] = cell_count
        
    results['r'+str(r+1)] = cell_count_doses

results['doses'] = egf_doses

results = results.set_index('doses')

results = results.sort_index()

results['total'] = results.sum(1).values

#%%

ti = 46.31

output_palmetto_main = os.path.join('/media','arnab','bigchungus','projects','ccle_egf_drs','SPARCED','output','egf_drs')

dir_doses = os.listdir(os.path.join(output_palmetto_main,os.listdir(output_palmetto_main)[0]))

egf_doses = [float(str(x).split('_')[1]) for x in dir_doses]

results1 = pd.DataFrame()

for r in range(4):
    
    cell_count_doses = np.zeros(len(dir_doses))
    
    for d in range(len(dir_doses)):
        
        cell_count = cellcount_dir(os.path.join(output_palmetto_main,'egf_drs_r'+str(r+1),dir_doses[d]),ti)
        cell_count_doses[d] = cell_count
        
    results1['r'+str(r+1)] = cell_count_doses

results1['doses'] = egf_doses

results1 = results1.set_index('doses')

results1 = results1.sort_index()

results1['total'] = results1.sum(1).values


#%% set 2 - egf doses 0.0010-0.01 nM

ti = 46.31

output_palmetto_main = os.path.join('/media','arnab','bigchungus','projects','ccle_egf_drs','SPARCED','output','egf_drs')

dir_doses = os.listdir(os.path.join(output_palmetto_main,os.listdir(output_palmetto_main)[4]))

egf_doses = [float(str(x).split('_')[1]) for x in dir_doses]

results2 = pd.DataFrame()

for r in range(4,8):
    
    cell_count_doses = np.zeros(len(dir_doses))
    
    for d in range(len(dir_doses)):
        
        cell_count = cellcount_dir(os.path.join(output_palmetto_main,'egf_drs_r'+str(r+1),dir_doses[d]),ti)
        cell_count_doses[d] = cell_count
        
    results2['r'+str(r+1)] = cell_count_doses

results2['doses'] = egf_doses

results2 = results2.set_index('doses')

results2 = results2.sort_index()

results2['total'] = results2.sum(1).values

#%%