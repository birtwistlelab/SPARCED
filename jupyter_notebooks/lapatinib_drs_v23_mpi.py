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
parser.add_argument('--cellpop', metavar='cellpop', help='starting cell population', default = 3)
parser.add_argument('--td',metavar='td', help='cell line doubling time (hrs) ', default = 48)
parser.add_argument('--sim_name',metavar='sim_name', help='insert exp name', default = 'lapatinib_drs_test_mpi')
parser.add_argument('--mb_tr',metavar='mb_tr',help='Mb trough upper limit (nM)', default = 2.0)
parser.add_argument('--exp_time',metavar='exp_time',help='Enter experiment time in hours', default = 72.0)
parser.add_argument('--diag',metavar='diag',help='Toggle diagnostic mode', default = 'off')
parser.add_argument('--perturb',metavar='perturb',help='Specify perturbed species', default = 'lapatinib')
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


diag = str(args.diag)

#%%

species_all = list(model.getStateIds())

solver = model.getSolver()          # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

cell_pop = int(args.cellpop)


#%% preincubate

STIMligs = [0.0,0,0,0,0,0,0.0]

perturb = str(args.perturb)


dose = float(args.dose)*10e2

sp_input = pd.read_csv(os.path.join(wd,'initializer','species_au565.txt'),sep='\t',header=None,index_col=0,squeeze=True)
# sp_input['lapatinib'] = dose


species_initializations = np.array(sp_input)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
species_initializations[155:162] = STIMligs


dose_egf = float(args.egf)

if perturb == 'egf':
    output_dose = os.path.join(output_path,str(perturb)+'_'+str(float(dose_egf)))
    
else:
    output_dose = os.path.join(output_path,str(perturb)+'_'+str(float(args.dose)))

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
    
    print("Running gen0 cell(%d) on rank %d" %(cell_n+1, MY_RANK))
    
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
mb_tr = float(args.mb_tr)

def find_dp(xoutS,tout,species_all=species_all):
    data = xoutS[:,list(species_all).index('Mb')]
    p,_ = find_peaks(data,height=30)
    b = (np.diff(np.sign(np.diff(data))) > 0).nonzero()[0] + 1
    
    if len(b)!=0:
        b = np.array(b)[data[b]<mb_tr]
    
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
    
    if len(dp_all)!=0:
        dp_all_actual = list(np.array(dp_all)[data[dp_all]<mb_tr])
        dp_all = dp_all_actual

    
    return(dp_all)

#%% initiate gen1

exp_time = float(args.exp_time)

th = exp_time + 3.0


cellpop_g1 = cell_pop

output_dir = output_dose

def cell_g1(cell_n):
    
    print("Running gen1 cell(%d) on rank %d" %(cell_n+1, MY_RANK))
    
    cell_name = 'g1_c'+str(cell_n+1)
    
    # s_preinc_i = np.loadtxt(os.path.join(output_dose,'c'+str(cell_n+1)+'_preinc.txt'),delimiter='\t')
    x_s_g0 = np.loadtxt(os.path.join(output_dir,'g0_c'+str(cell_n+1)+'_xoutS.txt'),delimiter='\t')
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

    
with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
    if executor is not None:
        gen1_begins = time.time()
        print("Submitting gen1 stage from rank", MY_RANK)
        res = executor.map(cell_g1, range(cellpop_g1), unordered=True)
        num_done = len(list(res))
        gen1_ends = time.time()
        gen1_runtime = float(gen1_ends - gen1_begins)
        print("Finished", num_done,"gen1 tasks in",gen1_runtime,"seconds.")
        
#%% analyze gen1


def read_cell_g1(cell_n):
    
    gx_cx = 'g1_c'+str(int(cell_n+1))    
    
    outputs_ls = os.listdir(output_dir)
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),outputs_ls))
    xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
    tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    xoutG_file = list(filter(lambda x:x.endswith('xoutG.txt'),outputs_gc))[0]
    
    xoutS_g1 = np.loadtxt(os.path.join(output_dir,xoutS_file),delimiter='\t')
    xoutG_g1 = np.loadtxt(os.path.join(output_dir,xoutG_file),delimiter='\t')
    tout_g1 = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    
    xoutS_g0 = np.loadtxt(os.path.join(output_dir,'g0_c'+str(int(cell_n+1))+'_xoutS.txt'),delimiter='\t')
    tout_g0 = np.loadtxt(os.path.join(output_dir,'g0_c'+str(int(cell_n+1))+'_tout.txt'),delimiter='\t')
    tp_g0 = int(np.loadtxt(os.path.join(output_dir,'g0_c'+str(int(cell_n+1))+'_tp.txt'),delimiter='\t'))
    
    
    tneg_g0_min = max(tout_g0[:tp_g0]) - 16*3600
    
    tneg_idx_start = np.where(tout_g0[:tp_g0]>tneg_g0_min)[0][0]
    
    tout_g0_neg = tout_g0[:tp_g0][tneg_idx_start:tp_g0] - tout_g0[tp_g0]
    
    xoutS_new = np.concatenate((xoutS_g0[tneg_idx_start:tp_g0],xoutS_g1),axis=0)
    
    tout_new = np.concatenate((tout_g0_neg,tout_g1),axis=0)
   
    cb_peaks, _ = find_peaks(xoutS_new[:, list(species_all).index('Mb')],height=30)
    
    results = {}
    
    if diag == 'off':
    
        xoutS_lite = np.array(list(itertools.islice(xoutS_g1,0,(len(xoutS_g1)-1),20)))
        xoutG_lite = np.array(list(itertools.islice(xoutG_g1,0,(len(xoutG_g1)-1),20)))
        tout_lite = np.array(list(itertools.islice(tout_g1,0,(len(tout_g1)-1),20)))
    
    if len(cb_peaks)>0:
        
        dp_all = find_dp_all(xoutS_new)

        dp = np.nan

        if len(dp_all)>0:
            dp_idx = np.where(tout_new[dp_all]>0)[0][0]

            dp = dp_all[dp_idx]
  
        if ~np.isnan(dp):
            
            parp_dp = float(xoutS_new[dp,list(species_all).index('PARP')])
            cparp_dp = float(xoutS_new[dp,list(species_all).index('cPARP')])
            
            if parp_dp > cparp_dp:
            
                
                tdp_g2_cell = tout_new[dp]/3600
                
                sp_g2_cell = xoutS_new[dp]
                
                lin_g2_cell = 'c'+str(int(cell_n+1))
                
                results['cell'] = int(cell_n+1)
                results['dp'] = dp
                results['th_g2'] = th- tdp_g2_cell    
                results['lin'] = lin_g2_cell
                
                np.savetxt(os.path.join(output_dir,gx_cx+'_ic.txt'),sp_g2_cell,delimiter='\t')
                
                dp1 = np.where(tout_g1 == tout_new[dp])[0][0]
                
                if diag == 'off':
                
                    xoutS_lite = np.array(list(itertools.islice(xoutS_g1,0,(dp1+1),20)))
                    xoutG_lite = np.array(list(itertools.islice(xoutG_g1,0,(dp1+1),20)))
                    tout_lite = np.array(list(itertools.islice(tout_g1,0,(dp1+1),20)))


    if diag == 'off':
        np.savetxt(os.path.join(output_dir,xoutS_file),xoutS_lite,delimiter='\t')
        np.savetxt(os.path.join(output_dir,xoutG_file),xoutG_lite,delimiter='\t')
        np.savetxt(os.path.join(output_dir,tout_file),tout_lite,delimiter='\t')
    
    if diag == 'on':
        np.savetxt(os.path.join(output_dir,'g99_c'+str(cell_n+1)+'_xoutS.txt'),xoutS_new,delimiter='\t')
        # np.savetxt(os.path.join(output_dir,xoutG_file),xoutG_lite,delimiter='\t')
        np.savetxt(os.path.join(output_dir,'g99_c'+str(cell_n+1)+'_tout.txt'),tout_new,delimiter='\t')
            
            
            
    return results


with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
    if executor is not None:
        gen1a_begins = time.time()
        print("Starting gen1 analysis from rank", MY_RANK)
        results_g1 = list(executor.map(read_cell_g1, range(cellpop_g1), unordered=True))
        num_done = len(results_g1)
        gen1a_ends = time.time()
        gen1a_runtime = float(gen1a_ends - gen1a_begins)
        print("Finished", num_done,"gen1 analysis in",gen1a_runtime,"seconds.")
        
# if MY_RANK == 0:        
#     results_actual = np.array(results_g1)[np.where(results_g1)[0]]
    
#     if len(results_actual) != 0:
        
#         th_g2 = [r['th_g2'] for r in results_actual]
        
#         lin_g2 = [r['lin'] for r in results_actual]
    
#         th_g2 = [[th_g2[i]]+[th_g2[i]] for i in range(len(th_g2))]
#         th_g2 = [item for sublist in th_g2 for item in sublist]
        
#         lin_g2 = [[lin_g2[i]]+[lin_g2[i]] for i in range(len(lin_g2))]
#         lin_g2 = [item for sublist in lin_g2 for item in sublist]


# th_g2 = MPI.COMM_WORLD.bcast(th_g2, root = 0)
# lin_g2 = MPI.COMM_WORLD.bcast(lin_g2, root = 0)
#%% temp -  will need work
  

results_g1_all = None

if MY_RANK == 0:
    results_g1_all = results_g1
    

results_g1_all = MPI.COMM_WORLD.bcast(results_g1_all, root = 0)
  
results_actual = np.array(results_g1_all)[np.where(results_g1_all)[0]]

if len(results_actual) != 0:
    
    th_g2 = [r['th_g2'] for r in results_actual]
    
    lin_g2 = [r['lin'] for r in results_actual]

    th_g2 = [[th_g2[i]]+[th_g2[i]] for i in range(len(th_g2))]
    th_g2 = [item for sublist in th_g2 for item in sublist]
    
    lin_g2 = [[lin_g2[i]]+[lin_g2[i]] for i in range(len(lin_g2))]
    lin_g2 = [item for sublist in lin_g2 for item in sublist]
    
else:
    sys.exit("No division event detected at gen 1")



# np.savetxt(os.path.join(output_dose,'th_g2.txt'),np.array(th_g2),delimiter='\t')
# np.savetxt(os.path.join(output_dose,'lin_g2.txt'),np.array(lin_g2),delimiter='\t',fmt='%s')

#%% gen n while loop


lin_gn0 = lin_g2

th_gn0 = th_g2
   
cellpop_gn0 = len(th_g2)

g = 2


def cell_gn(cell_n,lin_gc,th_gc):
    
    # print("Running gen(%d) cell(%d) on rank %d" %(g,cell_n+1, MY_RANK))
    print("Running gen(%d) cell(%d) (lin(%s)) for (%d) hrs on rank %d" %(g,cell_n+1,lin_gc,th_gc, MY_RANK))
    
    cell_name = 'g'+str(g)+'_c'+str(cell_n+1)+'_lin_'+str(lin_gc)
    
    c0 = int(str(lin_gc).split('c')[-1])
    
    sp0 = np.loadtxt(os.path.join(output_dir,'g'+str(g-1)+'_c'+str(c0)+'_ic.txt'),delimiter='\t')
    
    xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th_gc,sp0,Vn,Vc,model,wd,omics_input,genereg_input)
    
    tout_all = tout_all + (th-th_gc)*3600

   
    
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutS.txt'),xoutS_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_xoutG.txt'),xoutG_all,delimiter='\t')
    np.savetxt(os.path.join(output_dir,str(cell_name)+'_tout.txt'),tout_all,delimiter='\t')



def read_cell_gn(cell_n):
    
    gx_cx = 'g'+str(g)+'_c'+str(int(cell_n+1))    
    
    outputs_ls = os.listdir(os.path.join(wd,'output',output_dir))
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),outputs_ls))
    xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
    tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    xoutG_file = list(filter(lambda x:x.endswith('xoutG.txt'),outputs_gc))[0]
    
    xoutS_all = np.loadtxt(os.path.join(output_dir,xoutS_file),delimiter='\t')
    xoutG_all = np.loadtxt(os.path.join(output_dir,xoutG_file),delimiter='\t')
    tout_all = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    
    if diag == 'off':
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
                results['lin'] = lin_gn_cell
                
                                
                np.savetxt(os.path.join(output_dir,gx_cx+'_ic.txt'),sp_gn_cell,delimiter='\t')
                
            if diag == 'off':

                xoutS_lite = np.array(list(itertools.islice(xoutS_all,0,(dp+1),20)))
                xoutG_lite = np.array(list(itertools.islice(xoutG_all,0,(dp+1),20)))
                tout_lite = np.array(list(itertools.islice(tout_all,0,(dp+1),20)))   

                

    if diag == 'off':
        np.savetxt(os.path.join(output_dir,xoutS_file),xoutS_lite,delimiter='\t')
        np.savetxt(os.path.join(output_dir,xoutG_file),xoutG_lite,delimiter='\t')
        np.savetxt(os.path.join(output_dir,tout_file),tout_lite,delimiter='\t')
            
            
    return results

while cellpop_gn0 > 0:
    
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            gen_n_begins = time.time()
            print("Submitting gen",g,"stage from rank", MY_RANK)
            res = executor.map(cell_gn, range(cellpop_gn0),lin_gn0,th_gn0, unordered=True)
            num_done = len(list(res))
            gen_n_ends = time.time()
            gen_n_runtime = float(gen_n_ends - gen_n_begins)
            print("Finished", num_done,"gen"+str(g), "tasks in",gen_n_runtime,"seconds.")
    
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            gen_na_begins = time.time()
            print("Starting gen",g ,"analysis from rank", MY_RANK)
            results_gn = list(executor.map(read_cell_gn, range(cellpop_gn0), unordered=True))
            num_done = len(results_gn)
            gen_na_ends = time.time()
            gen_na_runtime = float(gen_na_ends - gen_na_begins)
            print("Finished", num_done,"gen"+str(g), "analysis in",gen_na_runtime,"seconds.")
            
    results_gn_all = None

    if MY_RANK == 0:
        results_gn_all = results_gn
        
    
    results_gn_all = MPI.COMM_WORLD.bcast(results_gn_all, root = 0)
      
    results_gn_actual = np.array(results_gn_all)[np.where(results_gn_all)[0]]
    

        
    th_gn = [r['th_gn'] for r in results_gn_actual]
    
    lin_gn = [r['lin'] for r in results_gn_actual]

    th_gn = [[th_gn[i]]+[th_gn[i]] for i in range(len(th_gn))]
    th_gn = [item for sublist in th_gn for item in sublist]
    
    lin_gn = [[lin_gn[i]]+[lin_gn[i]] for i in range(len(lin_gn))]
    lin_gn = [item for sublist in lin_gn for item in sublist]
    
    cellpop_gn = len(th_gn)
    
    cellpop_gn0 = cellpop_gn
    
    if cellpop_gn0 > 0:
        
        np.savetxt(os.path.join(output_dose,'th_g'+str(g)+'.txt'),np.array(th_gn),delimiter='\t')
        np.savetxt(os.path.join(output_dose,'lin_g'+str(g)+'.txt'),np.array(lin_gn),delimiter='\t',fmt='%s')
        
        print("Division event detected at gen(%d)" %(g))
        g += 1
        
        lin_gn0 = lin_gn

        th_gn0 = th_gn
    else:
         print("No division event detected at gen(%d)" %(g))
        
sys.exit("Finishing simulation...")
        
#%% run gen2 manually
# th = 72
# lin_g2 = np.loadtxt(os.path.join(test_dir,'lin_g2.txt'),delimiter='\t',dtype=str)
# th_g2 = np.loadtxt(os.path.join(test_dir,'th_g2.txt'),delimiter='\t')
# cellpop_g2 = len(lin_g2)

# output_dir = test_dir

# g = 2

# for c in range(cellpop_g2):
#     cell_gn(c, lin_g2[c], th_g2[c])
    
    
#%% debug cell_gn

# def cell_gn(cell_n,lin_gc,th_gc):
    
#     print("Running gen(%d) cell(%d) on rank %d" %(g,cell_n+1, MY_RANK))
    
#     cell_name = 'g'+str(g)+'_c'+str(cell_n+1)+'_lin_'+str(lin_gc)
    
#     c0 = int(str(lin_gc).split('c')[-1])
    
#     sp0 = np.loadtxt(os.path.join(output_dir,'g'+str(g-1)+'_c'+str(c0)+'_ic.txt'),delimiter='\t')
    
#     xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th_gc,sp0,Vn,Vc,model,wd,omics_input,genereg_input)
    
#     tout_all = tout_all + (th-th_gc[cell_n])*3600

# cell_n = 0

# lin_gc = lin_g2[cell_n]
# th_gc = th_g2[cell_n]

# cell_name = 'g'+str(g)+'_c'+str(cell_n+1)+'_lin_'+str(lin_gc)

# c0 = int(str(lin_gc).split('c')[-1])
# sp0 = np.loadtxt(os.path.join(output_dir,'g'+str(g-1)+'_c'+str(c0)+'_ic.txt'),delimiter='\t')

# xoutS_all, xoutG_all, tout_all, flagA = RunSPARCED(flagD,th_gc,sp0,Vn,Vc,model,wd,omics_input,genereg_input)

# tout_all_test = tout_all + (th-th_gc[cell_n])*3600
