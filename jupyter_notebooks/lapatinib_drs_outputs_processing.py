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
import random

from mpi4py import MPI

MY_RANK = MPI.COMM_WORLD.Get_rank()


mpl.rcParams['figure.dpi'] = 300

parser = argparse.ArgumentParser(description='Input doses in uM')
parser.add_argument('--dose', metavar='dose', help='input laptinib dose in uM', default = 0.0)
parser.add_argument('--egf', metavar='egf', help='input EGF conc in nM', default = 100.0)
parser.add_argument('--cellpop', metavar='cellpop', help='starting cell population', default = 3)
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

STIMligs = [100.0,0,0,0,0,0,0.0]


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

#%% dirs 
# test_dir = '/media/arnab/Arnab/projects/sparced/sparced_initialization/sparced/output/lapatinib_drs_test_50/lapatinib_0.0'

# test_dir = '/media/arnab/Arnab/projects/sparced/sparced_initialization/sparced/output/lapatinib_drs_test_51/lapatinib_0.0'

# test_dir = '/media/arnab/Arnab/projects/sparced/sparced_initialization/sparced/output/lapatinib_drs_test_54/lapatinib_0.0'

# test_dir = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_r3_b0/lapatinib_0.0'

r3b0_dir = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_r3_b0/lapatinib_0.0'

test_dir = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_test_56/lapatinib_0.0'

test_dir2 = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_test_57/lapatinib_0.0'

test_dir3 = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_test_58/lapatinib_0.0'



#%% aux functions
import re
def read_tout(output_dir,tout_file):
    tout = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    tmax = max(tout)
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

# timecourse_tmax = get_tmax(output_dir)   

def timecourse_gc(output_dir,species, gx_cx,get_tmax=get_tmax,read_tout=read_tout,tx='default'):
    
    outputs_ls = os.listdir(os.path.join(wd,output_dir))
    
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),outputs_ls))
    
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

def get_lin(output_dir,gx_cx):
    
    outputs_ls = os.listdir(output_dir)    
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),outputs_ls))
    output_gc = list(filter(lambda x:x.endswith('_tout.txt'),outputs_gc))
    
    lin = re.search('lin_(.*)_tout.txt',output_gc[0]).group(1)
    
    lin = lin.split('c')[1:]
    
    return lin
    
def get_desc_g1(output_dir,cn):
    
    outputs_ls = os.listdir(output_dir)
    outputs_tout = list(filter(lambda x:x.endswith('_tout.txt'),outputs_ls))
    outputs_lin = list(filter(lambda x: 'lin_c'+str(cn) in x, outputs_tout))
    
    
    # g_all = list(filter(re.compile('g\d+').match, outputs_lin))
    
    gx_all = [re.search('g\d+',outputs_lin[l]).group(0) for l in range(len(outputs_lin))]
    
    g_all = [re.search('\d+',gx_all[l]).group() for l in range(len(gx_all))]
    
    g_all = np.array(g_all).astype('int')
    
    g_max = max(g_all)
    
    output_gmax = list(filter(lambda x:x.startswith('g'+str(g_max)),outputs_lin))[0]
    
    lin_gmax = re.search('lin_(.*)_tout.txt',output_gmax).group(1).split('c')[1:]
    
    
    
    c_gmax = re.search('(g)'+str(g_max)+'(_c\d+)',output_gmax).group()
    
    lin_gmax.append(c_gmax.split('_c')[1])
    
    return lin_gmax

def get_desc_g1_2(output_dir,cn):
    outputs_ls = os.listdir(output_dir)
    outputs_tout = list(filter(lambda x:x.endswith('_tout.txt'),outputs_ls))
    outputs_g2 = list(filter(lambda x:x.startswith('g2_c'),outputs_tout))

    outputs_lin = list(filter(lambda x: 'lin_c'+str(cn)+'_' in x, outputs_g2))
    desc_2 = [re.search('g\d+_c(.*)_lin',outputs_lin[m]).group(1) for m in range(len(outputs_lin))]
    
    desc_2.sort()
    
    return desc_2
    
    

    
def timecourse_lin(output_dir,species, g1_cx,legend='off',get_tmax=get_tmax,read_tout=read_tout,tx='default'):
    
    
    lin_g1 = get_desc_g1(output_dir,g1_cx)
    
    if 'timecourse_tmax' not in globals():
        tmax = get_tmax(output_dir)
    else:
        tmax = timecourse_tmax
    
    outputs_ls = os.listdir(os.path.join(wd,output_dir))
    
    x_t_lin = []
    tout_lin = []
    
    x_t_max = []
    
    
    
    for k in range(len(lin_g1)):
        gx_cx = 'g'+str(k+1)+'_c'+str(lin_g1[k])
    
        outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),outputs_ls))
    
        xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
        tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    


        x_s = np.loadtxt(os.path.join(output_dir,xoutS_file),delimiter='\t')
        tout_all = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')

    
        x_t = x_s[:, list(species_all).index(species)]
        
        x_t_lin.append(x_t)
        tout_lin.append(tout_all)
        x_t_max.append(max(x_t))
        
    ymax = max(np.array(x_t_max))
    
    gens = []
    for l in range(len(x_t_lin)):
        
        plt.plot(tout_lin[l]/3600, x_t_lin[l])
        plt.ylabel('species: '+str(species))
        plt.xlabel('time(h)')
        gens.append('gen '+str(l+1))
    plt.ylim(0, ymax*1.25)
    if type(tx)==str:
        plt.xlim(0, round(max(tout_all)/3600))
        plt.xlim(0, round(tmax/3600))
        
    elif type(tx)==int or type(tx)==float:
        plt.xlim(0,tx)
    
    plt.title("lineage: g1_c"+str(g1_cx))
    if legend == 'on':
        plt.legend(gens)

    plt.show   

def get_desc_gn(output_dir,gx_cx,g):
    
    gc = int(gx_cx.split('_')[0].split('g')[1])
    
    if gc == 1:
        
        desc_n = get_desc_g1_2(output_dir,int(gx_cx.split('c')[1]))
        
        desc_n.sort()
        
    else:
        outputs_ls = os.listdir(output_dir)
        
        outputs_tout = list(filter(lambda x:x.endswith('_tout.txt'),outputs_ls))
        
        p = get_lin(output_dir,gx_cx)
        
        lin_str = ''
        
        for k in range(len(p)):
            lin_str = lin_str+'c'+str(p[k])
        
        cn = gx_cx.split('c')[1]
        
        lin_str = lin_str+'c'+str(cn)
        
        desc_g_all = list(filter(lambda x:x.startswith('g'+str(g)),outputs_tout))
        
        desc_g = list(filter(lambda x:'lin_'+str(lin_str) in x,desc_g_all))
        
        desc_n = [re.search('g\d+_c(.*)_lin',desc_g[m]).group(1) for m in range(len(desc_g))]
        
        desc_n.sort()
    
    
    
    return desc_n


def tout_line(output_dir,gx_cx):
    outputs_ls = os.listdir(output_dir)
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx),outputs_ls))
    tout_file = list(filter(lambda x:x.endswith('_tout.txt'),outputs_gc))[0]
    
    tout = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    
    return (tout[0]/3600,tout[-1]/3600)


def lin_tree(output_dir,g1_cx):
    outputs_ls = os.listdir(output_dir)
    outputs_tout = list(filter(lambda x:x.endswith('_tout.txt'),outputs_ls))
    outputs_lin = list(filter(lambda x: 'lin_c'+str(g1_cx)+'c' in x or 'lin_c'+str(g1_cx)+'_' in x, outputs_tout))
    
    outputs_lin.sort()
    
    gx_all = [re.search('g\d+',outputs_lin[l]).group(0) for l in range(len(outputs_lin))]
    
    g_all = [re.search('\d+',gx_all[l]).group() for l in range(len(gx_all))]
    
    g_all = np.array(g_all).astype('int')
    
    g_max = int(max(g_all))
    
    h = 10
    
    inc = np.array([-1.15,1.15])
    
    h1 = h + inc
    
    red = random.random()
    green = random.random()
    blue = random.random()
            
    color = (red,green,blue)
    
    
    plt.plot((tout_line(test_dir,'g1_c'+str(g1_cx))[0],tout_line(test_dir,'g1_c'+str(g1_cx))[1]),(h,h), '-', c=color)
    plt.plot((tout_line(test_dir,'g1_c'+str(g1_cx))[1],tout_line(test_dir,'g1_c'+str(g1_cx))[1]),(h1[0],h1[1]), '-', c=color)
    
    g = 2
    c_des = get_desc_gn(output_dir,'g1_c'+str(g1_cx),2)
    

    hx = 0.72
    
    while g <= g_max:
        
        red = random.random()
        green = random.random()
        blue = random.random()
            
        color = (red,green,blue)
        
        hx = hx*hx
        
        
        c_des_new = []
        
        for cc in range(len(c_des)):
            
            # red = random.random()
            # green = random.random()
            # blue = random.random()
            
            # color = (red,green,blue)
            
            x1 = tout_line(test_dir,'g'+str(g)+'_c'+str(c_des[cc]))[0]
            x2 = tout_line(test_dir,'g'+str(g)+'_c'+str(c_des[cc]))[1]
            plt.plot((x1,x2),(h1[cc],h1[cc]), '-', c=color)
            
            if g != g_max:
                plt.plot((x2,x2),(h1[cc]+inc[0]*hx,h1[cc]+inc[1]*hx),'-', c=color)
         
            c_des_new.append(get_desc_gn(output_dir,'g'+str(g)+'_c'+str(c_des[cc]),g+1))
        
        c_des_new = np.array(c_des_new).flatten()
        
        c_des = list(c_des_new)
        
        h_new = [h1[i] + inc*hx for i in range(len(h1))]
        
        h1 = np.array(h_new).flatten()
        
        g += 1
    
    plt.title('lineage tree: g1_c'+str(g1_cx))
    plt.yticks([])
    plt.xlabel('time(hours)')
    plt.show()
            
    

    
    
    # g_all = list(filter(re.compile('g\d+').match, outputs_lin))
    

    
#%%

def cellpop_hist(output_dir):
    
    outputs_ls = os.listdir(output_dir)
    outputs_tout = list(filter(lambda x:x.endswith('_tout.txt'),outputs_ls))
    outputs_actual = list(filter(lambda x:not (x.startswith('g0_c') or x.startswith('g99_c')),outputs_tout))
    outputs_actual.sort()
    
    gx_all = [re.search('g\d+',outputs_actual[l]).group(0) for l in range(len(outputs_actual))]
    
    g_all = [re.search('\d+',gx_all[l]).group() for l in range(len(gx_all))]
    
    g_all = np.array(g_all).astype('int')
    
    g_max = int(max(g_all))
    
    cellpop_g = pd.Series()
    
    for g in range(g_max):
        if g < 5:
            cellpop_g['gen'+str(g+1)] = sum(g_all==(g+1))
        
    cellpop_g.plot(kind='bar',ylim=(0,500))

#%% debug - run gen1 analysis manually



#%% debug find_dp
def read_tmin(output_dir,tout_file):
    tout = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    tmin = min(tout)
    return (tmin)
def get_tmin(output_dir,read_tmin=read_tmin):
    
    outputs_ls = os.listdir(output_dir)
    tout_files = list(filter(lambda x:x.endswith('tout.txt'),outputs_ls))
                      
    results_tmin = []                 
                  
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results_tout = [executor.submit(read_tmin, output_dir, tout_file) for tout_file in tout_files]

    for f in concurrent.futures.as_completed(results_tout):
        results_tmin.append(f.result())
        
    results_tmin = np.array(results_tmin)
    
    tmin = min(results_tmin)
    
    return(tmin)

timecourse_tmin = get_tmin(output_dir)

    
def timecourse_gc2(output_dir,species, gx_cx,get_tmax=get_tmax,read_tout=read_tout,tx='default'):
    
    outputs_ls = os.listdir(os.path.join(wd,output_dir))
    
    outputs_gc = list(filter(lambda x:x.startswith(gx_cx+'_'),outputs_ls))
    
    xoutS_file = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_gc))[0]
    tout_file = list(filter(lambda x:x.endswith('tout.txt'),outputs_gc))[0]
    
    if 'timecourse_tmax' not in globals():
        tmax = get_tmax(output_dir)
    else:
        tmax = timecourse_tmax
        
    # if 'timecourse_tmin' not in globals():
    tmin = get_tmin(output_dir)
    # else:
    #     tmin = timecourse_tmin



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
        plt.xlim(round(tmin/3600), round(tmax/3600))
        
    elif type(tx)==int or type(tx)==float:
        plt.xlim(round(tmin/3600),tx)
        
    plt.title(gx_cx)

    plt.show
    
#%% debug find_dp - part 2
timecourse_gc2(r3b0_dir,'Mb','g99_c2')

xoutS_r3b0 = np.loadtxt(os.path.join(r3b0_dir,'g99_c2_xoutS.txt'),delimiter='\t')
tout_r3b0 = np.loadtxt(os.path.join(r3b0_dir,'g99_c2_tout.txt'),delimiter='\t')

data_r3b0 = xoutS_r3b0[:,list(species_all).index('Mb')]

p,_=find_peaks(data_r3b0,height=30)
b = (np.diff(np.sign(np.diff(data_r3b0))) > 0).nonzero()[0] + 1
dp_all = []
for i in range(len(p)):
    b2 = np.where(b>p[i])[0]
    if len(b2)!=0:
        dp_all.append(b[b2[0]])
    



def find_dp_all(xoutS,species_all=species_all):
    data = xoutS[:,list(species_all).index('Mb')]
    p,_ = find_peaks(data,height=30)
    b = (np.diff(np.sign(np.diff(data))) > 0).nonzero()[0] + 1
    
    dp_all = []
    for i in range(len(p)):
        b2 = np.where(b>p[i])[0]
        if len(b2)!=0:
            dp_all.append(b[b2[0]])
            
    dp_all_actual = list(np.array(dp_all)[data[dp_all]<1])
    
    return(dp_all_actual)

def find_dp(xoutS,tout,species_all=species_all):
    data = xoutS[:,list(species_all).index('Mb')]
    p,_ = find_peaks(data,height=30)
    b = (np.diff(np.sign(np.diff(data))) > 0).nonzero()[0] + 1
    
    if sum(b>p[0]) > 0:
    
        dp = int(b[b>p[0]][0])
        
    else:
        dp = np.nan
    
    return(dp)

#%% debug find_dp - part 3
p,_ = find_peaks(data_r3b0,height=30)
b = (np.diff(np.sign(np.diff(data_r3b0))) > 0).nonzero()[0] + 1

if len(b)!=0:
    b = np.array(b)[data_r3b0[b]<1]
if sum(b>p[0]) > 0:
    
    dp = int(b[b>p[0]][0])
        
else:
    dp = np.nan
#%%

test_dir63 = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_test_63/lapatinib_0.0'


xoutS_t63 = np.loadtxt(os.path.join(test_dir63,'g99_c20_xoutS.txt'),delimiter='\t')
tout_t63= np.loadtxt(os.path.join(test_dir63,'g99_c20_tout.txt'),delimiter='\t')

data_t63= xoutS_t63[:,list(species_all).index('Mb')]

p,_=find_peaks(data_t63,height=30)
b = (np.diff(np.sign(np.diff(data_t63))) > 0).nonzero()[0] + 1
dp_all = []
for i in range(len(p)):
    b2 = np.where(b>p[i])[0]
    if len(b2)!=0:
        dp_all.append(b[b2[0]])

dp_all_actual = list(np.array(dp_all)[data_t63[dp_all]<1])


#%%

def read_toutmax(output_dir,tout_file):
    tout = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    tmax = max(tout)
    return (tmax)



def get_tmax(output_dir,read_tout=read_tout):
    
    outputs_ls = os.listdir(output_dir)
    tout_files = list(filter(lambda x:x.endswith('tout.txt'),outputs_ls))
                      
    results_tmax = []                 
                  
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results_tout = [executor.submit(read_toutmax, output_dir, tout_file) for tout_file in tout_files]

    for f in concurrent.futures.as_completed(results_tout):
        results_tmax.append(f.result())
        
    results_tmax = np.array(results_tmax)
    
    tmax = max(results_tmax)
    
    return(tmax)

# timecourse_tmax = get_tmax(output_dir) 


def read_tout(output_dir,tout_file):
    tout = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    return (tout)

def read_xout(output_dir,xout_file):
    # print(xout_file)
    
    xoutS = np.loadtxt(os.path.join(output_dir,xout_file),delimiter='\t')
    if np.shape(np.shape(xoutS))[0] == 2:
        xparp = xoutS[:,list(species_all).index('PARP')]
        xcparp = xoutS[:,list(species_all).index('cPARP')]
        
    
        
        xout_result = np.concatenate((xparp.reshape((1,len(xparp))),xcparp.reshape((1,len(xcparp)))))
    elif np.shape(np.shape(xoutS))[0] == 1:
        xparp = xoutS[list(species_all).index('PARP')]
        xcparp = xoutS[list(species_all).index('cPARP')]
        
        xout_result = np.array([xparp,xcparp])
        
        xout_result = xout_result.reshape(2,1)
    # tmax = max(tout)
    else:
        xout_result = np.nan
    
    return (xout_result)

def read_xout_sp(output_dir,xout_file,sp):
    print(xout_file+'\n')
    xoutS = np.loadtxt(os.path.join(output_dir,xout_file),delimiter='\t')
    if np.shape(np.shape(xoutS))[0] == 2:
        xsp = xoutS[:,list(species_all).index(sp)]
    elif np.shape(np.shape(xoutS))[0] == 1:
        xsp = xoutS[list(species_all).index(sp)]
    else:
        xsp = np.nan
    return xsp
    

def drs_outputs(output_dir):
    outputs_ls = os.listdir(os.path.join(output_dir))
    
    outputs_ls = list(filter(lambda x:'g0_c' not in x,outputs_ls))
    outputs_ls = list(filter(lambda x:'g99_c' not in x,outputs_ls))
    
    outputs_ls.sort()
    
    outputs_tout = list(filter(lambda x:x.endswith('tout.txt'),outputs_ls))
    outputs_tout.sort()
    outputs_xout = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_ls))
    outputs_xout.sort()  
    n_cells = len(outputs_tout)
  
    # tmax = get_tmax(output_dir)

    
    results_tout = []
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
         results_tout_map = executor.map(read_tout, [output_dir]*len(outputs_tout), outputs_tout)
    
    for result in results_tout_map:
        results_tout.append(result)
        
    results_xout = []
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
          results_xout_map = executor.map(read_xout, [output_dir]*len(outputs_xout), outputs_xout)
      
    for result in results_xout_map:
        results_xout.append(result)
        
    
    tmax = 72*3600
    
    tout_starts = []
    for i in range(len(results_tout)):
        if results_tout[i].size>1:
            tout_starts.append(results_tout[i][0])
        elif results_tout[i].size == 1:
            tout_starts.append(results_tout[i])
         
    # tout_starts = [results_tout[i][0] for i in range(len(results_tout))]
    tout_starts = np.array(tout_starts)
    
    tout_ends = []
    for i in range(len(results_tout)):
        if results_tout[i].size>1:
            tout_ends.append(results_tout[i][-1])
        elif results_tout[i].size == 1:
            tout_ends.append(tmax)

    
    timepoints_all = np.concatenate((tout_starts,tout_ends))
    
    timepoints_all = np.unique(timepoints_all)
   
    cells_all = np.zeros((n_cells,len(timepoints_all)))
    

    for c in range(n_cells):
        s = np.where(timepoints_all == tout_starts[c])[0][0]
        e = np.where(timepoints_all == tout_ends[c])[0][0]
        
        cells_all[c,s:e+1] = np.ones(e+1-s)
          
    tout_deaths = np.ones(n_cells)*np.nan
  
    for c in range(n_cells):
    
        xs1 = results_xout[c]
    
        flagA = np.where(xs1[1]>xs1[0])[0]
        if np.shape(xs1)[1]>1:
            if len(flagA)!=0:
                td_idx = flagA[0]
                td = results_tout[c][td_idx]
                tout_deaths[c] = td
        elif np.shape(xs1)[1] == 1:
            if len(flagA)!=0:
                td = results_tout[c]
                tout_deaths[c] = td

        
        
    for k in range(len(tout_deaths)):
        if ~np.isnan(tout_deaths[k]):
            if tout_deaths[k] not in timepoints_all:
                timepoints_all = np.append(timepoints_all,tout_deaths[k])
                timepoints_all.sort()
                ip_idx = np.where(timepoints_all > tout_deaths[k])[0][0]
                cells_all = np.insert(cells_all, ip_idx, copy.deepcopy(cells_all[:,ip_idx-1]), axis = 1)
                cells_all[k,ip_idx:] = 0.0
            elif tout_deaths[k] in timepoints_all:
                ip_idx = np.where(timepoints_all == tout_deaths[k])[0][0]
                cells_all[k,ip_idx:] = 0.0
                
    timecourse_cellpop = [sum(cells_all[:,t]) for t in range(len(timepoints_all))]
    
    return timecourse_cellpop, timepoints_all, tout_deaths
                
#%%
test_dir = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_r4/lapatinib_0.0'        



#%%
drs_dir = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_r4'

ls_drs_dir = os.listdir(drs_dir)

ls_drs_dir.sort()

doses = [float(dose.split('_')[1]) for dose in ls_drs_dir]

doses.sort()

result_doses = {}

#%%
for dose_dir in ls_drs_dir:
    tcp, tpa, tds = drs_outputs(os.path.join(drs_dir,dose_dir))
    result_dose = {}
    result_dose['timecourse_cellpop'] = tcp
    result_dose['timepoints_all'] = tpa
    result_dose['tout_deaths'] = tds
    
    result_doses[dose_dir] = result_dose

#%%
import pickle
f = open(os.path.join(wd,"result_doses.pkl"),"wb")
pickle.dump(result_doses,f)
f.close()

#%%
file_to_read = open(os.path.join(wd,'result_doses.pkl'),'rb')
result_dic = pickle.load(file_to_read)

#%%

for dose_dir in ls_drs_dir:
    tcp = result_doses[dose_dir]['timecourse_cellpop']
    tps = result_doses[dose_dir]['timepoints_all']
    
    plt.plot(tps/3600,tcp)


plt.title('cell population over time')
plt.xlabel('time (hrs)')
plt.ylabel('number of cells')
plt.legend(ls_drs_dir)
plt.show()

#%%
import pickle

file_to_read = open(os.path.join(wd,'result_doses.pkl'),'rb')
result_doses_output = pickle.load(file_to_read)
result_doses = result_doses_output

#%%


doses = [float(dose.split('_')[1]) for dose in ls_drs_dir]

doses.sort()

doses = np.array(doses)

# tp = 71.5 # timepoints of interest(hours)

tp = 60

cp_t = []

for d in range(len(doses)):

    tp_dose_all = result_doses[ls_drs_dir[d]]['timepoints_all']
    cellpop_t = result_doses[ls_drs_dir[d]]['timecourse_cellpop']
    
    tp_s = tp*3600
    
    if tp_s in tp_dose_all:
        cp_d = cellpop_t[np.where(tp_dose_all == tp_s)[0][0]]
    else:
        a1 = cellpop_t[np.where(tp_dose_all<tp_s)[0][-1]]
        a2 = cellpop_t[np.where(tp_dose_all>tp_s)[0][0]]
               
        cp_d = (a1+a2)/2
    cp_t.append(cp_d)
    
cp_t = np.array(cp_t)

        
plt.scatter(doses,cp_t/cp_t[0])
plt.xscale('log')
plt.xlim(1e-4,1e2)

plt.title('Lapatinib response at '+str(tp)+' hrs')
plt.ylabel('Relative cell population')
plt.xlabel('lapatinib dose (uM)')

plt.show()

#%% relative species concentrations vs dose

sp = 'ppAKT'

tp_h = 71.5

tp_s = tp_h*3600

x_t_doses = []

for d in range(len(doses)):

    output_dir = os.path.join(drs_dir,ls_drs_dir[d])


# output_dir = test_dir
    
    outputs_ls = os.listdir(os.path.join(wd,output_dir))
    
    
    outputs_ls = list(filter(lambda x:'g0_c' not in x,outputs_ls))
    outputs_ls = list(filter(lambda x:'g99_c' not in x,outputs_ls))
    
    outputs_ls.sort()
    
    outputs_tout = list(filter(lambda x:x.endswith('tout.txt'),outputs_ls))
    outputs_tout.sort()
    outputs_xout = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_ls))
    outputs_xout.sort()
    
    n_cells = len(outputs_tout)
    
    results_tout = []
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
         # results_tout_tpe = [executor.submit(read_tout, output_dir, tout_file) for tout_file in outputs_tout]
         results_tout_map = executor.map(read_tout, [output_dir]*len(outputs_tout), outputs_tout)
    
    for result in results_tout_map:
        results_tout.append(result)
        
    tmax = 72*3600
        
    tout_starts = []
    for i in range(len(results_tout)):
        if results_tout[i].size>1:
            tout_starts.append(results_tout[i][0])
        elif results_tout[i].size == 1:
            tout_starts.append(results_tout[i])
         
    # tout_starts = [results_tout[i][0] for i in range(len(results_tout))]
    tout_starts = np.array(tout_starts)
    
    tout_ends = []
    for i in range(len(results_tout)):
        if results_tout[i].size>1:
            tout_ends.append(results_tout[i][-1])
        elif results_tout[i].size == 1:
            tout_ends.append(tmax)
    
    
    # tout_ends = [results_tout[i][-1] for i in range(len(results_tout))]
    tout_ends = np.array(tout_ends)
        
    cells2include = []
    
    for c in range(n_cells):
        if results_tout[c].size > 1:
            if tp_s > tout_starts[c] or tp_s == tout_starts[c]:
                if tp_s < tout_ends[c] or tp_s == tout_ends[c]:
                    cells2include.append(c)
    
    cells2include = np.array(cells2include)
    
    results_xout_sp = []
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
          # results_xout_tpe = [executor.submit(read_xout, output_dir, xout_file) for xout_file in outputs_xout]
          results_xout_sp_map = executor.map(read_xout_sp, [output_dir]*len(cells2include), list(np.array(outputs_xout)[cells2include]),[sp]*len(cells2include))
    
    # for f in concurrent.futures.as_completed(results_xout_tpe):
    #       results_xout.append(f.result())
    
    for result in results_xout_sp_map:
        results_xout_sp.append(result)
        
    x_t_cells = []
    
    for l in range(len(cells2include)):
        tout_c = results_tout[cells2include[l]]
        if tp_s not in tout_c:
            tp1 = np.where(tp_s>tout_c)[0][-1]
            tp2 = np.where(tp_s<tout_c)[0][0]
            x_t_cells.append((results_xout_sp[l][tp1]+results_xout_sp[l][tp2])/2)
        elif tp_s in tout_c:
            tp_idx = np.where(tout_c == tp_s)[0][0]
            x_t_cells.append(results_xout_sp[l][tp_idx])
# kk = [results_tout[cells2include[i]].size for i in range(len(cells2include))]
# ll = [results_xout_sp[i].size for i in range(len(cells2include))]
    x_t_doses.append(x_t_cells)

#%%0

x_med = [np.median(x_t_doses[d]) for d in range(len(doses))]

plt.scatter(doses,x_med/x_med[0])

plt.xscale('log')
plt.xlim(1e-4,1e2)

plt.title('Relative median '+sp+' at '+str(tp_h)+' hrs')
plt.ylabel(sp)
plt.xlabel('lapatinib dose (uM)')

plt.show()





#%%

# test_dir63 = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_test_63/lapatinib_0.0'
# test_dir = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_r4/lapatinib_0.0025'  
test_dir = '/media/arnab/bigchungus/projects/ccle_egf_drs/SPARCED/output/lapatinib_drs_r4/lapatinib_8.0'  


output_dir = test_dir

outputs_ls = os.listdir(os.path.join(wd,output_dir))


outputs_ls = list(filter(lambda x:'g0_c' not in x,outputs_ls))
outputs_ls = list(filter(lambda x:'g99_c' not in x,outputs_ls))

outputs_ls.sort()

outputs_tout = list(filter(lambda x:x.endswith('tout.txt'),outputs_ls))
outputs_tout.sort()
outputs_xout = list(filter(lambda x:x.endswith('xoutS.txt'),outputs_ls))
outputs_xout.sort()

n_cells = len(outputs_tout)

tout_all = []

# checklist = []

# for i in range(999):
#     rr = random.randint(0,len(outputs_xout)-1)
#     checklist.append(outputs_xout[rr][:-9] == outputs_tout[rr][:-9])

# tmax = get_tmax(output_dir)



def read_tout(output_dir,tout_file):
    tout = np.loadtxt(os.path.join(output_dir,tout_file),delimiter='\t')
    # tmax = max(tout)
    return (tout)

results_tout = []

with concurrent.futures.ThreadPoolExecutor() as executor:
     # results_tout_tpe = [executor.submit(read_tout, output_dir, tout_file) for tout_file in outputs_tout]
     results_tout_map = executor.map(read_tout, [output_dir]*len(outputs_tout), outputs_tout)


# for f in concurrent.futures.as_completed(results_tout_tpe):
#      results_tout.append(f.result())

for result in results_tout_map:
    results_tout.append(result) 

tmax = 72*3600





#%%       
def read_xout(output_dir,xout_file):
    print(xout_file+' \n')
    
    xoutS = np.loadtxt(os.path.join(output_dir,xout_file),delimiter='\t')
    if np.shape(np.shape(xoutS))[0] == 2:
        xparp = xoutS[:,list(species_all).index('PARP')]
        xcparp = xoutS[:,list(species_all).index('cPARP')]
        
    
        
        xout_result = np.concatenate((xparp.reshape((1,len(xparp))),xcparp.reshape((1,len(xcparp)))))
    elif np.shape(np.shape(xoutS))[0] == 1:
        xparp = xoutS[list(species_all).index('PARP')]
        xcparp = xoutS[list(species_all).index('cPARP')]
        
        xout_result = np.array([xparp,xcparp])
        
        xout_result = xout_result.reshape(2,1)
    # tmax = max(tout)
    else:
        xout_result = np.nan
    
    return (xout_result)

# x1 = np.loadtxt(os.path.join(output_dir,outputs_xout[45]),delimiter='\t')

# x2 = x1[:,list(species_all).index('PARP')]
# x3 = x1[:,list(species_all).index('cPARP')]
#%% debug

# xout_file = 'g6_c9_lin_c6c20c69c109c121_xoutS.txt'
xout_file = 'g1_c18_xoutS.txt'


xoutS = np.loadtxt(os.path.join(output_dir,xout_file),delimiter='\t')
xparp = xoutS[list(species_all).index('PARP')]
xcparp = xoutS[list(species_all).index('cPARP')]
xout_result = np.array([xparp,xcparp])

xout_result = xout_result.reshape(2,1)

xout_file2 = 'g5_c830_lin_c87c194c370c689_xoutS.txt'
xs = np.loadtxt(os.path.join(output_dir,xout_file2),delimiter='\t')


#%%
results_xout = []

with concurrent.futures.ThreadPoolExecutor() as executor:
      # results_xout_tpe = [executor.submit(read_xout, output_dir, xout_file) for xout_file in outputs_xout]
      results_xout_map = executor.map(read_xout, [output_dir]*len(outputs_xout), outputs_xout)

# for f in concurrent.futures.as_completed(results_xout_tpe):
#       results_xout.append(f.result())

for result in results_xout_map:
    results_xout.append(result)
    
# results_xout_tuple = tuple(results_xout)
     
#%%
tmax = 72*3600

tout_starts = []
for i in range(len(results_tout)):
    if results_tout[i].size>1:
        tout_starts.append(results_tout[i][0])
    elif results_tout[i].size == 1:
        tout_starts.append(results_tout[i])
     
# tout_starts = [results_tout[i][0] for i in range(len(results_tout))]
tout_starts = np.array(tout_starts)

tout_ends = []
for i in range(len(results_tout)):
    if results_tout[i].size>1:
        tout_ends.append(results_tout[i][-1])
    elif results_tout[i].size == 1:
        tout_ends.append(tmax)


# tout_ends = [results_tout[i][-1] for i in range(len(results_tout))]
tout_ends = np.array(tout_ends)

timepoints_all = np.concatenate((tout_starts,tout_ends))

timepoints_all = np.unique(timepoints_all)

# timepoints_all = np.concatenate((np.array(tout_starts).reshape((1,len(tout_starts))),np.array(tout_ends)).reshape((1,len(tout_ends))),1)

cells_all = np.zeros((n_cells,len(timepoints_all)))

# c = 4

for c in range(n_cells):
    s = np.where(timepoints_all == tout_starts[c])[0][0]
    e = np.where(timepoints_all == tout_ends[c])[0][0]
    
    cells_all[c,s:e+1] = np.ones(e+1-s)
    

timecourse_cellpop = [sum(cells_all[:,t]) for t in range(len(timepoints_all))]

plt.plot(timepoints_all/3600,timecourse_cellpop)
plt.title('cell population over time')
plt.xlabel('time (hrs)')
plt.ylabel('number of cells')

#%% 
import copy

results_xout_test = copy.deepcopy(results_xout)

# results_xout_test = [x[:] for x in results_xout] #results_xout.copy()



# results_xout_test = list(results_xout_tuple)

# results_xout_test = []

# for item in results_xout:
#     results_xout_test.append(item)

# results_xout_test.extend(results_xout)

cells_all_test = copy.deepcopy(cells_all)
timepoints_all_test = copy.deepcopy(timepoints_all)

#%% kill cells

k_cells = np.random.randint(0,n_cells,size=int(n_cells*.15))

k_cells = np.unique(k_cells)

k_cells.sort()

for k in range(len(k_cells)):
    
    x_k = results_xout_test[k_cells[k]]
    td = np.random.randint(0,np.shape(x_k)[1])
    x_k[1,td:] = max(x_k[0])*2
    
    results_xout_test[k_cells[k]] = x_k
    
#%%
tout_deaths = np.ones(n_cells)*np.nan


for c in range(n_cells):

    xs1 = results_xout_test[c]

    flagA = np.where(xs1[1]>xs1[0])[0]
    if np.shape(xs1)[1]>1:
        if len(flagA)!=0:
            td_idx = flagA[0]
            td = results_tout[c][td_idx]
            tout_deaths[c] = td
    elif np.shape(xs1)[1] == 1:
        if len(flagA)!=0:
            td = results_tout[c]
            tout_deaths[c] = td

 
#%%
# np.where(~np.isnan(tout_deaths))[0] == k_cells

flagA_cells = np.where(~np.isnan(tout_deaths))[0] 
cells_all_test = copy.deepcopy(cells_all)
timepoints_all_test = copy.deepcopy(timepoints_all)

for k in range(len(tout_deaths)):
    if ~np.isnan(tout_deaths[k]):
        if tout_deaths[k] not in timepoints_all_test:
            timepoints_all_test = np.append(timepoints_all_test,tout_deaths[k])
            timepoints_all_test.sort()
            ip_idx = np.where(timepoints_all_test > tout_deaths[k])[0][0]
            cells_all_test = np.insert(cells_all_test, ip_idx, copy.deepcopy(cells_all_test[:,ip_idx-1]), axis = 1)
            cells_all_test[k,ip_idx:] = 0.0
        elif tout_deaths[k] in timepoints_all_test:
            ip_idx = np.where(timepoints_all_test == tout_deaths[k])[0][0]
            cells_all_test[k,ip_idx:] = 0.0
            

#%%

timecourse_cellpop_2 = [sum(cells_all_test[:,t]) for t in range(len(timepoints_all_test))]

plt.plot(timepoints_all_test/3600,timecourse_cellpop_2)
plt.title('cell population over time')
plt.xlabel('time (hrs)')
plt.ylabel('number of cells')