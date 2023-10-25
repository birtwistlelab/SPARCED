# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 22:58:16 2023

@author: Arnab
"""
# import required libraries

import pickle
import copy
import re
import random
import itertools
import math

#%%
import os
import sys

import libsbml
import numpy as np
import pandas as pd
from scipy.stats import percentileofscore

from Bio import Phylo

from io import StringIO

from scipy.interpolate import interp1d
from scipy.stats import percentileofscore
from scipy.stats import gaussian_kde
import math
import seaborn as sns
import itertools
import plotly.figure_factory as ff
import plotly.io as pio

#%%
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.sans-serif'] = ['Arial']


#%%
cd = os.getcwd()
wd = os.path.dirname(cd)
sys.path.append(os.path.join(wd,'bin'))

sbml_file = "SPARCED.xml"


sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(os.path.join(wd,sbml_file))
sbml_model = sbml_doc.getModel()

species_all = [str(x.getId()) for x in list(sbml_model.getListOfSpecies())]



#%%

output_dir = 'D:\\projects\\sparced_mpi_pd\\sparced_5\\output\\testmpi_cellpop50\\trame_EC_0.0'
output_dir = 'D:\\projects\\sparced_mpi_pd\\sparced_5\\output\\local_test01\\trame_EC_0.0'
output_dir = 'D:\\projects\\sparced_mpi_pd\\sparced_5\\output\\local_test07\\trame_EC_0.0'

#%% cellpop_drs_v13

output_dir = 'D:\\projects\\sparced_mpi_pd\\sparced_5\\output\\local_test09\\lapat_EC_0.0'
output_dir = 'D:\\projects\\sparced_mpi_pd\\sparced_5\\output\\local_test09\\lapat_EC_1.0'



#%%
sim_name = 'drs_lapat_rep2'

drug = 'lapat_EC'

output_dir_main = os.path.join(wd,'output')

output_dir_sim = os.path.join(output_dir_main,sim_name)

dir_doses_all = os.listdir(output_dir_sim)

doses_all = [float(x.split('_')[-1]) for x in dir_doses_all]

doses_all.sort()

dir1 = os.path.join(output_dir_sim,drug+'_'+str(doses_all[4]))

# output_all = load_outputs(dir1)

#%%

# define class for reading dose response outputs

class drs_dict():
    def __init__(self,main,drug,rep,dose_level):
        self.main = main
        self.drug = drug
        self.rep = rep
        self.dose_level = dose_level
        self.path_exp = os.path.join(main,'drs_'+drug)
        self.path_reps = os.listdir(self.path_exp)
        self.reps = [int(list(filter(str.isdigit, x.split('_')[-1]))[0]) for x in self.path_reps]
        self.path_rep = os.path.join(self.path_exp,self.path_exp,'drs_'+str(drug)+'_rep'+str(rep))
        self.path_doses = os.listdir(self.path_rep)
        self.doses = [float(x.split('_')[-1]) for x in self.path_doses]
        self.doses.sort()
        
        self.path_dose = os.path.join(self.path_rep,str(drug)+'_EC_'+str(self.doses[dose_level]))
        
        dose_files = os.listdir(self.path_dose)
        self.results = {}
        
        for i,file in enumerate(dose_files):
            with open(os.path.join(self.path_dose,file),'rb') as f:
                self.results['output_g'+str(i+1)] = pickle.load(f)
        
        self.lin_all = {}
        self.n_gen = len(self.results)
        
        for g in range(2,self.n_gen+1):
            lin_g = {}
            n_cells = len(self.results['output_g'+str(g)])
            for c in range(1,n_cells+1):
                lin_g[str(c)] = self.results['output_g'+str(g)][str(c)]['output']['lin'].split('c')[1:]
                
            self.lin_all['g'+str(g)] = lin_g
        
    def get_desc(self,g1cx): #find all gen1 cell descendents
        desc_all = {}
    
        for g in range(2,self.n_gen+1):
            desc_g = []
            n_cells = len(self.lin_all['g'+str(g)])
            for c in range(1,n_cells+1):
                if self.lin_all['g'+str(g)][str(c)][0] == str(g1cx):
                    desc_g.append(c)
            desc_all['g'+str(g)] = desc_g
    
        return desc_all
    
    def timecourse_lin(self,g1cx,sp): #cell lineage timecourse for species
        plt.rcParams["lines.linewidth"]=2.0

        colors = ['blue','orange','green','red','purple','brown','pink','gray','olive','cyan']  
        
        # sp = 'cPARP'
        
        desc_all = self.get_desc(g1cx)
        
        tout_g1 = self.results['output_g1'][str(g1cx)]['output']['tout']
        xout_g1 = self.results['output_g1'][str(g1cx)]['output']['xoutS']
        xout_sp = xout_g1[:,list(species_all).index(sp)]
        
        xouts = []
        touts = []
        xouts.append(xout_sp)
        touts.append(tout_g1)
        
        plt.plot(tout_g1/3600,xout_sp,c=colors[0])
        
        for g in range(2,self.n_gen+1):
            
            cells_g = desc_all['g'+str(g)]
            
            if len(cells_g)>0:
                for c in cells_g:
                    tout_c = self.results['output_g'+str(g)][str(c)]['output']['tout']
                    xout_c = self.results['output_g'+str(g)][str(c)]['output']['xoutS']
                    xout_sp_c = xout_c[:,list(species_all).index(sp)]
                    xouts.append(xout_sp_c)
                    touts.append(tout_c)
                    plt.plot(tout_c/3600,xout_sp_c,colors[g-1])
        
        xmax = max(np.array([max(x) for x in xouts]))
        tmax = max(np.array([max(tp/3600) for tp in touts]))
        
        plt.xlim(0,tmax*1.15)
        plt.ylim(0,xmax*1.25)
        plt.xlabel('Time (hours)',font='Arial',fontsize=15)
        plt.ylabel(sp)
        
        plt.show()        
        
    def timecourse_lin_obs(self,g1cx,obs_formula): #cell lineage timecourse for observable
        plt.rcParams["lines.linewidth"]=2.0

        colors = ['blue','orange','green','red','purple','brown','pink','gray','olive','cyan']  
        
        # sp = 'cPARP'
        global x_obs
        desc_all = self.get_desc(g1cx)
        
        tout_g1 = self.results['output_g1'][str(g1cx)]['output']['tout']
        xs = self.results['output_g1'][str(g1cx)]['output']['xoutS']
        
    
        sp_obs = re.findall(r'[a-zA-Z]\w*',obs_formula)
        sp_obs = list(np.unique(sp_obs))
        
        if 'e' in sp_obs:
            sp_obs.remove('e')
        
        
        for i in range(len(sp_obs)):
            exec(f"{sp_obs[i]} = xs[:,species_all.index('{sp_obs[i]}')]", globals(), locals())
        
        exec(f"x_obs = {obs_formula}")
        # print(x_obs)
        # xout_sp = xout_g1[:,list(species_all).index(sp)]
        
        obs_all = []
        touts = []
        # print(locals()['x_obs'])
        obs_all.append(locals()['x_obs'])
        touts.append(tout_g1)
        
        plt.plot(tout_g1/3600,locals()['x_obs'],c=colors[0])
        
        for g in range(2,self.n_gen+1):
            
            cells_g = desc_all['g'+str(g)]
            
            if len(cells_g)>0:
                for c in cells_g:
                    tout_c = self.results['output_g'+str(g)][str(c)]['output']['tout']
                    xs = self.results['output_g'+str(g)][str(c)]['output']['xoutS']
                    for i in range(len(sp_obs)):
                        exec(f"{sp_obs[i]} = xs[:,species_all.index('{sp_obs[i]}')]")
                    
                    exec(f"x_obs = {obs_formula}")                
                    
    
                    obs_all.append(locals()['x_obs'])
                    touts.append(tout_c)
                    plt.plot(tout_c/3600,locals()['x_obs'],colors[g-1])
        
        xmax = max(np.array([max(x) for x in obs_all]))
        tmax = max(np.array([max(tp/3600) for tp in touts]))
        
        plt.xlim(0,tmax*1.15)
        plt.ylim(0,xmax*1.25)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel('Time (h)')
        # plt.ylabel(sp)
        
        plt.show()
        
    def pop_dyn(self): # cell population over time
        
        n_cells = sum([len(self.results[x]) for x in self.results.keys()])
        
        tout_starts = []
        tout_ends = []
        
        tout_deaths = np.ones(n_cells)*np.nan
        cell_idx = []
        cell_count = 0
        
    
        
        for g in range(self.n_gen):
            output_gen = self.results['output_g'+str(g+1)]
            n_cells_g = len(output_gen)
        
            for c in range(n_cells_g):
                cell_idx.append((g+1,c+1))
                xout = output_gen[str(c+1)]['output']['xoutS']
                if np.shape(np.shape(xout))[0] == 2:
                    xparp = xout[:,list(species_all).index('PARP')]
                    xcparp = xout[:,list(species_all).index('cPARP')]        
                    xout_result = np.concatenate((xparp.reshape((1,len(xparp))),xcparp.reshape((1,len(xcparp)))))
                elif np.shape(np.shape(xout))[0] == 1:
                    xparp = xout[list(species_all).index('PARP')]
                    xcparp = xout[list(species_all).index('cPARP')]
                    
                    xout_result = np.array([xparp,xcparp])
                    
                    xout_result = xout_result.reshape(2,1)
                
                xs1 = xout_result
        
                flagA = np.where(xs1[1]>xs1[0])[0]
        
        
                
                tout = output_gen[str(c+1)]['output']['tout']
                if tout.size>1:
                    tout_starts.append(tout[0])
                    tout_ends.append(tout[-1])
                    if len(flagA)!=0:
                        td_idx = flagA[0]
                        td = tout[td_idx]
                        tout_deaths[cell_count] = td
                        # print(td,g+1,c+1)
                        # tout_deaths[g+1,c+1] = td
                    
                elif tout.size == 1:
                    tout_starts.append(tout[0])
                    tout_ends.append(tout[0])
                    if len(flagA)!=0:
                        td = tout
                        tout_deaths[cell_count] = td
                        # print(td,g+1,c+1)
                        # tout_deaths[g+1,c+1] = td
                cell_count +=1
        
        timepoints_all = np.concatenate((tout_starts,tout_ends))
        
        timepoints_all = np.unique(timepoints_all)
        
        cells_all = np.zeros((n_cells,len(timepoints_all)))
        
        for c in range(n_cells):
            s = np.where(timepoints_all == tout_starts[c])[0][0]
            e = np.where(timepoints_all == tout_ends[c])[0][0]
            
            cells_all[c,s:e+1] = np.ones(e+1-s)
        
             
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

    def pop_dyn_obs(self,obs_formula): #population level observables
        
        n_cells = sum([len(self.results[x]) for x in self.results.keys()])
        
        tout_starts = []
        tout_ends = []
        
        tout_deaths = np.ones(n_cells)*np.nan
        cell_idx = []
        cell_count = 0
        
    
        
        for g in range(self.n_gen):
            output_gen = self.results['output_g'+str(g+1)]
            n_cells_g = len(output_gen)
        
            for c in range(n_cells_g):
                cell_idx.append((g+1,c+1))
                xout = output_gen[str(c+1)]['output']['xoutS']
                if np.shape(np.shape(xout))[0] == 2:
                    xparp = xout[:,list(species_all).index('PARP')]
                    xcparp = xout[:,list(species_all).index('cPARP')]        
                    xout_result = np.concatenate((xparp.reshape((1,len(xparp))),xcparp.reshape((1,len(xcparp)))))
                elif np.shape(np.shape(xout))[0] == 1:
                    xparp = xout[list(species_all).index('PARP')]
                    xcparp = xout[list(species_all).index('cPARP')]
                    
                    xout_result = np.array([xparp,xcparp])
                    
                    xout_result = xout_result.reshape(2,1)
                
                xs1 = xout_result
        
                flagA = np.where(xs1[1]>xs1[0])[0]
        
        
                
                tout = output_gen[str(c+1)]['output']['tout']
                if tout.size>1:
                    tout_starts.append(tout[0])
                    tout_ends.append(tout[-1])
                    if len(flagA)!=0:
                        td_idx = flagA[0]
                        td = tout[td_idx]
                        tout_deaths[cell_count] = td
                        # print(td,g+1,c+1)
                        # tout_deaths[g+1,c+1] = td
                    
                elif tout.size == 1:
                    tout_starts.append(tout[0])
                    tout_ends.append(tout[0])
                    if len(flagA)!=0:
                        td = tout
                        tout_deaths[cell_count] = td
                        # print(td,g+1,c+1)
                        # tout_deaths[g+1,c+1] = td
                cell_count +=1
        
        timepoints_all = np.concatenate((tout_starts,tout_ends))
        
        timepoints_all = np.unique(timepoints_all)
        
        cells_all = np.zeros((n_cells,len(timepoints_all)))
        
        for c in range(n_cells):
            s = np.where(timepoints_all == tout_starts[c])[0][0]
            e = np.where(timepoints_all == tout_ends[c])[0][0]
            
            cells_all[c,s:e+1] = np.ones(e+1-s)
        
             
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
        
        pop_obs = []
        
        
        sp_obs = re.findall(r'[a-zA-Z]\w*',obs_formula)
        sp_obs = list(np.unique(sp_obs))
        
        if 'e' in sp_obs:
            sp_obs.remove('e')
        
        for g in range(self.n_gen):
            output_gen = self.results['output_g'+str(g+1)]
            n_cells_g = len(output_gen)
        
            for c in range(n_cells_g):
                tout_c = output_gen[str(c+1)]['output']['tout']
                xs  = output_gen[str(c+1)]['output']['xoutS']

                if np.shape(np.shape(xout))[0] == 2:
                    xparp = xs[:,list(species_all).index('PARP')]
                    xcparp = xs[:,list(species_all).index('cPARP')]        
                    xout_result = np.concatenate((xparp.reshape((1,len(xparp))),xcparp.reshape((1,len(xcparp)))))
                elif np.shape(np.shape(xout))[0] == 1:
                    xparp = xs[list(species_all).index('PARP')]
                    xcparp = xs[list(species_all).index('cPARP')]
                    
                    xout_result = np.array([xparp,xcparp])
                    
                    xout_result = xout_result.reshape(2,1)
                
                xs1 = xout_result
        
                flagA = np.where(xs1[1]>xs1[0])[0]
                
                if tout_c.size>1:
                    tout_start = tout_c[0]
                    tout_end = tout_c[-1]
                    if len(flagA)!=0:
                        td_idx = flagA[0]
                        td = tout_c[td_idx]
                        tout_end = td
                        # print(td,g+1,c+1)
                        # tout_deaths[g+1,c+1] = td
                    
                elif tout_c.size == 1:
                    tout_start = tout[0]
                    tout_end = tout[0]
                    
                tp_fill = timepoints_all[np.logical_and(tout_start<=timepoints_all,timepoints_all<=tout_end)]
                # xs_obs = np.zeros(len(timepoints_all))
                xs_obs = np.ones(len(timepoints_all))*np.nan
                
                #xs_obs[np.isin(timepoints_all,[tp_fill])] = x_obs_fill

                        # print(td,g+1,c+1)
                        # tout_deaths[g+1,c+1] = td

                
                
                
                for i in range(len(sp_obs)):
                    exec(f"{sp_obs[i]} = xs[:,species_all.index('{sp_obs[i]}')]")
                exec(f"x_obs = {obs_formula}")
                
                cell_obs = {}
                
                # tout_obs = np.intersect1d(timepoints_all,tout_c)
                
                if tout_c.size>1:
                    obs_interp = interp1d(tout_c,locals()['x_obs'])
                    x_obs_fill = obs_interp(tp_fill)
                    
                elif tout_c.size == 1:
                    x_obs_fill = locals()['x_obs']

                
                # x_obs_fill = obs_interp(tp_fill)
                
                xs_obs[np.isin(timepoints_all,[tp_fill])] = x_obs_fill
                
                # cell_obs['x_obs'] = locals()['x_obs']
                cell_obs['x_obs'] = xs_obs
                cell_obs['tout'] = timepoints_all
                
                # pop_obs['g'+str(g+1)+'c'+str(c+1)] = cell_obs
                
                pop_obs.append(xs_obs)
        
        pop_obs_final = np.array(pop_obs)
            
        return cells_all, timepoints_all, tout_deaths, pop_obs_final      

        
    def get_desc_gc(self,gn,cn): # find daughter cell for any gen/cell
        desc_gn1 = []
        gn1 = gn+1
        if 'g'+str(gn1) in self.lin_all.keys():
            gn_cells = self.lin_all['g'+str(gn1)]
            for c in range(1,len(gn_cells)+1):
                if gn_cells[str(c)][gn-1] == str(cn):
                    desc_gn1.append(c)                
        return desc_gn1
        
    def get_len_gc(self,gn,cn): # cell life duration
        gc_xout = self.results['output_g'+str(gn)][str(cn)]['output']['xoutS']
        parp = gc_xout[:,species_all.index('PARP')]
        cparp = gc_xout[:,species_all.index('cPARP')]
        gc_tout = self.results['output_g'+str(gn)][str(cn)]['output']['tout']
        gc_start = gc_tout[0]
        gc_end = gc_tout[-1]
        
        flagA = np.where(cparp>parp)[0]
        if len(flagA)>0:
            t_death = gc_tout[flagA[0]]
            gc_end = t_death
        
        gc_len = gc_end - gc_start              
        return gc_len/3600
    
    def lin_tree_solo(self,g1cx): # single cell lineage tree
        plt.rcParams["lines.linewidth"]=2.0
        desc_all = self.get_desc(g1cx)
        g1_len = self.get_len_gc(1, g1cx)
        newick_str = 'g1c'+str(g1cx)+':'+str(g1_len)
        desc_g1 = self.get_desc_gc(1,g1cx)     
        desc_g1_len = [self.get_len_gc(2,x) for x in desc_g1]        
        newick_insert = '(g2c'+str(desc_g1[0])+':'+str(desc_g1_len[0])+',g2c'+str(desc_g1[1])+':'+str(desc_g1_len[1])+')'
        newick_new = newick_insert+newick_str
        for g in range(2,self.n_gen+1):
            gn_cells = desc_all['g'+str(g)]
            if len(gn_cells)>0:
                for c in gn_cells:
                    mother = 'g'+str(g)+'c'+str(c)
                    desc_gc = self.get_desc_gc(g,c)
                    if len(desc_gc)>0:
                        desc_gc_len = [self.get_len_gc(g+1,x) for x in desc_gc]
                        newick_insert = '(g'+str(g+1)+'c'+str(desc_gc[0])+':'+str(desc_gc_len[0])+',g'+str(g+1)+'c'+str(desc_gc[1])+':'+str(desc_gc_len[1])+')'
                        newick_insert_idx = newick_new.find(mother)
                        newick_update = newick_new[:newick_insert_idx]+newick_insert+newick_new[newick_insert_idx:]
                        newick_new = newick_update
                        
        tree = Phylo.read(StringIO(newick_new),"newick")
        colors = ['blue','orange','green','red','purple','brown','pink','gray','olive','cyan']  
        terminals = tree.get_terminals()
        gxcx = r'g(\d+)c'
        for ter in terminals:
            name = ter.name
            gen = int(re.search(gxcx,name).group(1))
            ter.color = colors[gen-1]
        nonterminals = tree.get_nonterminals()
        for nt in nonterminals:
            name = nt.name
            gen = int(re.search(gxcx,name).group(1))
            nt.color = colors[gen-1]         
                    
        # Phylo.draw(tree,xlabel=("h"),ylabel="c") 
        fig_solo, ax_solo = plt.subplots()
                
        Phylo.draw(tree, axes=ax_solo,do_show=False)
        ax_solo.set_xlabel("Time (hours)",fontsize=15,font='Arial')
        ax_solo.set_ylabel("")
        ax_solo.tick_params(axis='y',labelsize=15,left=False,labelleft=False)
        plt.show()
        
        
    def dendro(self): # cell population dendrogram
        plt.rcParams["lines.linewidth"]=0.5
        cells_n = len(self.results['output_g1'].keys())
        newick_pop = "("
        
        for g1cx in range(1,cells_n+1):
        
            desc_all = self.get_desc(g1cx)
            
            g1_len = self.results['output_g1'][str(g1cx)]['output']['tout'][-1]/3600
            
            newick_new = 'g1c'+str(g1cx)+':'+str(g1_len)
            
            desc_g1 = self.get_desc_gc(1,g1cx)
            
            desc_g1_len = [self.get_len_gc(2,x) for x in desc_g1]
            
            if len(desc_g1)>0:
            
                newick_insert = '(g2c'+str(desc_g1[0])+':'+str(desc_g1_len[0])+',g2c'+str(desc_g1[1])+':'+str(desc_g1_len[1])+')'
                
                newick_new = newick_insert+newick_new
        
        
            
                for g in range(2,self.n_gen+1):
                    gn_cells = desc_all['g'+str(g)]
                    if len(gn_cells)>0:
                        for c in gn_cells:
                            mother = 'g'+str(g)+'c'+str(c)
                            desc_gc = self.get_desc_gc(g,c)
                            if len(desc_gc)>0:
                                desc_gc_len = [self.get_len_gc(g+1,x) for x in desc_gc]
                                newick_insert = '(g'+str(g+1)+'c'+str(desc_gc[0])+':'+str(desc_gc_len[0])+',g'+str(g+1)+'c'+str(desc_gc[1])+':'+str(desc_gc_len[1])+')'
                                newick_insert_idx = newick_new.find(mother)
                                newick_update = newick_new[:newick_insert_idx]+newick_insert+newick_new[newick_insert_idx:]
                                newick_new = newick_update
            if g1cx < cells_n:            
                newick_pop = newick_pop + str(newick_new) + ","
            elif g1cx == cells_n:
                newick_pop = newick_pop + str(newick_new) + ")"
         
        #%
        colors = ['blue','orange','green','red','purple','brown','pink','gray','olive','cyan'] 
        tree_pop = Phylo.read(StringIO(newick_pop),"newick")
        terminals = tree_pop.get_terminals()
        gxcx = r'g(\d+)c'
        for ter in terminals:
            name = ter.name
            gen = int(re.search(gxcx,name).group(1))
            ter.color = colors[gen-1]
            ter.name = ""
        nonterminals = tree_pop.get_nonterminals()[1:]
        for nt in nonterminals:
            name = nt.name
            gen = int(re.search(gxcx,name).group(1))
            nt.color = colors[gen-1]
            nt.name = ""
        
        fig_pop, ax_pop = plt.subplots()
        
        
        Phylo.draw(tree_pop, axes=ax_pop,do_show=False)
        # ax_pop.set_title("Population dendrogram",fontsize=10)
        ax_pop.set_xlabel("Time (hours)",fontsize=15,font='Arial')
        ax_pop.set_ylabel("")
        ax_pop.set_xlim(0,72)
        ax_pop.tick_params(axis='y',labelsize=15,left=False,labelleft=False)
        
        
        plt.show()       
        

    
#%%
exp_title = 'in_silico_drs'
output_main = os.path.join(wd,'output',exp_title)

output_noegf = os.path.join(wd,'output','mcf10a_noegf')
output_nostim = os.path.join(wd,'output','mcf10a_nostim')

dict_noegf = drs_dict(output_noegf,'nerat',1,0)
dict_nostim = drs_dict(output_nostim,'nerat',1,0)

exp_drs2 = 'in_silico_drs2'
output_drs2 = os.path.join(wd,'output',exp_drs2)

exp_drs3 = 'in_silico_drs3'
output_drs3 = os.path.join(wd,'output',exp_drs3)

#%%
nerat1_5 = drs_dict(output_main,'nerat',1,5)
aaa1=nerat1_5.pop_dyn_obs(formula_egfr)[1]
aaa3=nerat1_5.pop_dyn_obs(formula_egfr)[3]

aaa56 = nerat1_5.pop_dyn_obs(formula_ppERK)[3]

obs_median_pop = np.array([np.nanmedian(aaa3[:,tp]) for tp in range(len(aaa1))])

obs_median_pop2 = np.array([np.nanmedian(aaa56[:,tp]) for tp in range(len(aaa1))])

plt.plot(aaa1/3600,obs_median_pop)

plt.plot(aaa1/3600,obs_median_pop2)

#%%

def pop_obs_plot(obs_array,timepoints):
    obs_mean = np.array([np.nanmean(obs_array[:,tp]) for tp in range(len(timepoints))])
    obs_sd = np.array([np.nanstd(obs_array[:,tp]) for tp in range(len(timepoints))])


    for cell in range(np.shape(obs_array)[0]):
        plt.plot(timepoints/3600,obs_array[cell,:],linewidth=0.5,color='grey')
    
    plt.plot(timepoints/3600,obs_mean,linewidth=2.0,color='red',label='Mean')
    plt.plot(timepoints/3600,obs_mean+obs_sd,linewidth=2.0,color='blue',label='+SD')
    plt.plot(timepoints/3600,obs_mean-obs_sd,linewidth=2.0,color='blue',label='-SD')
    
    
    ymax = np.nanmax(obs_array)*1.25
    
    plt.ylim(0,ymax)
    plt.xlim(0,74)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(bbox_to_anchor=(1.2,1.2),fontsize=20)
    
    plt.show()
    
#%% dose specific target engagement

for dl in range(len(doses_all)):
    dose_dict = drs_dict(output_main,'nerat',1,dl)
    cells,tps,tout_deaths,obs_array = dose_dict.pop_dyn_obs(formula_egfr)
    obs_median = np.array([np.nanmedian(obs_array[:,tp]) for tp in range(len(tps))])
    plt.plot(tps/3600,obs_median,label=doses_all[dl])
    
plt.ylim(0,1.2)
plt.xlim(0,74)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(bbox_to_anchor=(1.05,1.1),title='Doses (uM)',fontsize=15,title_fontsize=15)
plt.show()        

#%% dose specific erk activity

for dl in range(len(doses_all)):
    dose_dict = drs_dict(output_main,'nerat',1,dl)
    cells,tps,tout_deaths,obs_array = dose_dict.pop_dyn_obs(formula_ppERK)
    pop_obs_plot(obs_array,tps)

#%% palbociclib target engagement

for dl in range(len(doses_all)):
    dose_dict = drs_dict(output_main,'palbo',1,dl)
    cells,tps,tout_deaths,obs_array = dose_dict.pop_dyn_obs(formula_palbo_target)
    obs_median = np.array([np.nanmedian(obs_array[:,tp]) for tp in range(len(tps))])
    plt.plot(tps/3600,obs_median,label=doses_all[dl])
    
# plt.ylim(0,1.2)
plt.xlim(0,74)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(bbox_to_anchor=(1.05,1.1),title='Doses (uM)',fontsize=15,title_fontsize=15)
plt.show()   

#%%

obs_mean = np.array([np.nanmean(aaa56[:,tp]) for tp in range(len(aaa1))])
obs_sd = np.array([np.nanstd(aaa56[:,tp]) for tp in range(len(aaa1))])


for cell in range(np.shape(aaa56)[0]):
    plt.plot(aaa1/3600,aaa56[cell,:],linewidth=0.5,color='grey')

plt.plot(aaa1/3600,obs_mean,linewidth=2.0,color='red')
plt.plot(aaa1/3600,obs_mean+obs_sd,linewidth=2.0,color='blue')
plt.plot(aaa1/3600,obs_mean-obs_sd,linewidth=2.0,color='blue')


ymax = np.nanmax(aaa56)*1.25

plt.ylim(0,ymax)
plt.xlim(0,74)

plt.show()



#%% drs_PopDyn

# drug = 'alpel'

# drs_drug = {}

def drs_summarize(drug,dose_lvl,output_dir=output_main):
    
    drs_dose = {}
    for rep in range(10):
        print('now running...'+str(drug)+'...'+str(dose_lvl)+'...'+str(rep+1))
        drs_dict0 = drs_dict(output_dir,drug,rep+1,dose_lvl)
        pd,tp,td = drs_dict0.pop_dyn()
        drs_rep = {}
        drs_rep['cellpop'] = pd
        drs_rep['tout'] = tp
        drs_rep['t_death'] = td
        drs_dose['r'+str(rep+1)] = drs_rep
    return drs_dose

#%%

drs_all = {}

#%%
drug = 'trame'

#%%

drs_drug = {}

for d in range(10):
    
    drs_dose = drs_summarize(drug, d)
    drs_drug['d'+str(d)] = drs_dose
    
drs_all[drug] = drs_drug

#%%

pickle.dump(drs_all, open(os.path.join(wd,'output','in_silico_drs_summary','drs_summary.pkl'),'wb'))

#%% drs/rep-median


drs_median = {}

drugs = list(drs_all.keys())

for dr in drugs:
    
    drs_median_drug = {}
    
    for dl in range(10):
        
        drs_median_dose = {}
        
        tp_drs = [drs_all[dr]['d'+str(dl)]['r'+str(rep+1)]['tout'] for rep in range(10)]
        popdyn_reps0 = [drs_all[dr]['d'+str(dl)]['r'+str(rep+1)]['cellpop'] for rep in range(10)]
        tp_all = np.array(list(itertools.chain(*tp_drs)))
        tp_all = np.unique(tp_all)
        tp_max = min([tp_drs[x][-1] for x in range(len(tp_drs))])
        tp_max_idx = np.where(tp_all == tp_max)[0][0]
        tp_all = tp_all[:tp_max_idx+1]
        popdyn_reps = []
        for rep in range(10):
            interpolator = interp1d(tp_drs[rep],popdyn_reps0[rep])
            y_new = interpolator(tp_all)
            popdyn_reps.append(y_new)
            
        popdyn_med = np.median(popdyn_reps,axis=0)
        
        drs_median_dose['cellpop'] = popdyn_med
        drs_median_dose['tout'] = tp_all
        
        drs_median_drug['d'+str(dl)] = drs_median_dose
        
    drs_median[dr] = drs_median_drug

#%%

pickle.dump(drs_median, open(os.path.join(wd,'output','in_silico_drs_summary','drs_median.pkl'),'wb'))

#%%

with open(os.path.join(wd,'output','in_silico_drs_summary','drs_summary.pkl'),'rb') as f:
    drs_summary_full = pickle.load(f)

with open(os.path.join(wd,'output','in_silico_drs_summary','drs_median.pkl'),'rb') as f:
    drs_summary_median = pickle.load(f)

#%%
drs_median = drs_summary_median
drugs = list(drs_median.keys())
time_n = 72.0

drs_grcalc = pd.DataFrame(data=np.ones((45,6))*np.nan,columns=['cell_line','treatment','concentration','cell_count','cell_count__ctrl','cell_count__time0'])

drs_grcalc.loc[:,'cell_line'] = 'mcf10a'

k = 0
for dr,drug in enumerate(drugs):
    drs_grcalc.loc[:,'treatment'][k:k+9] = str(drug)
    drs_grcalc.loc[:,'concentration'][k:k+9] = doses_all[1:]
    
    drs_ctrl = drs_median[drug]['d0']
    
    tn_idx0 = np.where(drs_ctrl['tout']/3600>time_n)[0][0]
    
    x_ctrl = drs_ctrl['cellpop'][tn_idx0]
    
    
    drs_grcalc.loc[:,'cell_count__ctrl'][k:k+9] = x_ctrl
    
    x_t0 = []
    x_tn = []
    
    for dl in range(1,10):
        
        drs_dl = drs_median[drug]['d'+str(dl)]
        tn_idx = np.where(drs_dl['tout']/3600>time_n)[0][0]
        
        x1 = drs_dl['cellpop'][0]
        x2 = drs_dl['cellpop'][tn_idx]
        
        x_t0.append(x1)
        x_tn.append(x2)
        
    drs_grcalc.loc[:,'cell_count'][k:k+9] = x_tn
    drs_grcalc.loc[:,'cell_count__time0'][k:k+9] = x_t0
    # for dl in range(1,10):
    #     drs_grcalc[drs_grcalc['treatment']==str(drug)]['concentration'][int(dl-1)] = doses_all[dl]
    k = k+9

#%%

drs_grcalc.to_csv(os.path.join(wd,'output','in_silico_drs_summary','drs_grcalc.tsv'),sep='\t',index=False)


#%% drs_plots

drug = drugs[4]

for dl in range(10):
    
    dose = doses_all[dl]
    
    x_dl = drs_median[drug]['d'+str(dl)]['tout']/3600
    y_dl = drs_median[drug]['d'+str(dl)]['cellpop']
    
    plt.plot(x_dl,y_dl,label=str(dose))
    

plt.legend(bbox_to_anchor=(1.05,1.0))
plt.ylim(0,550)
plt.xlim(0,72)
plt.xlabel('Time (h)')
plt.ylabel('# of cells')
plt.title(drug)
plt.show()

#%%

pickle.dump(drs_median, open(os.path.join(wd,'output','in_silico_drs_summary','drs_median.pkl'),'wb'))

#%%

with open(os.path.join(wd,'output','in_silico_drs_summary','drs_summary.pkl'),'rb') as f:
    drs_summary_full = pickle.load(f)

#%%

def gr_calc_row (drug,time_h,dl,rep,drs_summary_dict,cell_line='mcf10a_sim'):
    
    dose = doses_all[dl]
    
    dose_dict = drs_summary_dict[str(drug)]['d'+str(dl)]['r'+str(rep+1)]
    
    cellpop = dose_dict['cellpop']
    tout = dose_dict['tout']
    
    interpolator = interp1d(tout,cellpop)
    
    cell_count = interpolator(time_h*3600)
    
    new_row = {}
    new_row['cell_line'] = cell_line
    new_row['treatment'] = str(drug)
    new_row['treatment_duration__hrs'] = time_h
    new_row['concentration'] = dose
    new_row['cell_count'] = cell_count
    
    
    return new_row





#%%


time_h = [48,72]

drs_grcalc2 = pd.DataFrame(data=None,columns=['cell_line','treatment','treatment_duration__hrs','concentration','cell_count'])



for drug in drs_summary_full.keys():
    for dl in range(len(drs_summary_full[str(drug)].keys())):
        for rep in range(len(drs_summary_full[str(drug)]['d'+str(dl)].keys())):
            new_row0 = gr_calc_row(drug,0,dl,rep,drs_summary_full)
            drs_grcalc2 = drs_grcalc2.append(new_row0, ignore_index=True)
            for t in time_h:
                new_row = gr_calc_row(drug,t,dl,rep,drs_summary_full)
                drs_grcalc2 = drs_grcalc2.append(new_row, ignore_index=True)
            
    
drs_grcalc2.to_csv(os.path.join(wd,'output','in_silico_drs_summary','drs_grcalc2.tsv'),sep='\t',index=False)


#%%


def gr_calc_row3 (drug,time_h,dl,rep,drs_summary_dict,cell_line='mcf10a_sim'):
    
    dose = doses_all[dl]
    
    dose_dict = drs_summary_dict[str(drug)]['d'+str(dl)]['r'+str(rep+1)]
    # ctrl_dict = drs_summary_dict[str(drug)]['d0']['r'+str(rep+1)]
    
    cellpop = dose_dict['cellpop']
    tout = dose_dict['tout']
    
    # cellpop_ctrl = ctrl_dict['cellpop']
    # tout_ctrl = ctrl_dict['tout']
    
    interpolator = interp1d(tout,cellpop)
    # interpolator_ctrl = interp1d(tout_ctrl,cellpop_ctrl)
    
    cell_count = interpolator(time_h*3600)
    
    # cell_count_ctrl = interpolator_ctrl(time_h*3600)
    
    cell_count_time0 = dose_dict['cellpop'][0]
    
    drugs_all = list(drs_summary_dict.keys())
 

    drugs_all = list(drs_summary_dict.keys())
    
    ctrl_pool = []
    
    for dr in drugs_all:
        for rp in range(1,len(drs_summary_dict[dr]['d0'].keys())+1):
            rp_dict = drs_summary_dict[dr]['d0']['r'+str(rp)]
            rp_cellpop = np.array(rp_dict['cellpop'])
            rp_tout = np.array(rp_dict['tout'])
            
            rp_interp = interp1d(rp_tout,rp_cellpop)
            rp_cellcount = float(rp_interp(time_h*3600))
            ctrl_pool.append(rp_cellcount)
            
    dose_pool = []
    
    for rp in range(1,len(drs_summary_dict[drug]['d'+str(dl)].keys())+1):
        rp_dose_dict = drs_summary_dict[drug]['d'+str(dl)]['r'+str(rp)]
        rp_dose_cellpop = np.array(rp_dose_dict['cellpop'])
        rp_dose_tout = np.array(rp_dose_dict['tout'])
        
        rp_dose_interp = interp1d(rp_dose_tout,rp_dose_cellpop)
        rp_dose_cellcount = float(rp_dose_interp(time_h*3600))
        dose_pool.append(rp_dose_cellcount)    
 
    cell_count_percentile = percentileofscore(dose_pool, cell_count)
    cell_count_ctrl = np.percentile(ctrl_pool,cell_count_percentile)


    new_row = {}
    new_row['cell_line'] = cell_line
    new_row['agent'] = str(drug)
    new_row['timepoint'] = time_h
    new_row['concentration'] = dose
    new_row['cell_count'] = cell_count
    new_row['cell_count__ctrl'] = cell_count_ctrl
    new_row['cell_count__time0'] = cell_count_time0 
    
    
    return new_row

#%%
time_h = 72
dl = 5
rep = 2
drug = 'trame'
drs_summary_dict = drs_summary_full


drugs_all = list(drs_summary_dict.keys())

ctrl_pool = []

for dr in drugs_all:
    for rp in range(1,len(drs_summary_dict[dr]['d0'].keys())+1):
        rp_dict = drs_summary_dict[dr]['d0']['r'+str(rp)]
        rp_cellpop = np.array(rp_dict['cellpop'])
        rp_tout = np.array(rp_dict['tout'])
        
        rp_interp = interp1d(rp_tout,rp_cellpop)
        rp_cellcount = float(rp_interp(time_h*3600))
        ctrl_pool.append(rp_cellcount)
        
dose_pool = []

for rp in range(1,len(drs_summary_dict[drug]['d'+str(dl)].keys())+1):
    rp_dose_dict = drs_summary_dict[drug]['d'+str(dl)]['r'+str(rp)]
    rp_dose_cellpop = np.array(rp_dose_dict['cellpop'])
    rp_dose_tout = np.array(rp_dose_dict['tout'])
    
    rp_dose_interp = interp1d(rp_dose_tout,rp_dose_cellpop)
    rp_dose_cellcount = float(rp_dose_interp(time_h*3600))
    dose_pool.append(rp_dose_cellcount)


dose_dict = drs_summary_dict[str(drug)]['d'+str(dl)]['r'+str(rep+1)]

cellpop = dose_dict['cellpop']
tout = dose_dict['tout']
interpolator = interp1d(tout,cellpop)

cell_count = float(interpolator(time_h*3600))
cell_count_time0 = dose_dict['cellpop'][0]

cell_count_percentile = percentileofscore(dose_pool, cell_count)
cell_count_ctrl = np.percentile(ctrl_pool,cell_count_percentile)


#%%

plt.hist(ctrl_pool)
plt.show()

#%%



percentileofscore(ctrl_pool, ctrl_pool[41])

#%%

time_h = [48,72]

drs_grcalc3 = pd.DataFrame(data=None,columns=['cell_line','agent','timepoint','concentration','cell_count','cell_count__ctrl','cell_count__time0'])



for drug in drs_summary_full.keys():
    for dl in range(1,len(drs_summary_full[str(drug)].keys())):
        for rep in range(len(drs_summary_full[str(drug)]['d'+str(dl)].keys())):
            for t in time_h:
                new_row = gr_calc_row3(drug,t,dl,rep,drs_summary_full)
                drs_grcalc3 = drs_grcalc3.append(new_row, ignore_index=True)


drs_grcalc3.to_csv(os.path.join(wd,'output','in_silico_drs_summary','drs_grcalc3.tsv'),sep='\t',index=False)


#%%

drs_all_gr = pd.read_csv(os.path.join(wd,'output','in_silico_drs_summary','drs_grcalc3_grc.tsv'),sep='\t')



#%% load exp data

# drs_gr_exp = pd.read_csv(os.path.join(wd,'output','in_silico_drs_summary','grvalues_merged.csv'),sep=',')

dir_exp = os.path.join(wd,'output','in_silico_drs_summary','mcf10a_drs_exp')

drs_gr_exp = {}

grv_c1a1_data = pd.read_csv(os.path.join(dir_exp,'GRvalues_center1_scientistA_2017.csv'),sep=',')

grv_c1a1 = pd.DataFrame()

grv_c1a1['cell_line'] = grv_c1a1_data['cell_line']
grv_c1a1['agent'] = grv_c1a1_data['agent']
grv_c1a1['concentration'] = grv_c1a1_data['concentration']
grv_c1a1['cell_count'] = grv_c1a1_data['cell_count']
grv_c1a1['cell_count__ctrl'] = grv_c1a1_data['cell_count__ctrl']
grv_c1a1['cell_count__time0'] = grv_c1a1_data['cell_count__time0']
grv_c1a1['GRvalue'] = grv_c1a1_data['GRvalue']


#
grv_c1a2_data = pd.read_csv(os.path.join(dir_exp,'GRvalues_center1_scientistA_2019.csv'),sep=',')
grv_c1a2 = pd.DataFrame()
grv_c1a2['cell_line'] = grv_c1a2_data['cell_line']
grv_c1a2['agent'] = grv_c1a2_data['agent']
grv_c1a2['timepoint'] = grv_c1a2_data['timepoint']
grv_c1a2['cell_count'] = grv_c1a2_data['cell_count']
grv_c1a2['cell_count__ctrl'] = grv_c1a2_data['cell_count__ctrl']
grv_c1a2['cell_count__time0'] = grv_c1a2_data['cell_count__time0']
grv_c1a2['GRvalue'] = grv_c1a2_data['GRvalue']

#
#%%

def read_exp_data(file_name):
    dir_exp = os.path.join(wd,'output','in_silico_drs_summary','mcf10a_drs_exp')
    exp_data = pd.read_csv(os.path.join(dir_exp,file_name),sep=',')
    grv_df = pd.DataFrame()
    grv_df['cell_line'] = exp_data['cell_line']
    
    drugs = exp_data['agent'].values
    
    for d in range(len(drugs)):
        if '/' in drugs[d]:
            drug = str(drugs[d]).split('/')[0]
            drugs[d] = drug
    
    # grv_df['agent'] = exp_data['agent']
    grv_df['agent'] = drugs
    grv_df['concentration'] = exp_data['concentration']
    try:
        grv_df['timepoint'] = exp_data['timepoint']
    except:
        grv_df['timepoint'] = np.ones(len(exp_data['cell_line'].values))*np.nan
    grv_df['cell_count'] = exp_data['cell_count']
    grv_df['cell_count__ctrl'] = exp_data['cell_count__ctrl']
    grv_df['cell_count__time0'] = exp_data['cell_count__time0']
    grv_df['GRvalue'] = exp_data['GRvalue']
    
    return grv_df

#%%
grv_c1a1 = read_exp_data('GRvalues_center1_scientistA_2017.csv')
grv_c1a2 = read_exp_data('GRvalues_center1_scientistA_2019.csv')
grv_c1b = read_exp_data('GRvalues_center1_scientistB.csv')
grv_c1c = read_exp_data('GRvalues_center1_scientistC.csv')
grv_c2 = read_exp_data('GRvalues_center2.csv')
grv_c3 = read_exp_data('GRvalues_center3.csv')
grv_c4 = read_exp_data('GRvalues_center4.csv')
grv_c5 = read_exp_data('GRvalues_center5.csv')

#%% plot exp data
drugs_exp = ['Alpelisib','Neratinib','Trametinib','Palbociclib']
grv_df = grv_c1a1
drug = drugs_exp[3]

grv_drug = grv_df[grv_df['agent']==drug]

doses_drug = np.unique(grv_drug['concentration'].values)

# dose_lvl_count = [sum(grv_drug['concentration'].values==d) for d in dose_lvls]
x_values = doses_drug

y_values_all = [grv_drug[grv_drug['concentration']==xval]['GRvalue'].values for xval in x_values]

y_values_median = [np.median(ys) for ys in y_values_all]

y_min = [np.min(ys) for ys in y_values_all]

y_max = [np.max(ys) for ys in y_values_all]

y_err_min = [y_values_median[dl] - y_min[dl] for dl in range(len(y_values_median))]

y_err_max = [y_max[dl] - y_values_median[dl] for dl in range(len(y_values_median))]

yerror = [y_err_min,y_err_max]

plt.errorbar(x_values,y_values_median,yerr=yerror,fmt='o-',capsize=5)
plt.xscale('log')
# plt.ylim(min(y_min)*.5,max(y_max)*1.25)
plt.ylim(0,max(y_max)*1.25)
plt.title(str(drug))
plt.ylabel('GR value')
plt.xlabel('Concentration (uM)')
plt.show()

#%% combine and map exp doses to sim

drugs_exp = ['Alpelisib','Neratinib','Trametinib','Palbociclib']
grv_exp = {'grv_c1a1':grv_c1a1,'grv_c1a2':grv_c1a2,'grv_c1b':grv_c1b,'grv_c1c':grv_c1c,'grv_c2':grv_c2,'grv_c3':grv_c3,'grv_c4':grv_c4,'grv_c5':grv_c5}


grv_exp_df = pd.DataFrame()

for key in grv_exp.keys():
    grv_df = grv_exp[key]
    grv_df['source'] = np.ones(np.shape(grv_df)[0])*np.nan
    grv_df['source'] = key
    grv_exp_df = grv_exp_df.append(grv_df, ignore_index=True)
    

grv_exp_compare = pd.DataFrame()

# grv_c2[math.isclose(grv_c2['concentration'],doses_all[3],abs_tol=1e-4)]

for dose_sim in doses_all[1:]:
    for row_n in range(np.shape(grv_exp_df)[0]):
        grv_exp_row = grv_exp_df.iloc[row_n,:].copy()
        if grv_exp_row['agent'] in drugs_exp:
            if math.isclose(float(grv_exp_row['concentration']),dose_sim,abs_tol=1e-4):
                grv_exp_row['dose_sim'] = dose_sim
                grv_exp_compare = grv_exp_compare.append(grv_exp_row,ignore_index=True)

#%% gr compare

drug_compare = drugs_exp[3]

drug_sim = drug_compare.lower()[:5]

grv_exp_drug = grv_exp_compare[grv_exp_compare['agent']==drug_compare]

x_values = doses_all[1:]

y_values_all_exp = [grv_exp_drug[grv_exp_drug['dose_sim']==xval]['GRvalue'].values for xval in x_values]

y_values_median_exp = [np.median(ys) for ys in y_values_all_exp]

if sum(np.isnan(y_values_median_exp))>0:
    
    nan_idx = np.where(np.isnan(y_values_median_exp))[0][0]
    
    x_values_actual = x_values[:nan_idx]
    x_values = x_values_actual
    
    y_values_all_exp = [grv_exp_drug[grv_exp_drug['dose_sim']==xval]['GRvalue'].values for xval in x_values]

    y_values_median_exp = [np.median(ys) for ys in y_values_all_exp]
    
    


y_min_exp = [np.min(ys) for ys in y_values_all_exp]

y_max_exp = [np.max(ys) for ys in y_values_all_exp]

y_err_min_exp = [y_values_median_exp[dl] - y_min_exp[dl] for dl in range(len(y_values_median_exp))]

y_err_max_exp = [y_max_exp[dl] - y_values_median_exp[dl] for dl in range(len(y_values_median_exp))]

yerror_exp = [y_err_min_exp,y_err_max_exp]


#%

tp = 72

drs_filtered = drs_all_gr[drs_all_gr['agent']==drug_sim]
drs_filtered = drs_filtered[drs_filtered['timepoint']==tp]

y_values_all_sim = [drs_filtered[drs_filtered['concentration']==xval]['GRvalue'].values for xval in x_values]

y_values_median_sim = [np.median(ys) for ys in y_values_all_sim]

y_min_sim = [np.min(ys) for ys in y_values_all_sim]

y_max_sim = [np.max(ys) for ys in y_values_all_sim]

y_err_min_sim = [y_values_median_sim[dl] - y_min_sim[dl] for dl in range(len(y_values_median_sim))]

y_err_max_sim = [y_max_sim[dl] - y_values_median_sim[dl] for dl in range(len(y_values_median_sim))]

yerror_sim = [y_err_min_sim,y_err_max_sim]



#%

plt.errorbar(x_values,y_values_median_exp,yerr=yerror_exp,fmt='o-',capsize=5,label='Experiment')
plt.errorbar(x_values,y_values_median_sim,yerr=yerror_sim,fmt='o-',capsize=5,label='Simulation')

plt.xscale('log')
# plt.ylim(min(y_min)*.5,max(y_max)*1.25)
plt.ylim(-0.5,1.3)
plt.title(str(drug_compare),font='Arial')
plt.ylabel('GR value',font='Arial')
plt.xlabel('Concentration (uM)',font='Arial')
plt.legend(loc='lower left')
plt.show()



#%%

drugs_exp = np.unique(drs_gr_exp['agent'].values)

drug_x = drugs_exp[0]

drs_exp_filtered = drs_gr_exp[drs_gr_exp['agent']==drug_x]

x_values = np.unique(drs_gr_exp['concentration'].values)
y_values_all = []


#%%

plt.plot(x_values,y_values_median)
plt.show()



#%% drs2 test

nerat1_0 = drs_dict(output_drs2,'nerat',1,0)

nerat1_9 = drs_dict(output_drs2,'nerat',1,9)

#%% drs2 plot

pd1, tp1, td1 = nerat1_0.pop_dyn()
pd9, tp9, td9 = nerat1_9.pop_dyn()

plt.plot(tp1/3600,pd1)
plt.plot(tp9/3600,pd9)

#%%
rep = 9
drug = 'nerat'

for dl in range(10):
    drs_dict0 = drs_dict(output_drs2,drug,rep,dl)
    pd,tp,td = drs_dict0.pop_dyn()
    plt.plot(tp/3600,pd,label=str(doses_all[dl]))
    
plt.ylim(0,200)
plt.xlim(0,72)
plt.xlabel('Time(h)')
plt.ylabel('Cell #')
plt.legend(bbox_to_anchor=(1.05,1.0))
plt.show()

#%% summarize drs2
drs2_all = {}

#%%

drug = 'trame'

#%%

drs_drug = {}

for d in range(10):
    
    drs_dose = drs_summarize(drug, d,output_drs2)
    drs_drug['d'+str(d)] = drs_dose
    
drs2_all[drug] = drs_drug

#%%

pickle.dump(drs2_all, open(os.path.join(wd,'output','in_silico_drs2_summary','drs_summary.pkl'),'wb'))

#%% gr-calc/drs2/nerat

drs_summary= drs2_all

def gr_calc_row3 (drug,time_h,dl,rep,drs_summary_dict,cell_line='mcf10a_sim'):
    
    dose = doses_all[dl]
    
    dose_dict = drs_summary_dict[str(drug)]['d'+str(dl)]['r'+str(rep+1)]
    # ctrl_dict = drs_summary_dict[str(drug)]['d0']['r'+str(rep+1)]
    
    cellpop = dose_dict['cellpop']
    tout = dose_dict['tout']
    
    # cellpop_ctrl = ctrl_dict['cellpop']
    # tout_ctrl = ctrl_dict['tout']
    
    interpolator = interp1d(tout,cellpop)
    # interpolator_ctrl = interp1d(tout_ctrl,cellpop_ctrl)
    
    cell_count = interpolator(time_h*3600)
    
    # cell_count_ctrl = interpolator_ctrl(time_h*3600)
    
    cell_count_time0 = dose_dict['cellpop'][0]
    
    # drugs_all = list(drs_summary_dict.keys())
 

    # drugs_all = list(drs_summary_dict.keys())
    drugs_all = ['nerat']
    
    ctrl_pool = []
    
    for dr in drugs_all:
        for rp in range(1,len(drs_summary_dict[dr]['d0'].keys())+1):
            rp_dict = drs_summary_dict[dr]['d0']['r'+str(rp)]
            rp_cellpop = np.array(rp_dict['cellpop'])
            rp_tout = np.array(rp_dict['tout'])
            
            rp_interp = interp1d(rp_tout,rp_cellpop)
            rp_cellcount = float(rp_interp(time_h*3600))
            ctrl_pool.append(rp_cellcount)
            
    dose_pool = []
    
    for rp in range(1,len(drs_summary_dict[drug]['d'+str(dl)].keys())+1):
        rp_dose_dict = drs_summary_dict[drug]['d'+str(dl)]['r'+str(rp)]
        rp_dose_cellpop = np.array(rp_dose_dict['cellpop'])
        rp_dose_tout = np.array(rp_dose_dict['tout'])
        
        rp_dose_interp = interp1d(rp_dose_tout,rp_dose_cellpop)
        rp_dose_cellcount = float(rp_dose_interp(time_h*3600))
        dose_pool.append(rp_dose_cellcount)    
 
    cell_count_percentile = percentileofscore(dose_pool, cell_count)
    cell_count_ctrl = np.percentile(ctrl_pool,cell_count_percentile)


    new_row = {}
    new_row['cell_line'] = cell_line
    new_row['agent'] = str(drug)
    new_row['timepoint'] = time_h
    new_row['concentration'] = dose
    new_row['cell_count'] = cell_count
    new_row['cell_count__ctrl'] = cell_count_ctrl
    new_row['cell_count__time0'] = cell_count_time0 
    
    
    return new_row

#%%



time_h = [48,72]
drs_summary_full = drs_summary

drs2_grcalc3 = pd.DataFrame(data=None,columns=['cell_line','agent','timepoint','concentration','cell_count','cell_count__ctrl','cell_count__time0'])



# drug = 'nerat'
# for dl in range(1,len(drs_summary_full[str(drug)].keys())):
#     for rep in range(len(drs_summary_full[str(drug)]['d'+str(dl)].keys())):
#         for t in time_h:
#             new_row = gr_calc_row3(drug,t,dl,rep,drs_summary_full)
#             drs_grcalc3_nerat = drs_grcalc3_nerat.append(new_row, ignore_index=True)

for drug in drs_summary_full.keys():
    for dl in range(1,len(drs_summary_full[str(drug)].keys())):
        for rep in range(len(drs_summary_full[str(drug)]['d'+str(dl)].keys())):
            for t in time_h:
                new_row = gr_calc_row3(drug,t,dl,rep,drs_summary_full)
                drs2_grcalc3 = drs2_grcalc3.append(new_row, ignore_index=True)


drs2_grcalc3.to_csv(os.path.join(wd,'output','in_silico_drs2_summary','drs2_grcalc3.tsv'),sep='\t',index=False)

#%% drs2/nerat/gr-plot

drs2_nerat_gr = pd.read_csv(os.path.join(wd,'output','in_silico_drs2_summary','drs2_grcalc3_nerat_grc.tsv'),sep='\t')

drug = 'nerat'
tp = 72

# Sample data (replace this with your actual data)
x_values = doses_all[1:]

drs_filtered = drs2_nerat_gr[drs2_nerat_gr['agent']==drug]
drs_filtered = drs_filtered[drs_filtered['timepoint']==tp]

y_values_all = [drs_filtered[drs_filtered['concentration']==xval]['GRvalue'].values for xval in x_values]

y_values_median = [np.median(ys) for ys in y_values_all]

y_min = [np.min(ys) for ys in y_values_all]

y_max = [np.max(ys) for ys in y_values_all]

y_err_min = [y_values_median[dl] - y_min[dl] for dl in range(len(y_values_median))]

y_err_max = [y_max[dl] - y_values_median[dl] for dl in range(len(y_values_median))]

yerror = [y_err_min,y_err_max]

plt.errorbar(x_values,y_values_median,yerr=yerror,fmt='o-',capsize=5)
plt.xscale('log')
# plt.ylim(min(y_min)*.5,max(y_max)*1.25)
plt.ylim(min(y_min)*1.1,max(y_max)*1.25)
plt.title(str(drug))
plt.ylabel('GR value')
plt.xlabel('Concentration (uM)')
plt.show()

#%% drs2 gr plot vs exp

drugs_exp = ['Alpelisib','Neratinib','Trametinib','Palbociclib']
grv_exp = {'grv_c1a1':grv_c1a1,'grv_c1a2':grv_c1a2,'grv_c1b':grv_c1b,'grv_c1c':grv_c1c,'grv_c2':grv_c2,'grv_c3':grv_c3,'grv_c4':grv_c4,'grv_c5':grv_c5}


grv_exp_df = pd.DataFrame()

for key in grv_exp.keys():
    grv_df = grv_exp[key]
    grv_df['source'] = np.ones(np.shape(grv_df)[0])*np.nan
    grv_df['source'] = key
    grv_exp_df = grv_exp_df.append(grv_df, ignore_index=True)
    

grv_exp_compare = pd.DataFrame()

# grv_c2[math.isclose(grv_c2['concentration'],doses_all[3],abs_tol=1e-4)]

for dose_sim in doses_all[1:]:
    for row_n in range(np.shape(grv_exp_df)[0]):
        grv_exp_row = grv_exp_df.iloc[row_n,:].copy()
        if grv_exp_row['agent'] in drugs_exp:
            if math.isclose(float(grv_exp_row['concentration']),dose_sim,abs_tol=1e-4):
                grv_exp_row['dose_sim'] = dose_sim
                grv_exp_compare = grv_exp_compare.append(grv_exp_row,ignore_index=True)


#%% gr compare/drs2
drs_all_gr = pd.read_csv(os.path.join(wd,'output','in_silico_drs2_summary','drs2_grcalc3_grc.tsv'),sep='\t')

plt.rcParams["lines.linewidth"] = 2

drug_compare = drugs_exp[3]

drug_sim = drug_compare.lower()[:5]

grv_exp_drug = grv_exp_compare[grv_exp_compare['agent']==drug_compare]

x_values = doses_all[1:]

y_values_all_exp = [grv_exp_drug[grv_exp_drug['dose_sim']==xval]['GRvalue'].values for xval in x_values]

y_values_median_exp = [np.median(ys) for ys in y_values_all_exp]

if sum(np.isnan(y_values_median_exp))>0:
    
    nan_idx = np.where(np.isnan(y_values_median_exp))[0][0]
    
    x_values_actual = x_values[:nan_idx]
    x_values = x_values_actual
    
    y_values_all_exp = [grv_exp_drug[grv_exp_drug['dose_sim']==xval]['GRvalue'].values for xval in x_values]

    y_values_median_exp = [np.median(ys) for ys in y_values_all_exp]
    
    


y_min_exp = [np.min(ys) for ys in y_values_all_exp]

y_max_exp = [np.max(ys) for ys in y_values_all_exp]

y_err_min_exp = [y_values_median_exp[dl] - y_min_exp[dl] for dl in range(len(y_values_median_exp))]

y_err_max_exp = [y_max_exp[dl] - y_values_median_exp[dl] for dl in range(len(y_values_median_exp))]

yerror_exp = [y_err_min_exp,y_err_max_exp]


#%

tp = 72

drs_filtered = drs_all_gr[drs_all_gr['agent']==drug_sim]
drs_filtered = drs_filtered[drs_filtered['timepoint']==tp]

y_values_all_sim = [drs_filtered[drs_filtered['concentration']==xval]['GRvalue'].values for xval in x_values]

y_values_median_sim = [np.median(ys) for ys in y_values_all_sim]

y_min_sim = [np.min(ys) for ys in y_values_all_sim]

y_max_sim = [np.max(ys) for ys in y_values_all_sim]

y_err_min_sim = [y_values_median_sim[dl] - y_min_sim[dl] for dl in range(len(y_values_median_sim))]

y_err_max_sim = [y_max_sim[dl] - y_values_median_sim[dl] for dl in range(len(y_values_median_sim))]

yerror_sim = [y_err_min_sim,y_err_max_sim]

y_min = min(min(y_min_exp),min(y_min_sim))
y_max = max(max(y_max_exp),max(y_max_sim))

#%

plt.errorbar(x_values,y_values_median_exp,yerr=yerror_exp,fmt='o-',capsize=5,label='Experiment')
plt.errorbar(x_values,y_values_median_sim,yerr=yerror_sim,fmt='o-',capsize=5,label='Simulation')

plt.xscale('log')
plt.ylim(y_min*1.25,y_max*1.25)
# plt.ylim(-1,1.8)
plt.title(str(drug_compare),font='Arial')
plt.ylabel('GR value',font='Arial')
plt.xlabel('Concentration (uM)',font='Arial')
plt.legend(loc='lower left')
plt.show()

#%% drs2/rep-median

drs_all = drs2_all


drs_median = {}

# drugs = list(drs_all.keys())

drugs = ['nerat']

for dr in drugs:
    
    drs_median_drug = {}
    
    for dl in range(10):
        
        drs_median_dose = {}
        
        tp_drs = [drs_all[dr]['d'+str(dl)]['r'+str(rep+1)]['tout'] for rep in range(10)]
        popdyn_reps0 = [drs_all[dr]['d'+str(dl)]['r'+str(rep+1)]['cellpop'] for rep in range(10)]
        tp_all = np.array(list(itertools.chain(*tp_drs)))
        tp_all = np.unique(tp_all)
        tp_max = min([tp_drs[x][-1] for x in range(len(tp_drs))])
        tp_max_idx = np.where(tp_all == tp_max)[0][0]
        tp_all = tp_all[:tp_max_idx+1]
        popdyn_reps = []
        for rep in range(10):
            interpolator = interp1d(tp_drs[rep],popdyn_reps0[rep])
            y_new = interpolator(tp_all)
            popdyn_reps.append(y_new)
            
        popdyn_med = np.median(popdyn_reps,axis=0)
        
        drs_median_dose['cellpop'] = popdyn_med
        drs_median_dose['tout'] = tp_all
        
        drs_median_drug['d'+str(dl)] = drs_median_dose
        
    drs_median[dr] = drs_median_drug


#%%

drug = 'nerat'

for dl in range(10):
    data = drs_median[drug]['d'+str(dl)]
    xvals = data['tout']
    yvals = data['cellpop']
    plt.plot(xvals/3600,yvals,label=str(doses_all[dl]))
    
plt.xlim(0,72)
plt.ylim(0,200)
plt.legend(bbox_to_anchor=(1.05,1.0))
plt.show()

#%% drs3 summarize
drs3_all = {}

#%%
drug = 'trame'

drs_drug = {}

for d in range(10):
    
    drs_dose = drs_summarize(drug, d,output_drs3)
    drs_drug['d'+str(d)] = drs_dose
    
drs3_all[drug] = drs_drug

#%%

pickle.dump(drs3_all, open(os.path.join(wd,'output','in_silico_drs3_summary','drs_summary.pkl'),'wb'))

#%% drs3-summary/rep-median
drs_all = drs3_all

drs_median = {}

drugs = list(drs_all.keys())

for dr in drugs:
    
    drs_median_drug = {}
    
    for dl in range(10):
        
        drs_median_dose = {}
        
        tp_drs = [drs_all[dr]['d'+str(dl)]['r'+str(rep+1)]['tout'] for rep in range(10)]
        popdyn_reps0 = [drs_all[dr]['d'+str(dl)]['r'+str(rep+1)]['cellpop'] for rep in range(10)]
        tp_all = np.array(list(itertools.chain(*tp_drs)))
        tp_all = np.unique(tp_all)
        tp_max = min([tp_drs[x][-1] for x in range(len(tp_drs))])
        tp_max_idx = np.where(tp_all == tp_max)[0][0]
        tp_all = tp_all[:tp_max_idx+1]
        popdyn_reps = []
        for rep in range(10):
            interpolator = interp1d(tp_drs[rep],popdyn_reps0[rep])
            y_new = interpolator(tp_all)
            popdyn_reps.append(y_new)
            
        popdyn_med = np.median(popdyn_reps,axis=0)
        
        drs_median_dose['cellpop'] = popdyn_med
        drs_median_dose['tout'] = tp_all
        
        drs_median_drug['d'+str(dl)] = drs_median_dose
        
    drs_median[dr] = drs_median_drug

#%%
pickle.dump(drs_median, open(os.path.join(wd,'output','in_silico_drs3_summary','drs_median.pkl'),'wb'))

#%% drs3-median plots

drug = drugs[1]

for dl in range(10):
    
    dose = doses_all[dl]
    
    x_dl = drs_median[drug]['d'+str(dl)]['tout']/3600
    y_dl = drs_median[drug]['d'+str(dl)]['cellpop']
    
    plt.plot(x_dl,y_dl,label=str(dose))
    

plt.legend(bbox_to_anchor=(1.05,1.0))
plt.ylim(0,550)
plt.xlim(0,72)
plt.xlabel('Time (h)')
plt.ylabel('# of cells')
plt.title(drug)
plt.show()


#%%

drs_summary= drs3_all

def gr_calc_row3 (drug,time_h,dl,rep,drs_summary_dict,cell_line='mcf10a_sim'):
    
    dose = doses_all[dl]
    
    dose_dict = drs_summary_dict[str(drug)]['d'+str(dl)]['r'+str(rep+1)]
    # ctrl_dict = drs_summary_dict[str(drug)]['d0']['r'+str(rep+1)]
    
    cellpop = dose_dict['cellpop']
    tout = dose_dict['tout']
    
    # cellpop_ctrl = ctrl_dict['cellpop']
    # tout_ctrl = ctrl_dict['tout']
    
    interpolator = interp1d(tout,cellpop)
    # interpolator_ctrl = interp1d(tout_ctrl,cellpop_ctrl)
    
    cell_count = interpolator(time_h*3600)
    
    # cell_count_ctrl = interpolator_ctrl(time_h*3600)
    
    cell_count_time0 = dose_dict['cellpop'][0]
    
    # drugs_all = list(drs_summary_dict.keys())
 

    # drugs_all = list(drs_summary_dict.keys())
    drugs_all = ['nerat']
    
    ctrl_pool = []
    
    for dr in drugs_all:
        for rp in range(1,len(drs_summary_dict[dr]['d0'].keys())+1):
            rp_dict = drs_summary_dict[dr]['d0']['r'+str(rp)]
            rp_cellpop = np.array(rp_dict['cellpop'])
            rp_tout = np.array(rp_dict['tout'])
            
            rp_interp = interp1d(rp_tout,rp_cellpop)
            rp_cellcount = float(rp_interp(time_h*3600))
            ctrl_pool.append(rp_cellcount)
            
    dose_pool = []
    
    for rp in range(1,len(drs_summary_dict[drug]['d'+str(dl)].keys())+1):
        rp_dose_dict = drs_summary_dict[drug]['d'+str(dl)]['r'+str(rp)]
        rp_dose_cellpop = np.array(rp_dose_dict['cellpop'])
        rp_dose_tout = np.array(rp_dose_dict['tout'])
        
        rp_dose_interp = interp1d(rp_dose_tout,rp_dose_cellpop)
        rp_dose_cellcount = float(rp_dose_interp(time_h*3600))
        dose_pool.append(rp_dose_cellcount)    
 
    cell_count_percentile = percentileofscore(dose_pool, cell_count)
    cell_count_ctrl = np.percentile(ctrl_pool,cell_count_percentile)


    new_row = {}
    new_row['cell_line'] = cell_line
    new_row['agent'] = str(drug)
    new_row['timepoint'] = time_h
    new_row['concentration'] = dose
    new_row['cell_count'] = cell_count
    new_row['cell_count__ctrl'] = cell_count_ctrl
    new_row['cell_count__time0'] = cell_count_time0 
    
    
    return new_row

#%%

time_h = [48,72]
drs_summary_full = drs_summary

drs3_grcalc3 = pd.DataFrame(data=None,columns=['cell_line','agent','timepoint','concentration','cell_count','cell_count__ctrl','cell_count__time0'])



# drug = 'nerat'
# for dl in range(1,len(drs_summary_full[str(drug)].keys())):
#     for rep in range(len(drs_summary_full[str(drug)]['d'+str(dl)].keys())):
#         for t in time_h:
#             new_row = gr_calc_row3(drug,t,dl,rep,drs_summary_full)
#             drs_grcalc3_nerat = drs_grcalc3_nerat.append(new_row, ignore_index=True)

for drug in drs_summary_full.keys():
    for dl in range(1,len(drs_summary_full[str(drug)].keys())):
        for rep in range(len(drs_summary_full[str(drug)]['d'+str(dl)].keys())):
            for t in time_h:
                new_row = gr_calc_row3(drug,t,dl,rep,drs_summary_full)
                drs3_grcalc3 = drs3_grcalc3.append(new_row, ignore_index=True)


drs3_grcalc3.to_csv(os.path.join(wd,'output','in_silico_drs3_summary','drs3_grcalc3.tsv'),sep='\t',index=False)

#%% drs3 gr plot vs exp

drugs_exp = ['Alpelisib','Neratinib','Trametinib','Palbociclib']
grv_exp = {'grv_c1a1':grv_c1a1,'grv_c1a2':grv_c1a2,'grv_c1b':grv_c1b,'grv_c1c':grv_c1c,'grv_c2':grv_c2,'grv_c3':grv_c3,'grv_c4':grv_c4,'grv_c5':grv_c5}


grv_exp_df = pd.DataFrame()

for key in grv_exp.keys():
    grv_df = grv_exp[key]
    grv_df['source'] = np.ones(np.shape(grv_df)[0])*np.nan
    grv_df['source'] = key
    grv_exp_df = grv_exp_df.append(grv_df, ignore_index=True)
    

grv_exp_compare = pd.DataFrame()

# grv_c2[math.isclose(grv_c2['concentration'],doses_all[3],abs_tol=1e-4)]

for dose_sim in doses_all[1:]:
    for row_n in range(np.shape(grv_exp_df)[0]):
        grv_exp_row = grv_exp_df.iloc[row_n,:].copy()
        if grv_exp_row['agent'] in drugs_exp:
            if math.isclose(float(grv_exp_row['concentration']),dose_sim,abs_tol=1e-4):
                grv_exp_row['dose_sim'] = dose_sim
                grv_exp_compare = grv_exp_compare.append(grv_exp_row,ignore_index=True)

#%% gr compare/drs3
drs_all_gr = pd.read_csv(os.path.join(wd,'output','in_silico_drs3_summary','drs3_grcalc3_grc.tsv'),sep='\t')

plt.rcParams["lines.linewidth"] = 2

drug_compare = drugs_exp[3]

drug_sim = drug_compare.lower()[:5]

grv_exp_drug = grv_exp_compare[grv_exp_compare['agent']==drug_compare]

x_values = doses_all[1:]

y_values_all_exp = [grv_exp_drug[grv_exp_drug['dose_sim']==xval]['GRvalue'].values for xval in x_values]

y_values_median_exp = [np.median(ys) for ys in y_values_all_exp]

if sum(np.isnan(y_values_median_exp))>0:
    
    nan_idx = np.where(np.isnan(y_values_median_exp))[0][0]
    
    x_values_actual = x_values[:nan_idx]
    x_values = x_values_actual
    
    y_values_all_exp = [grv_exp_drug[grv_exp_drug['dose_sim']==xval]['GRvalue'].values for xval in x_values]

    y_values_median_exp = [np.median(ys) for ys in y_values_all_exp]
    
    


y_min_exp = [np.min(ys) for ys in y_values_all_exp]

y_max_exp = [np.max(ys) for ys in y_values_all_exp]

y_err_min_exp = [y_values_median_exp[dl] - y_min_exp[dl] for dl in range(len(y_values_median_exp))]

y_err_max_exp = [y_max_exp[dl] - y_values_median_exp[dl] for dl in range(len(y_values_median_exp))]

yerror_exp = [y_err_min_exp,y_err_max_exp]


#%

tp = 72

drs_filtered = drs_all_gr[drs_all_gr['agent']==drug_sim]
drs_filtered = drs_filtered[drs_filtered['timepoint']==tp]

y_values_all_sim = [drs_filtered[drs_filtered['concentration']==xval]['GRvalue'].values for xval in x_values]

y_values_median_sim = [np.median(ys) for ys in y_values_all_sim]

y_min_sim = [np.min(ys) for ys in y_values_all_sim]

y_max_sim = [np.max(ys) for ys in y_values_all_sim]

y_err_min_sim = [y_values_median_sim[dl] - y_min_sim[dl] for dl in range(len(y_values_median_sim))]

y_err_max_sim = [y_max_sim[dl] - y_values_median_sim[dl] for dl in range(len(y_values_median_sim))]

yerror_sim = [y_err_min_sim,y_err_max_sim]

y_min = min(min(y_min_exp),min(y_min_sim))
y_max = max(max(y_max_exp),max(y_max_sim))

#%

plt.errorbar(x_values,y_values_median_exp,yerr=yerror_exp,fmt='o-',capsize=5,label='Experiment')
plt.errorbar(x_values,y_values_median_sim,yerr=yerror_sim,fmt='o-',capsize=5,label='Simulation')

plt.xscale('log')
plt.ylim(y_min*1.25,y_max*1.25)
# plt.ylim(-1,1.8)
plt.title(str(drug_compare),font='Arial')
plt.ylabel('GR value',font='Arial')
plt.xlabel('Concentration (uM)',font='Arial')
plt.legend(loc='lower left')
plt.show()


#%% Determine cell cycle progress

lapat1_0 = drs_dict(output_main,'lapat',1,0)
lapat1_0.timecourse_lin(11,'Md')

xout_test = lapat1_0.results['output_g1']['11']['output']['xoutS']

tout_test = lapat1_0.results['output_g1']['11']['output']['tout']

sp_md = xout_test[:,species_all.index('Md')]
sp_me = xout_test[:,species_all.index('Me')]
sp_ma = xout_test[:,species_all.index('Ma')]
sp_mb = xout_test[:,species_all.index('Mb')]

cyc_eab = np.array([sp_me/30,sp_ma/25,sp_mb/43])

cc = np.array([np.mean(cyc_eab[:,k]) for k in range(np.shape(cyc_eab)[1])])

cc_xmpl = np.loadtxt(os.path.join(wd,'input_files','cc_eab_mean.txt'),delimiter='\t')
tout_xmpl = np.loadtxt(os.path.join(wd,'input_files','cc_eab_tout.txt'),delimiter='\t')

plt.plot(tout_test/3600,cc)
plt.plot(tout_xmpl/3600,cc_xmpl)
plt.show()

#%%
diff = 84030

# diff = 299670

plt.plot(tout_xmpl/3600,cc_xmpl)
plt.plot((tout_test+diff)/3600,cc)
plt.show()



#%%
cc_interpolator = interp1d(tout_xmpl,cc_xmpl)

tout_test2 = tout_test+diff

cc_diff = cc-cc_interpolator(tout_test2)


#%%

cc_interpolator = interp1d(tout_xmpl,cc_xmpl)

# tout_test2 = tout_test+diff

# cc_diff = cc-cc_interpolator(tout_test2)

len_xmpl = int(len(tout_xmpl)*0.6)

t_diff_start = 0
t_diff_step = 30

t_diff_range = np.arange(t_diff_start, t_diff_start+t_diff_step*len_xmpl,t_diff_step)

cc_diff_sums = []

for t_diff in t_diff_range:
    tout_test2 = tout_test+t_diff
    cc_diff = cc-cc_interpolator(tout_test2)
    
    cc_diff_actual = [abs(diff) for diff in cc_diff]
    
    cc_diff_sums.append(sum(cc_diff_actual))


diff = t_diff_range[np.where(np.array(cc_diff_sums)==min(cc_diff_sums))[0]][0]



#%% find troughs

b = (np.diff(np.sign(np.diff(cc_xmpl))) > 0).nonzero()[0]

tout_xmpl_dps = tout_xmpl[b[cc_xmpl[b] < 0.2]]

tout_xmpl_dps_h = tout_xmpl_dps/3600

# tout_xmpl_dps = tout_xmpl_dps[tout_xmpl_dps>15*3600]

# tout_xmpl_dps = tout_xmpl[(np.diff(np.sign(np.diff(cc_xmpl))) > 0).nonzero()[0]]

# tout_xmpl_dps = tout_xmpl_dps[tout_xmpl_dps>15*3600]

tout_xmpl_cps = np.append(np.zeros(1),tout_xmpl_dps)

# diff = 84030

t_div_end = np.where(tout_xmpl_cps>=diff)[0][0]

t_div_start = t_div_end - 1

cc_progress = (diff-tout_xmpl_cps[t_div_start])/(tout_xmpl_cps[t_div_end]-tout_xmpl_cps[t_div_start])

#%% load dose dict

palbo1_5 = drs_dict(output_main,'palbo',1,5)
palbo1_6 = drs_dict(output_main,'palbo',1,6)

#%% cc_progress loop

# dose_dict = lapat1_0
dose_dict = palbo1_5

cc_progress_dose = []

cc_division = np.zeros(100)

t_diff_start = 0
t_diff_step = 30

len_xmpl = int(len(tout_xmpl)*0.4)

t_diff_range = np.arange(t_diff_start, t_diff_start+t_diff_step*len_xmpl,t_diff_step)

cc_interpolator = interp1d(tout_xmpl,cc_xmpl)

b = (np.diff(np.sign(np.diff(cc_xmpl))) > 0).nonzero()[0]

tout_xmpl_dps = tout_xmpl[b[cc_xmpl[b] < 0.2]]

tout_xmpl_dps_h = tout_xmpl_dps/3600

tout_xmpl_cps = np.append(np.zeros(1),tout_xmpl_dps)

for g1c in range(100):
    c_idx = g1c+1
    
    
    xout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['xoutS']
    
    tout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['tout']
    
    # sp_mb = xout_test[:,species_all.index('Mb')]
    
    sp_md = xout_test[:,species_all.index('Md')]
    sp_me = xout_test[:,species_all.index('Me')]
    sp_ma = xout_test[:,species_all.index('Ma')]
    sp_mb = xout_test[:,species_all.index('Mb')]
    
    cyc_eab = np.array([sp_me/30,sp_ma/25,sp_mb/43])

    cc = np.array([np.mean(cyc_eab[:,k]) for k in range(np.shape(cyc_eab)[1])])
 
    
    cc_diff_sums = []
    
    for t_diff in t_diff_range:
        tout_test2 = tout_test+t_diff
        # mb_diff = sp_mb-mb_interpolator(tout_test2)
        
        # mb_diff_actual = [abs(diff) for diff in mb_diff]
        
        # mb_diff_sums.append(sum(mb_diff_actual))
        
        cc_diff = cc-cc_interpolator(tout_test2)
    
        cc_diff_actual = [abs(diff) for diff in cc_diff]
        
        cc_diff_sums.append(sum(cc_diff_actual))


    diff = t_diff_range[np.where(np.array(cc_diff_sums)==min(cc_diff_sums))[0]][0]


    # diff = t_diff_range[np.where(np.array(mb_diff_sums)==min(mb_diff_sums))[0]][0]
    
    t_div_end = np.where(tout_xmpl_cps>=diff)[0][0]

    t_div_start = t_div_end - 1
    
    cell_flagA = xout_test[:,species_all.index('cPARP')][0]/xout_test[:,species_all.index('PARP')][0]
    
    if cell_flagA >= 1:
        
        cc_progress = np.nan
    
    elif diff/3600 <= 12:
        cc_progress = 0

    else:
        cc_progress = (diff-tout_xmpl_cps[t_div_start])/(tout_xmpl_cps[t_div_end]-tout_xmpl_cps[t_div_start])
    
    # plt.plot(tout_xmpl/3600,cc_xmpl)
    # plt.plot((tout_test+diff)/3600,cc)
    # plt.title(str(c_idx)+', cc_progress = '+str(cc_progress))
    # plt.show()
    cc_progress_dose.append(cc_progress)
    if len(dose_dict.get_desc(c_idx)['g2'])>0:
        cc_division[g1c] = 1.0
    
#%%

plt.scatter(cc_progress_dose,cc_division)

#%% all doses

for dl in range(len(doses_all)):
    

    
    dose_dict = drs_dict(output_main,'palbo',1,dl)
    
    cc_progress_dose = []
    
    cc_division = np.zeros(100)
    
    t_diff_start = 0
    t_diff_step = 30
    
    len_xmpl = int(len(tout_xmpl)*0.4)
    
    t_diff_range = np.arange(t_diff_start, t_diff_start+t_diff_step*len_xmpl,t_diff_step)
    
    cc_interpolator = interp1d(tout_xmpl,cc_xmpl)
    
    b = (np.diff(np.sign(np.diff(cc_xmpl))) > 0).nonzero()[0]
    
    tout_xmpl_dps = tout_xmpl[b[cc_xmpl[b] < 0.2]]
    
    tout_xmpl_dps_h = tout_xmpl_dps/3600
    
    tout_xmpl_cps = np.append(np.zeros(1),tout_xmpl_dps)
    
    for g1c in range(100):
        c_idx = g1c+1
        
        
        xout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['xoutS']
        
        tout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['tout']
        
        # sp_mb = xout_test[:,species_all.index('Mb')]
        
        sp_md = xout_test[:,species_all.index('Md')]
        sp_me = xout_test[:,species_all.index('Me')]
        sp_ma = xout_test[:,species_all.index('Ma')]
        sp_mb = xout_test[:,species_all.index('Mb')]
        
        cyc_eab = np.array([sp_me/30,sp_ma/25,sp_mb/43])
    
        cc = np.array([np.mean(cyc_eab[:,k]) for k in range(np.shape(cyc_eab)[1])])
     
        
        cc_diff_sums = []
        
        for t_diff in t_diff_range:
            tout_test2 = tout_test+t_diff
            # mb_diff = sp_mb-mb_interpolator(tout_test2)
            
            # mb_diff_actual = [abs(diff) for diff in mb_diff]
            
            # mb_diff_sums.append(sum(mb_diff_actual))
            
            cc_diff = cc-cc_interpolator(tout_test2)
        
            cc_diff_actual = [abs(diff) for diff in cc_diff]
            
            cc_diff_sums.append(sum(cc_diff_actual))
    
    
        diff = t_diff_range[np.where(np.array(cc_diff_sums)==min(cc_diff_sums))[0]][0]
    
    
        # diff = t_diff_range[np.where(np.array(mb_diff_sums)==min(mb_diff_sums))[0]][0]
        
        t_div_end = np.where(tout_xmpl_cps>=diff)[0][0]
    
        t_div_start = t_div_end - 1
        
        cell_flagA = xout_test[:,species_all.index('cPARP')][0]/xout_test[:,species_all.index('PARP')][0]
        
        if cell_flagA >= 1:
            
            cc_progress = np.nan
        
        elif diff/3600 <= 12:
            cc_progress = 0
    
        else:
            cc_progress = (diff-tout_xmpl_cps[t_div_start])/(tout_xmpl_cps[t_div_end]-tout_xmpl_cps[t_div_start])
        
        # plt.plot(tout_xmpl/3600,cc_xmpl)
        # plt.plot((tout_test+diff)/3600,cc)
        # plt.title(str(c_idx)+', cc_progress = '+str(cc_progress))
        # plt.show()
        cc_progress_dose.append(cc_progress)
        if len(dose_dict.get_desc(c_idx)['g2'])>0:
            cc_division[g1c] = 1.0
        elif cc_progress>0 and len(dose_dict.get_desc(c_idx)['g2'])==0:
            print('dose_lvl: '+str(dl)+', cell: '+str(c_idx))
    
    plt.scatter(cc_progress_dose,cc_division)
    plt.show()
    
    
#%% all reps/dose specific

dl = 0
# div_gen = 2

cc_progress_dl = []
g1_division_dl = []
g2_division_dl = []
    
for rep in range(1,11):
    
    dose_dict = drs_dict(output_main,'palbo',rep,dl)
    
    cc_progress_rep = []
    
    g1_division = np.zeros(100)
    g2_division = np.zeros(100)
    
    t_diff_start = 0
    t_diff_step = 30
    
    len_xmpl = int(len(tout_xmpl)*0.4)
    
    t_diff_range = np.arange(t_diff_start, t_diff_start+t_diff_step*len_xmpl,t_diff_step)
    
    cc_interpolator = interp1d(tout_xmpl,cc_xmpl)
    
    b = (np.diff(np.sign(np.diff(cc_xmpl))) > 0).nonzero()[0]
    
    tout_xmpl_dps = tout_xmpl[b[cc_xmpl[b] < 0.2]]
    
    tout_xmpl_dps_h = tout_xmpl_dps/3600
    
    tout_xmpl_cps = np.append(np.zeros(1),tout_xmpl_dps)
    
    for g1c in range(100):
        c_idx = g1c+1
        
        
        xout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['xoutS']
        
        tout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['tout']
        
        # sp_mb = xout_test[:,species_all.index('Mb')]
        
        sp_md = xout_test[:,species_all.index('Md')]
        sp_me = xout_test[:,species_all.index('Me')]
        sp_ma = xout_test[:,species_all.index('Ma')]
        sp_mb = xout_test[:,species_all.index('Mb')]
        
        cyc_eab = np.array([sp_me/30,sp_ma/25,sp_mb/43])
    
        cc = np.array([np.mean(cyc_eab[:,k]) for k in range(np.shape(cyc_eab)[1])])
     
        
        cc_diff_sums = []
        
        for t_diff in t_diff_range:
            tout_test2 = tout_test+t_diff
            # mb_diff = sp_mb-mb_interpolator(tout_test2)
            
            # mb_diff_actual = [abs(diff) for diff in mb_diff]
            
            # mb_diff_sums.append(sum(mb_diff_actual))
            
            cc_diff = cc-cc_interpolator(tout_test2)
        
            cc_diff_actual = [abs(diff) for diff in cc_diff]
            
            cc_diff_sums.append(sum(cc_diff_actual))
    
    
        diff = t_diff_range[np.where(np.array(cc_diff_sums)==min(cc_diff_sums))[0]][0]
    
    
        # diff = t_diff_range[np.where(np.array(mb_diff_sums)==min(mb_diff_sums))[0]][0]
        
        t_div_end = np.where(tout_xmpl_cps>=diff)[0][0]
    
        t_div_start = t_div_end - 1
        
        cell_flagA = xout_test[:,species_all.index('cPARP')][0]/xout_test[:,species_all.index('PARP')][0]
        
        if cell_flagA >= 1:
            
            cc_progress = np.nan
        
        elif diff/3600 <= 12:
            cc_progress = 0
    
        else:
            cc_progress = (diff-tout_xmpl_cps[t_div_start])/(tout_xmpl_cps[t_div_end]-tout_xmpl_cps[t_div_start])
        
        # plt.plot(tout_xmpl/3600,cc_xmpl)
        # plt.plot((tout_test+diff)/3600,cc)
        # plt.title(str(c_idx)+', cc_progress = '+str(cc_progress))
        # plt.show()
        cc_progress_rep.append(cc_progress)
        if 'g2' in dose_dict.get_desc(c_idx).keys():
            if len(dose_dict.get_desc(c_idx)['g2'])>0:
                g1_division[g1c] = 1.0
        if 'g3' in dose_dict.get_desc(c_idx).keys():
            if len(dose_dict.get_desc(c_idx)['g3'])>0:
                g2_division[g1c] = 1.0
        # elif cc_progress>0 and len(dose_dict.get_desc(c_idx)['g2'])==0:
        #     print('dose_lvl: '+str(dl)+', cell: '+str(c_idx))
    
    cc_progress_dl.append(cc_progress_rep)
    g1_division_dl.append(g1_division)
    g2_division_dl.append(g2_division)
    
    
    
    # plt.scatter(cc_progress_rep,cc_division)
    # plt.show()


cpd = [item for sublist in cc_progress_dl for item in sublist]
cdd_g1 = [item for sublist in g1_division_dl for item in sublist]
cdd_g2 = [item for sublist in g2_division_dl for item in sublist]

palbo_cd = pd.DataFrame()

palbo_cd['cc_progress'] = cpd
palbo_cd['g1_div'] = cdd_g1
palbo_cd['g2_div'] = cdd_g2

palbo_cd.to_csv(os.path.join(wd,'output','in_silico_drs_summary','palbo_cc','palbo_cd_dl_'+str(dl)+'.txt'),sep='\t')

#%%

palbo10_3 = drs_dict(output_main,'palbo',10,3)



#%%

cpd2 = [item for sublist in cc_progress_dl for item in sublist]
cdd2 = [item for sublist in cc_division_dl for item in sublist]

cdp2_div = np.array(cpd2)[np.array(cdd2)==1]
cdp2_nondiv = np.array(cpd2)[np.array(cdd2)==0]

plt.hist(cdp2_div,density=True)
plt.hist(cdp2_nondiv,density=True)

pd.Series(cdp2_div).plot(kind='density')

pd.Series(cdp2_div).plot.kde(bw_method=1)

#%%
dl = 5
palbo_cd_dl = pd.read_csv(os.path.join(wd,'output','in_silico_drs_summary','palbo_cc','palbo_cd_dl_'+str(dl)+'.txt'),sep='\t',index_col=0,header=0)

cdp_g1div = palbo_cd_dl['cc_progress'][palbo_cd_dl['g1_div']==1.0]
cdp_g2div = palbo_cd_dl['cc_progress'][palbo_cd_dl['g2_div']==1.0]

cdp_g1nondiv = palbo_cd_dl['cc_progress'][palbo_cd_dl['g1_div']==0.0]
cdp_g1nondiv = cdp_g1nondiv[~np.isnan(cdp_g1nondiv)]
sns.kdeplot(cdp_g1div,label='dividing')
# sns.kdeplot(cdp_g2div,label='g2')
sns.kdeplot(cdp_g1nondiv,label='non-dividing')
plt.legend()
plt.show()

#%%
sns.kdeplot(cdp2_div)
# sns.kdeplot(cdp2_nondiv)
plt.xlabel('Cell cycle progress')
# plt.ylabel('Probability density of cell division (gen '+str(div_gen)+')')
plt.title('Palbociclib: '+str(doses_all[dl])+' uM')
plt.show()

#%%

sns.kdeplot(cdp2_nondiv)

ax = sns.kdeplot(cdp_g1div)

kde_lines = ax.get_lines()[-1]
kde_x, kde_y = kde_lines.get_data()

kde_interpolator = interp1d(kde_x,kde_y)

cc_progress_range = np.linspace(0,1,101)

cc_progress_intervals = np.linspace(0,1,11)

cc_progress_density = kde_interpolator(cc_progress_range)

cc_progress_area = np.trapz(cc_progress_density,dx=0.01)

cc_progress_range[np.logical_and(cc_progress_range>=0.1, cc_progress_range<0.2)]

for k in range(len(cc_progress_intervals[:-1])):
    lower = float(cc_progress_intervals[k])
    upper = float(cc_progress_intervals[k+1])
    interval = cc_progress_range[np.logical_and(cc_progress_range>=lower, cc_progress_range<upper)]
    density_interval = kde_interpolator(interval)
    cc_progress_prob_interval = np.trapz(density_interval,dx=0.01)/cc_progress_area
    print(cc_progress_prob_interval)




#%% div probability estimate

def kde_function(cdp_array):
    ax = sns.kdeplot(cdp_array)
    kde_lines = ax.get_lines()[-1]
    kde_x,kde_y = kde_lines.get_data()
    kde_max = max(kde_x)
    kde_interpolator = interp1d(kde_x,kde_y)
    cc_progress_range = np.linspace(0,1,101)
    cc_progress_range = cc_progress_range[cc_progress_range<kde_max]
    cc_progress_density = kde_interpolator(cc_progress_range)
    cc_progress_area = np.trapz(cc_progress_density,dx=0.01)
    
    return kde_interpolator,cc_progress_area,kde_max

#%%


dl = 1
palbo_cd_dl = pd.read_csv(os.path.join(wd,'output','in_silico_drs_summary','palbo_cc','palbo_cd_dl_'+str(dl)+'.txt'),sep='\t',index_col=0,header=0)

cdp_g1div = palbo_cd_dl['cc_progress'][palbo_cd_dl['g1_div']==1.0]
cdp_g2div = palbo_cd_dl['cc_progress'][palbo_cd_dl['g2_div']==1.0]

cdp_g1nondiv = palbo_cd_dl['cc_progress'][palbo_cd_dl['g1_div']==0.0]
cdp_g1nondiv = cdp_g1nondiv[~np.isnan(cdp_g1nondiv)]

cdp_g2nondiv = palbo_cd_dl['cc_progress'][palbo_cd_dl['g2_div']==0.0]
cdp_g2nondiv = cdp_g2nondiv[~np.isnan(cdp_g2nondiv)]

n_g1div = len(cdp_g1div)
n_g1nondiv = len(cdp_g1nondiv)
n_g2div = len(cdp_g2div)
n_g2nondiv = len(cdp_g2nondiv)

g1div_intp,g1div_area,g1div_kdemax = kde_function(cdp_g1div)
g1nondiv_intp,g1nondiv_area,g1nondiv_kdemax = kde_function(cdp_g1nondiv)

g2div_intp,g2div_area,g2div_kdemax = kde_function(cdp_g2div)
g2nondiv_intp,g2nondiv_area,g2nondiv_kdemax = kde_function(cdp_g2nondiv)



# cc_progress_range = np.linspace(0,1,101)
cc_progress_range = np.linspace(0,1,1001)

# cc_progress_intervals = np.linspace(0,1,11)
cc_progress_intervals = np.linspace(0,1,101)

cc_progress_density = kde_interpolator(cc_progress_range)

cc_progress_area = np.trapz(cc_progress_density,dx=0.01)

cc_progress_g1df = pd.DataFrame(data=np.ones(len(cc_progress_range))*np.nan,index=cc_progress_range)
cc_progress_g2df = pd.DataFrame(data=np.ones(len(cc_progress_range))*np.nan,index=cc_progress_range)

cc_progress_g1df = cc_progress_g1df[cc_progress_g1df.index<min(g1div_kdemax,g1nondiv_kdemax)]
cc_progress_g2df = cc_progress_g2df[cc_progress_g2df.index<min(g2div_kdemax,g2nondiv_kdemax)]

for k in range(len(cc_progress_intervals[:-1])):
    lower = float(cc_progress_intervals[k])
    upper = float(cc_progress_intervals[k+1])
    interval = cc_progress_range[np.logical_and(cc_progress_range>=lower, cc_progress_range<upper)]
    interval = interval[interval<min(g1div_kdemax,g1nondiv_kdemax)]
    
    g2_interval = cc_progress_range[np.logical_and(cc_progress_range>=lower, cc_progress_range<upper)]
    g2_interval = g2_interval[g2_interval<min(g2div_kdemax,g2nondiv_kdemax)]

    
    g1div_density_interval = g1div_intp(interval)
    g1nondiv_density_interval = g1nondiv_intp(interval)
    g1div_frac_interval = np.trapz(g1div_density_interval,dx=0.01)/g1div_area
    g1nondiv_frac_interval = np.trapz(g1nondiv_density_interval,dx=0.01)/g1nondiv_area
    
    g2div_density_interval = g2div_intp(g2_interval)
    g2nondiv_density_interval = g2nondiv_intp(g2_interval)
    g2div_frac_interval = np.trapz(g2div_density_interval,dx=0.01)/g2div_area
    g2nondiv_frac_interval = np.trapz(g2nondiv_density_interval,dx=0.01)/g2nondiv_area    
    
    
    
    g1div_prob_interval = n_g1div*g1div_frac_interval/(n_g1div*g1div_frac_interval+n_g1nondiv*g1nondiv_frac_interval)
    
    g2div_prob_interval = n_g2div*g2div_frac_interval/(n_g2div*g2div_frac_interval+n_g2nondiv*g2nondiv_frac_interval)
    
    for i in interval:
        cc_progress_g1df.loc[i,0] = g1div_prob_interval
        
    for i in g2_interval:
        cc_progress_g2df.loc[i,0] = g2div_prob_interval

    
plt.figure()
plt.plot(cc_progress_g1df.index,cc_progress_g1df[0],label='Generation 1')
plt.plot(cc_progress_g2df.index,cc_progress_g2df[0],label='Generation 2')
plt.title('Palbociclib: '+str(doses_all[dl])+' uM',fontsize=15)
plt.xlabel('Fractional cell cycle progress at the time of drug dose',fontsize=15)
plt.ylabel('Cell division probability',fontsize=15)
plt.ylim(0,1.2)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend()
plt.show()

# plt.figure()
# plt.plot(cc_progress_g2df.index,cc_progress_g2df[0],label='gen 2')
# plt.title('Palbociclib: '+str(doses_all[dl])+' uM')
# plt.xlabel('Cell cycle progress at t=0')
# plt.ylabel('Gen2 division probability')
# plt.ylim(0,1.1)
# plt.show()

#%%




#%% cc progress single cell

c_idx = 59
dose_dict = palbo1_5

   
xout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['xoutS']

tout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['tout']

# sp_mb = xout_test[:,species_all.index('Mb')]

sp_md = xout_test[:,species_all.index('Md')]
sp_me = xout_test[:,species_all.index('Me')]
sp_ma = xout_test[:,species_all.index('Ma')]
sp_mb = xout_test[:,species_all.index('Mb')]

cyc_eab = np.array([sp_me/30,sp_ma/25,sp_mb/43])

cc = np.array([np.mean(cyc_eab[:,k]) for k in range(np.shape(cyc_eab)[1])])
 

cc_diff_sums = []

len_xmpl = int(len(tout_xmpl)*0.8)

t_diff_range = np.arange(t_diff_start, t_diff_start+t_diff_step*len_xmpl,t_diff_step)

for t_diff in t_diff_range:
    tout_test2 = tout_test+t_diff
    # mb_diff = sp_mb-mb_interpolator(tout_test2)
    
    # mb_diff_actual = [abs(diff) for diff in mb_diff]
    
    # mb_diff_sums.append(sum(mb_diff_actual))
    
    cc_diff = cc-cc_interpolator(tout_test2)

    cc_diff_actual = [abs(diff) for diff in cc_diff]
    
    cc_diff_sums.append(sum(cc_diff_actual))


diff = t_diff_range[np.where(np.array(cc_diff_sums)==min(cc_diff_sums))[0]][0]


# diff = t_diff_range[np.where(np.array(mb_diff_sums)==min(mb_diff_sums))[0]][0]

t_div_end = np.where(tout_xmpl_cps>=diff)[0][0]

t_div_start = t_div_end - 1


if diff == 0:
    cc_progress = 0
else:
    cc_progress = (diff-tout_xmpl_cps[t_div_start])/(tout_xmpl_cps[t_div_end]-tout_xmpl_cps[t_div_start])

print(cc_progress)

plt.plot(tout_xmpl/3600,cc_xmpl)
plt.plot((tout_test+diff)/3600,cc)
plt.title(str(c_idx)+', cc_progress = '+str(cc_progress))
plt.show()
#%% cc_progress based on Mb

lapat1_0 = drs_dict(output_main,'lapat',1,0)

g1c = 11

xout_test = lapat1_0.results['output_g1'][str(g1c)]['output']['xoutS']

tout_test = lapat1_0.results['output_g1'][str(g1c)]['output']['tout']

sp_mb = xout_test[:,species_all.index('Mb')]

mb_xmpl = np.loadtxt(os.path.join(wd,'input_files','cc_mb.txt'),delimiter='\t')
tout_xmpl = np.loadtxt(os.path.join(wd,'input_files','cc_mb_tout.txt'),delimiter='\t')


# plt.plot(tout_test/3600,sp_mb)
# plt.plot(tout_xmpl/3600,mb_xmpl)
# plt.show()

# diff = 84030

# diff = 299670

# plt.plot(tout_xmpl/3600,mb_xmpl)
# plt.plot((tout_test+diff)/3600,sp_mb)
# plt.show()

mb_interpolator = interp1d(tout_xmpl,mb_xmpl)


len_xmpl = int(len(tout_xmpl)*0.6)

t_diff_start = 0
t_diff_step = 30

t_diff_range = np.arange(t_diff_start, t_diff_start+t_diff_step*len_xmpl,t_diff_step)

mb_diff_sums = []

for t_diff in t_diff_range:
    tout_test2 = tout_test+t_diff
    mb_diff = sp_mb-mb_interpolator(tout_test2)
    
    mb_diff_actual = [abs(diff) for diff in mb_diff]
    
    mb_diff_sums.append(sum(mb_diff_actual))


diff = t_diff_range[np.where(np.array(mb_diff_sums)==min(mb_diff_sums))[0]][0]

#%%
b = (np.diff(np.sign(np.diff(mb_xmpl))) > 0).nonzero()[0]

tout_xmpl_dps = tout_xmpl[b[mb_xmpl[b] < 8]]

tout_xmpl_dps = np.delete(tout_xmpl_dps,0)

tout_xmpl_dps = np.delete(tout_xmpl_dps,0)

tout_xmpl_cps = np.append(np.zeros(1),tout_xmpl_dps)

tout_xmpl_cps_h = tout_xmpl_cps/3600

t_div_end = np.where(tout_xmpl_cps>=diff)[0][0]

t_div_start = t_div_end - 1

cc_progress = (diff-tout_xmpl_cps[t_div_start])/(tout_xmpl_cps[t_div_end]-tout_xmpl_cps[t_div_start])

#%%

# palbo1_7 = drs_dict(output_main,'palbo',1,7)
# palbo1_3 = drs_dict(output_main,'palbo',1,3)
palbo1_4 = drs_dict(output_main,'palbo',1,4)
# palbo1_5 = drs_dict(output_main,'palbo',1,5)
# dose_dict = lapat1_0
# dose_dict = palbo1_7
dose_dict = palbo1_4

mb_xmpl = np.loadtxt(os.path.join(wd,'input_files','cc_mb.txt'),delimiter='\t')
tout_xmpl = np.loadtxt(os.path.join(wd,'input_files','cc_mb_tout.txt'),delimiter='\t')
mb_interpolator = interp1d(tout_xmpl,mb_xmpl)

b = (np.diff(np.sign(np.diff(mb_xmpl))) > 0).nonzero()[0]

tout_xmpl_dps = tout_xmpl[b[mb_xmpl[b] < 8]]

tout_xmpl_dps = np.delete(tout_xmpl_dps,0)

tout_xmpl_dps = np.delete(tout_xmpl_dps,0)

tout_xmpl_cps = np.append(np.zeros(1),tout_xmpl_dps)

tout_xmpl_cps_h = tout_xmpl_cps/3600


len_xmpl = int(len(tout_xmpl)*0.4)

t_diff_start = 0
t_diff_step = 30

t_diff_range = np.arange(t_diff_start, t_diff_start+t_diff_step*len_xmpl,t_diff_step)

for g1c in range(100):
    c_idx = g1c+1
    
    
    xout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['xoutS']
    
    tout_test = dose_dict.results['output_g1'][str(c_idx)]['output']['tout']
    
    sp_mb = xout_test[:,species_all.index('Mb')]
 
    
    mb_diff_sums = []
    
    for t_diff in t_diff_range:
        tout_test2 = tout_test+t_diff
        mb_diff = sp_mb-mb_interpolator(tout_test2)
        
        mb_diff_actual = [abs(diff) for diff in mb_diff]
        
        mb_diff_sums.append(sum(mb_diff_actual))

    

    diff = t_diff_range[np.where(np.array(mb_diff_sums)==min(mb_diff_sums))[0]][0]
    
    t_div_end = np.where(tout_xmpl_cps>=diff)[0][0]

    t_div_start = t_div_end - 1

    cc_progress = (diff-tout_xmpl_cps[t_div_start])/(tout_xmpl_cps[t_div_end]-tout_xmpl_cps[t_div_start])
    
    plt.plot(tout_xmpl/3600,mb_xmpl)
    plt.plot((tout_test+diff)/3600,sp_mb)
    plt.title(str(c_idx)+', cc_progress = '+str(cc_progress))
    plt.show()







#%%




#%%
lapat1_0 = drs_dict(output_main,'lapat',1,0)
lapat1_1 = drs_dict(output_main,'lapat',1,1)
lapat2_3 = drs_dict(output_main,'lapat',2,3)
lapat3_5 = drs_dict(output_main,'lapat',3,5)

palbo1_5 = drs_dict(output_main,'palbo',1,5)
#%% rep average population dynamics

lapat1_4 = drs_dict(output_main,'lapat',1,4)
lapat2_4 = drs_dict(output_main,'lapat',2,4)
lapat3_4 = drs_dict(output_main,'lapat',3,4)

trame1_0 = drs_dict(output_main,'trame',1,0)
trame1_1 = drs_dict(output_main,'trame',1,1)
trame1_5 = drs_dict(output_main,'trame',1,5)
trame1_9 = drs_dict(output_main,'trame',1,9)

#%%

a,b,c = lapat1_4.pop_dyn()
d,e,f = lapat2_4.pop_dyn()
g,h,i = lapat3_4.pop_dyn()

#%%
plt.plot(b/3600,a)
plt.plot(e/3600,d)
plt.plot(h/3600,g)


#%%
a,b,c = lapat3_5.pop_dyn()

plt.plot(b/3600,a)

#%%

for d in range(10):
    dict0 = drs_dict(output_main,'lapat',3,d)
    a,b,c = dict0.pop_dyn()
    plt.plot(b/3600,a)
    
plt.ylim(0,450)
plt.xlim(0,72)

plt.show()

#%% multi-rep/average

lapat1_4 = drs_dict(output_main,'lapat',1,4)
lapat2_4 = drs_dict(output_main,'lapat',2,4)
lapat3_4 = drs_dict(output_main,'lapat',3,4)
lapat4_4 = drs_dict(output_main,'lapat',4,4)
lapat5_4 = drs_dict(output_main,'lapat',5,4)

#%% investigate nerat

nerat1_0 = drs_dict(output_main,'nerat',1,0)
nerat1_5 = drs_dict(output_main,'nerat',1,5)
nerat1_9 = drs_dict(output_main,'nerat',1,9)

#%% drs2

drs2_nerat1_0 = drs_dict(output_drs2,'nerat',1,0)
drs2_nerat1_5 = drs_dict(output_drs2,'nerat',1,5)
drs2_nerat1_8 = drs_dict(output_drs2,'nerat',1,8)

#%%

for dl in range(len(doses_all)):
    dose_dict = drs_dict(output_drs2,'nerat',1,dl)
    cells,tps,tout_deaths,obs_array = dose_dict.pop_dyn_obs(formula_ppERK)
    pop_obs_plot(obs_array,tps)

#%% drs3
for dl in range(len(doses_all)):
    dose_dict = drs_dict(output_drs3,'nerat',1,dl)
    cells,tps,tout_deaths,obs_array = dose_dict.pop_dyn_obs(formula_ppERK)
    pop_obs_plot(obs_array,tps)


#%%

popdyn_drs = []
tp_drs = []

for rp in range(5):
    dict0 = drs_dict(output_main,'lapat',rp+1,4)
    
    a,b,c = dict0.pop_dyn()
    
    popdyn_drs.append(a)
    tp_drs.append(b)
    
#%%

for d in range(5):
    plt.plot(tp_drs[d]/3600,popdyn_drs[d])
    
plt.xlabel('Time(h)')
plt.ylabel('Cell #')
plt.ylim(0,550)
plt.xlim(0,72)
plt.show()




#%% common timepoints

tp_all = np.array(list(itertools.chain(*tp_drs)))

tp_all = np.unique(tp_all)

popdyn_drs_new = []

for rp in range(len(popdyn_drs)):
    interpolator = interp1d(tp_drs[rp],popdyn_drs[rp])
    y_new = interpolator(tp_all)
    popdyn_drs_new.append(y_new)

#%%

for d in range(5):
    plt.plot(tp_all/3600,popdyn_drs_new[d])
plt.xlim(0,72)
plt.ylim(0,550)
plt.show()

#%% drs - popdyn/rep average
import itertools
from scipy.interpolate import interp1d

drug = 'alpel'

path_reps = os.listdir(os.path.join(output_main,'drs_'+str(drug)))

rep_count = len(path_reps)

path_rep0 = os.path.join(output_main,'drs_'+str(drug),path_reps[0])

path_doses = os.listdir(path_rep0)

doses_all = [float(x.split('_')[-1]) for x in path_doses]

doses_all.sort()

dose_levels = len(doses_all)

popdyn_drs = []

# tp_drs = []

for dl in range(dose_levels):
    popdyn_reps0 = []
    tp_reps0 = []
    for rep in range(rep_count):
        dict0 = drs_dict(output_main,drug,rep+1,dl)
        a,b,c = dict0.pop_dyn()
        popdyn_reps0.append(a)
        tp_reps0.append(b)
    tp_all = np.array(list(itertools.chain(*tp_reps0)))
    tp_all = np.unique(tp_all)
    tp_max = min([tp_reps0[x][-1] for x in range(len(tp_reps0))])
    tp_max_idx = np.where(tp_all == tp_max)[0][0]
    tp_all = tp_all[:tp_max_idx+1]
    popdyn_reps = []
    for rep in range(rep_count):
        interpolator = interp1d(tp_reps0[rep],popdyn_reps0[rep])
        y_new = interpolator(tp_all)
        popdyn_reps.append(y_new)
        
    popdyn_avg = np.average(popdyn_reps,axis=0)
    popdyn_drs.append(popdyn_avg)

#%%



#%% plot obs

ObsMat = pd.read_csv(os.path.join(wd,'input_files','Observables.txt'),sep='\t',header=0,index_col=0)
Species_doc = pd.read_csv(os.path.join(wd,'input_files','Species.txt'),sep='\t',header=0,index_col=0)
Compartment_vol = pd.read_csv(os.path.join(wd,'input_files','Compartments.txt'),sep='\t',header=0,index_col=0)

sp_erk = list(np.array(species_all)[np.where(ObsMat.loc[:,'ERK'].values)[0]])

sp_pperk = list(np.array(sp_erk)[np.where(['ppERK' in sp_erk[k] for k in range(len(sp_erk))])[0]])

formula_num = ''

for sp in sp_pperk:
    sp_comp = Species_doc.loc[sp,'compartment']
    multiplier = float(ObsMat.loc[sp,'ERK'])*float(Compartment_vol.loc[sp_comp,'volume'])
    num_item = str(sp)+'*'+str(multiplier)
    formula_num = formula_num+'+'+num_item

formula_num = formula_num[1:]

formula_den = ''

for sp in sp_erk:
    sp_comp = Species_doc.loc[sp,'compartment']
    multiplier = float(ObsMat.loc[sp,'ERK'])*float(Compartment_vol.loc[sp_comp,'volume'])
    den_item = str(sp)+'*'+str(multiplier)
    formula_den = formula_den+'+'+den_item

formula_den = formula_den[1:]

formula_ppERK = '(' + formula_num + ')/(' + formula_den + ')'


##

sp_akt = list(np.array(species_all)[np.where(ObsMat.loc[:,'AKT'].values)[0]])

sp_ppakt = list(np.array(sp_akt)[np.where(['ppAKT' in sp_akt[k] for k in range(len(sp_akt))])[0]])

formula_num = ''

for sp in sp_ppakt:
    sp_comp = Species_doc.loc[sp,'compartment']
    multiplier = float(ObsMat.loc[sp,'AKT'])*float(Compartment_vol.loc[sp_comp,'volume'])
    num_item = str(sp)+'*'+str(multiplier)
    formula_num = formula_num+'+'+num_item

formula_num = formula_num[1:]

formula_den = ''

for sp in sp_akt:
    sp_comp = Species_doc.loc[sp,'compartment']
    multiplier = float(ObsMat.loc[sp,'AKT'])*float(Compartment_vol.loc[sp_comp,'volume'])
    den_item = str(sp)+'*'+str(multiplier)
    formula_den = formula_den+'+'+den_item

formula_den = formula_den[1:]

formula_ppAKT = '(' + formula_num + ')/(' + formula_den + ')'

#%% fractional EGFR inhibition

sp_EGFR = list(np.array(species_all)[np.where(ObsMat.loc[:,'E1'].values)[0]])

sp_EGFR_nerat = list(np.array(sp_EGFR)[np.where(['nerat' in sp_EGFR[k] for k in range(len(sp_EGFR))])[0]])


formula_num_egfr = ''

for sp in sp_EGFR_nerat:
    sp_comp = Species_doc.loc[sp,'compartment']
    multiplier = float(ObsMat.loc[sp,'E1'])*float(Compartment_vol.loc[sp_comp,'volume'])
    num_item = str(sp)+'*'+str(multiplier)
    formula_num_egfr = formula_num_egfr +'+'+num_item

formula_num_egfr = formula_num_egfr[1:]

formula_den_egfr = ''

for sp in sp_EGFR:
    sp_comp = Species_doc.loc[sp,'compartment']
    multiplier = float(ObsMat.loc[sp,'E1'])*float(Compartment_vol.loc[sp_comp,'volume'])
    den_item = str(sp)+'*'+str(multiplier)
    formula_den_egfr = formula_den_egfr+'+'+den_item

formula_den_egfr = formula_den_egfr[1:]

formula_egfr = '(' + formula_num_egfr + ')/(' + formula_den_egfr + ')'

#%% palbociclib target engagement

sp_Md = list(np.array(species_all)[np.where(ObsMat.loc[:,'Cd'].values)[0]])

sp_Md.remove('Cd')

sp_Md_drug = list(np.array(sp_Md)[np.where(['palbo' in sp_Md[k] for k in range(len(sp_Md))])])

formula_num_palbo_target = ''

for sp in sp_Md_drug:
    sp_comp = Species_doc.loc[sp,'compartment']
    multiplier = float(ObsMat.loc[sp,'Cd'])*float(Compartment_vol.loc[sp_comp,'volume'])
    num_item = str(sp)+'*'+str(multiplier)
    formula_num_palbo_target = formula_num_palbo_target +'+'+num_item
    
formula_num_palbo_target = formula_num_palbo_target[1:]
    
formula_den_palbo_target = ''

for sp in sp_Md:
    sp_comp = Species_doc.loc[sp,'compartment']
    multiplier = float(ObsMat.loc[sp,'Cd'])*float(Compartment_vol.loc[sp_comp,'volume'])
    den_item = str(sp)+'*'+str(multiplier)
    formula_den_palbo_target = formula_den_palbo_target +'+'+den_item

formula_den_palbo_target = formula_den_palbo_target[1:]

formula_palbo_target = '(' + formula_num_palbo_target + ')/(' + formula_den_palbo_target + ')'



    
#%% test
lapat1_0 = drs_dict(output_main,'lapat',1,0)

lapat1_1.timecourse_lin_obs(11,formula_ppERK)
lapat1_0.timecourse_lin_obs(11,formula_ppERK)
lapat1_0.timecourse_lin_obs(11,formula_ppAKT)

nerat1_5 = drs_dict(output_main,'nerat',1,5)
nerat1_5.timecourse_lin_obs(11,formula_egfr)

#%% palbociclib story

palbo1_0 = drs_dict(output_main,'palbo',1,0)
palbo1_5 = drs_dict(output_main,'palbo',1,5)
palbo1_9 = drs_dict(output_main,'palbo',1,9)
palbo2_9 = drs_dict(output_main,'palbo',2,9)
palbo1_1 = drs_dict(output_main,'palbo',1,1)
#%% control spread

popdyn_drs = []
tp_drs = []

drugs_all = ['alpel','lapat','nerat','palbo','trame']

for drug in drugs_all:
    for rep in range(5):
        dict0 = drs_dict(output_main,drug,rep+1,0)
        a,b,c = dict0.pop_dyn()
        popdyn_drs.append(a)
        tp_drs.append(b)

# for rep in range(2):
#     dict0 = drs_dict(output_main,'trame',rep+1,0)
#     a,b,c = dict0.pop_dyn()
#     popdyn_drs.append(a)
#     tp_drs.append(b)

#%%

tp_drs_all = [item for sublist in tp_drs for item in sublist]

tp_drs_all = np.unique(tp_drs_all)

tp_max = min([max(tp) for tp in tp_drs])

tp_drs_all = tp_drs_all[tp_drs_all<=tp_max]

popdyn_drs_uniform = []

for d in range(len(tp_drs)):
    
    popdyn_intp = interp1d(tp_drs[d],popdyn_drs[d])
    
    popdyn_new = popdyn_intp(tp_drs_all)
    
    popdyn_drs_uniform.append(popdyn_new)
    
    # plt.plot(tp_drs[d]/3600,popdyn_drs[d])
# plt.xlabel('Time(h)')
# plt.ylabel('Cell #')
# plt.ylim(0,600)
# plt.xlim(0,72)
# plt.show()

popdyn_drs_uniform = np.array(popdyn_drs_uniform)

popdyn_median = [np.median(popdyn_drs_uniform[:,tp]) for tp in range(np.shape(popdyn_drs_uniform)[1])]
popdyn_std = [np.std(popdyn_drs_uniform[:,tp]) for tp in range(np.shape(popdyn_drs_uniform)[1])]

popdyn_median_intp = interp1d(popdyn_median,tp_drs_all)


for p in range(len(popdyn_drs_uniform)):
    
    plt.plot(tp_drs_all/3600,popdyn_drs_uniform[p,:],linewidth=0.5,color='grey')

plt.plot(tp_drs_all/3600,popdyn_median,linewidth=2.0,color='red',label='Median')
plt.axhline(y=200,xmin=0,xmax=48.47/73,linestyle='--')
plt.axvline(x=48.47,ymin=0,ymax=200/500,linestyle='--')
# plt.plot(tp_drs_all/3600,popdyn_mean+popdyn_std,linewidth=2.0,color='blue',label='+/- SD')
# plt.plot(tp_drs_all/3600,popdyn_mean-popdyn_std,linewidth=2.0,color='blue')
plt.xlabel('Time (hours)',fontsize=15)
plt.ylabel('Number of cells',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0,73)
plt.ylim(0,500)
plt.legend(fontsize=15)
plt.show()





#%%

output_files = os.listdir(output_dir)

output_all = {}

for i,file in enumerate(output_files):
    with open(os.path.join(output_dir,file), 'rb') as f:
        output_all['output_g'+str(i+1)] = pickle.load(f)


def load_outputs(dir_dose):
    files = os.listdir(dir_dose)
    output_loaded = {}
    for i,file in enumerate(files):
        with open(os.path.join(dir_dose,file),'rb') as f:
            output_loaded['output_g'+str(i+1)] = pickle.load(f)
    return output_loaded
            
    
output_all = load_outputs(output_dir)

#%%



#%%

for d in range(len(cellpop_sim)):
    plt.plot(timepoints_sim[d]/3600,cellpop_sim[d],label=str(doses_all[d]))



# plt.ylim(0,max(max(cellpop_sim))*1.25)
plt.ylim(0,550)
plt.xlim(0,72)
plt.xlabel('Time (h)')
plt.ylabel('Cell population')
plt.legend(bbox_to_anchor=(1.05,1.0))
plt.show()

#%%

dir1 = os.path.join(output_dir_sim,drug+'_'+str(doses_all[4]))

output1 = load_outputs(dir1)

timecourse_lin(3,drug)



#%% obs parsing

# formula_obs = ["cPARP/PARP","ppERK/(ERK+ppERK)","ppAKT/(AKT**2+ppAKT**3)","pcFos"]
# formula_obs = ["BIM/tBid","BIM/(BIM+pBIM)"]
ObsMat = pd.read_csv(os.path.join(wd,'input_files','Observables.txt'),sep='\t',index_col=0,header=0)

np.array(species_all)[np.where(ObsMat.loc[:,'ERK'].values>0)[0]]


formula_obs_file = 'formula_obs_0.txt'
formula_obs = pd.read_csv(os.path.join(wd,'input_files',formula_obs_file),sep='\t',index_col=None,header=None,squeeze=True)



sp_obs = []
for formula in formula_obs:
    sp_f = re.findall(r'[a-zA-Z]\w*',formula)
    sp_obs.append(sp_f)
sp_obs = [item for sublist in sp_obs for item in sublist]
sp_obs = list(np.unique(sp_obs))

xs = output_all['output_g1']['3']['output']['xoutS']
# xs = rdata_perturb['x']

# xout_sp_obs = {}

for i in range(len(sp_obs)):
    exec(f"{sp_obs[i]} = xs[:,species_all.index('{sp_obs[i]}')]")
    # xout_sp_obs[sp_obs[i]] = xs[:,species_all.index(sp_obs[i])]


obs_cell = pd.Series()



# obs_0 = []
for o,obs in enumerate(formula_obs):
    exec(f"obs_{o} = {obs}")
    exec(f"obs_cell['obs_{o}'] = obs_{o}")
# for i in range(len(sp_obs)):
#     exec(f"{sp_obs[i]} = xs[:,species_all.index('{sp_obs[i]}')]")
#     # exec(f"{sp_obs[i]} = xs[:,species_all.index('{sp_obs[i]}')].reshape((1,len(xs)))",globals(),locals())
    
# obs_k = pd.Series()
# # obs_0 = []
# for o,obs in enumerate(formula_obs):
#     exec(f"obs_{o} = {obs}")
#     exec(f"obs_k['obs_{o}'] = obs_{o}")



#%% Figure - stochastic parameter perturbation

#k1813

# k1813_vals = [0.004,0.0004,0.0003,0.0002]
# fracs_div = [0.125,0.01525,0.0078125,0]

k1813_vals = [0.0004,0.004,0.0003,0.0002]
fracs_div = [0.01525,0.125,0.0078125,0]

plt.scatter(k1813_vals[0],fracs_div[0],c='red',marker='D',s=50,label='Default')
plt.scatter(k1813_vals[1:],fracs_div[1:],c='blue')

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel('Parameter value (basal RasGTP formation)',fontsize=15)
plt.ylabel('Probability of cell division without EGF',fontsize=15)
plt.xscale('log')
plt.xlim(1e-4,1e-2)
plt.legend(bbox_to_anchor=(0.9,0.2))
plt.show()

#%%

#mek activation

multiplier_vals = [1,2,3,5]
frac_div_egf = [0.2734,0.5781,0.65625,0.84375]
frac_div_noegf = [0.0,0.0234,0.039,0.03125]

plt.scatter(multiplier_vals,frac_div_egf,c='red',marker='D',s=50,label='+EGF')
plt.scatter(multiplier_vals,frac_div_noegf,c='blue',label='-EGF')
plt.legend()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
# plt.yscale('log')
plt.xlabel('Rate multiplier',fontsize=15)
plt.ylabel('Probability of cell division',fontsize=15)
plt.show()







#%%