# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 22:58:16 2023

@author: Arnab
"""
# import required libraries
import pickle
import copy
import re
import os
import sys
import libsbml
import numpy as np
from Bio import Phylo
from io import StringIO
from scipy.interpolate import interp1d
import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.sans-serif'] = ['Arial']

# Define the working and current directories
cd = os.getcwd()
wd = os.path.dirname(cd)
sys.path.append(os.path.join(wd,'bin'))

sbml_file = "SPARCED.xml" # SBML model name


sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(os.path.join(wd,sbml_file))
sbml_model = sbml_doc.getModel()

species_all = [str(x.getId()) for x in list(sbml_model.getListOfSpecies())]

# define class for reading dose response outputs
class drs_dict(): ###JRHUGGI: drs is a little obfusicated in meaning, I suggest something easier to understand.###
    def __init__(self,main,drug,rep,dose_level): 
        self.main = main
        self.drug = drug
        self.rep = rep
        self.dose_level = dose_level
        self.path_exp = os.path.join(main,'drs_'+drug)
        self.path_reps = os.listdir(self.path_exp) #lists drug directory within main
        self.reps = [int(list(filter(str.isdigit, x.split('_')[-1]))[0]) for x in self.path_reps]
        self.path_rep = os.path.join(self.path_exp,self.path_exp,'drs_'+str(drug)+'_rep'+str(rep))
        self.path_doses = os.listdir(self.path_rep)
        self.doses = [float(x.split('_')[-1]) for x in self.path_doses]
        self.doses.sort()
        ###JRHUGGI: The above lines are convoluted as hell. I suggest you use a dictionary to store all information, such as a JSON or XML dictionary file.###
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

        obs_all = []
        touts = []
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
