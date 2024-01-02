#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys 
import jdata as jd

os.chdir('/home/jrhuggi/projects/foreman/SPARCED/unit_tests/Fig3H/results/')

data = jd.load('Fig3H.json')


# In[2]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[3]:


__file__ = 'Untitled1.ipynb'
# os.chdir('/home/jrhuggi/projects/foreman/SPARCED/unit_tests/src/scripts/')
sys.path.append('/home/jrhuggi/projects/foreman/SPARCED/unit_tests/src/')
from plotting_unit_tests import PlottingUnitTests
from observable_calc import ObservableCalculator



# In[5]:


reorganized_stoch_cellPop = {}
for cell in data:
    for condition in data[cell]:
        reorganized_stoch_cellPop[condition] = {}
        
for cell in data:
    for condition in data[cell]:
        reorganized_stoch_cellPop[condition][cell] = data[cell][condition]
        


# In[39]:


tod = ObservableCalculator.collect_the_dead(data) #gather dead cells from experiment

percent_dead_cells = []

for condition in tod:
    # print(F"{condition}: {len(tod[condition])}")
    ratio = (len(tod[condition]) / len(data)) * 100
    # print(ratio)
    percent_dead_cells.append(ratio)

doses = [38.5, 3.85, 1.925, 0.0385, 0.1925, 19.25,   0.385]

percent_alive_cells = 100 - np.array(percent_dead_cells)

concatenated_alive_doses = []
concatenated_alive_doses.append(percent_alive_cells)
concatenated_alive_doses.append(doses)
# concatenated_alive_doses = np.array(concatenated_alive_doses)
combined_data = list(zip(doses, percent_alive_cells))

sorted_data = sorted(combined_data, key=lambda x: x[0])

sorted_doses, sorted_percent_alive_cells = zip(*sorted_data)


# In[38]:


dosesngperml = np.array(sorted_doses)*2.597402597402597e+01
print(np.log10(np.array(dosesngperml)))
trailexpdoses = [2.34375,4.6875,9.375,18.75,37.5,75.0,150.0,300.0] # From Bouhaddou2018 model paper
trailexppdeath = [89.79,93.53,88.40,79.53,58.08,36.34,20.18,21.11] # From Bouhaddou2018 model paper

plt.figure(figsize=(7, 4))
plt.plot(np.log10(np.array(trailexpdoses)),np.array(trailexppdeath),marker='*') # Experimental
plt.plot(np.log10(np.array(dosesngperml)),sorted_percent_alive_cells,marker='o') # Simulations
plt.xlabel('TRAIL (ng/mL)', multialignment='center')
plt.ylabel('% Surviving cells')
plt.grid(True)
# plt.xlim(0, 72)
plt.ylim(0, 120)
plt.xticks(np.arange(0,4,step=1),('10e0', '10e1', '10e2', '10e3'))
plt.legend(['Experimental','Simulations'])
plt.savefig('cPARP_Fig3H.png')


# In[ ]:




