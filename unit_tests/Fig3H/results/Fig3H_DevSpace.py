#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys 
import jdata as jd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src/')

from observable_calc import *

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3H/results/')
data = jd.load('Fig3H.json')

doses = [0.0385, 0.1925, 0.385, 1.925, 3.85, 19.25, 38.5]
alive_ratio = CellDeathMetrics(data, 'cPARP_total').alive_ratio()
dosesngperml = np.array(doses)*2.597402597402597e+01
# print(alive
trailexpdoses = [2.34375,4.6875,9.375,18.75,37.5,75.0,150.0,300.0] # From Bouhaddou2018 model paper
trailexppdeath = [89.79,93.53,88.40,79.53,58.08,36.34,20.18,21.11] # From Bouhaddou2018 model paper

plt.figure(figsize=(7, 4))
plt.plot(np.log10(np.array(trailexpdoses)),np.array(trailexppdeath),marker='*', color = 'black', linewidth=4, label='Experiment') # Experimental
plt.plot(np.log10(np.array(dosesngperml)),alive_ratio,marker='o', color='red', linewidth=4, label='Simulation') # Simulations
# plt.xlabel('TRAIL (ng/mL)', multialignment='center')
# plt.ylabel('% Surviving cells')
plt.grid(True)
# plt.xlim(0, 72)
plt.ylim(0, 120)
plt.xticks(np.arange(0,4,step=1),('10e0', '10e1', '10e2', '10e3'), weight='bold', fontsize=16)
ax = plt.gca()
ax.set_yticklabels(ax.get_yticks(), weight='bold', fontsize=16)
legend_properties = {'weight':'bold'}
plt.legend(['Experimental','Simulations'], frameon = False, bbox_to_anchor=[1.05, 1], loc='upper left', fontsize=16, prop=legend_properties)
plt.tight_layout()
plt.savefig('cPARP_Fig3H.png')




