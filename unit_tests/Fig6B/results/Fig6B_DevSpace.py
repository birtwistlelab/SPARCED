import os
import sys
import jdata as jd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig6B/results')

sys.path.insert(0, '/home/jonah/Desktop/SPARCED/unit_tests/src')

data = jd.load('Fig6B.json')


import matplotlib.pyplot as plt

fig, axes = plt.subplots(5, 1, figsize=(7, 7))
for i, condition in enumerate(data):
    for cell in data[condition]:
        axes[i].plot(data[condition][cell]['cycA_Cdk2_total']['toutS']/3600, \
            data[condition][cell]['cycA_Cdk2_total']['xoutS'], linewidth=4)


    # axes[i].set_title(entry)
    axes[i].set_ylim(0, 50)
    axes[i].axvline(x=24, color='black', linestyle='--', label='24-hour point')  # Add vertical dashed line
    axes[i].set_xticklabels(axes[i].get_xticks(), fontsize=16, weight='bold')
    axes[i].set_yticklabels(axes[i].get_yticks(), fontsize=16, weight='bold')

plt.tight_layout()
# plt.subplots_adjust(hspace=-0.00005)
fig.savefig('Fig6B.png', dpi=300)
# Create the proliferation bar plots
num_cells = 30
CellsInSphase = {}
labels = ['Sim', 'Exp']
exp_data = [36.4773, 3.899, 49.2645, 1.55545, 0.81705]
cond_labels = ['EGF', 'INS', 'EGF+INS', 'EGF+INS+MEKi', 'EGF+INS+AKTi']
fig, axes = plt.subplots(5, 1, figsize=(2, 7))
for i, condition in enumerate(data):
    CellsInSphase[condition] = {}
    CellsInSphase[condition]['Ratio'] = []
    for cell in data[condition]:

        if any (value > 20 for value in data[condition][cell]['Sphase']['xoutS']):
            # print(cell)
            CellsInSphase[condition]['Ratio'].append(cell)
    # print(CellsInSphase[condition]['Ratio'])
        
    
    ratio = (len(CellsInSphase[condition]['Ratio'])/num_cells)*100
    CellsInSphase[condition]['Ratio'] = []
    CellsInSphase[condition]['Ratio'].append(ratio)
    CellsInSphase[condition]['Ratio'].append(exp_data[i])
    dc_sem = np.sqrt(ratio*(100-ratio)/num_cells)
    axes[i].bar(labels, CellsInSphase[condition]['Ratio'], yerr=dc_sem, color = ['blue', 'grey'])
    axes[i].set_ylim(0, 60)
    axes[i].set_yticklabels(axes[i].get_yticks(), fontsize=16, weight='bold')
    axes[i].set_ylabel(cond_labels[i], weight='bold')
    axes[i].set_xticks([])
plt.tight_layout()
fig.savefig('SPhase_barplot.png')