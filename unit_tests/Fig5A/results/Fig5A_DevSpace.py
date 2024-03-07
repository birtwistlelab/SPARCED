import os 
import sys
import jdata as jd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig5A/results')

data = jd.load('Fig5A.json')
sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src')
colors = ['blue', 'red', 'yellow', 'purple', 'green']

from observable_calc import CellDeathMetrics

fig, axes = plt.subplots(5, 1, figsize=(7, 7))
for i, condition in enumerate(data):
    for cell in data[condition]:
        axes[i].plot(data[condition][cell]['cPARP_total']['toutS']/3600, \
            data[condition][cell]['cPARP_total']['xoutS'], linewidth=4, color=colors[i])
        
    axes[i].set_ylim(0, 400)
    axes[i].set_xlim(0, 72)
    # axes[0].axvline(x=24, color='black', linestyle='--', label='24-hour point')  # Add vertical dashed line
    axes[i].set_xticklabels(axes[i].get_xticks(), fontsize=16, weight='bold')
    axes[i].set_yticklabels(axes[i].get_yticks(), fontsize=16, weight='bold')

plt.tight_layout()
fig.savefig('Fig5A.png', dpi=300)

dead_cells = CellDeathMetrics(data , 'cPARP_total').time_to_death()
print(dead_cells)
dead_cells_24to72 = {}
fig , axes = plt.subplots(5, 1, figsize=(2.5, 7))
for i, condition in enumerate(dead_cells):
    dead_cells_24to72[condition] = {}
    dead_cells_24to72[condition]['24'] = []
    dead_cells_24to72[condition]['48'] = []
    dead_cells_24to72[condition]['72'] = []
    for cell in range(30):
        if dead_cells[condition][cell] is None:
            dead_cells[condition][cell] = 100
        if dead_cells[condition][cell] <= 24:
            dead_cells_24to72[condition]['24'].append(1)
        if dead_cells[condition][cell] <= 48:
            dead_cells_24to72[condition]['48'].append(1)
        if dead_cells[condition][cell] <= 72:
            dead_cells_24to72[condition]['72'].append(1)
    print(dead_cells_24to72[condition])
    dead_cells_24to72[condition]['24'] = len(dead_cells_24to72[condition]['24'])/len(data[condition])*100
    dead_cells_24to72[condition]['48'] = len(dead_cells_24to72[condition]['48'])/len(data[condition])*100
    dead_cells_24to72[condition]['72'] = len(dead_cells_24to72[condition]['72'])/len(data[condition])*100
    dc_sem = np.sqrt(dead_cells_24to72[condition]['72']*(100-dead_cells_24to72[condition]['72'])/len(data[condition]))
    axes[i].bar(dead_cells_24to72[condition].keys(), dead_cells_24to72[condition].values(), yerr=dc_sem, color=colors[i])
    axes[i].set_ylim(0, 100)
    axes[i].set_yticklabels(axes[i].get_yticks(), fontsize=16, weight='bold')
    # axes[i].set_xlabel(condition, weight='bold')
    if i < 4:
        axes[i].set_xticks([])
    else:
        axes[4].set_xticklabels([24, 48, 72], fontsize=16, weight='bold')

plt.tight_layout()
fig.savefig('Fig5A_barplot.png')