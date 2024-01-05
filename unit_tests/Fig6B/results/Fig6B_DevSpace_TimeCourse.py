import os
import sys
import jdata as jd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig6B/scripts')

sys.path.insert(0, '/home/jonah/Desktop/SPARCED/unit_tests/src')

data = jd.load('../results/Fig6B.json')


import matplotlib.pyplot as plt

# Assuming 'data' is your dictionary containing the entries
# 'data' should have the structure: {entry1: {cell1: {...}, cell2: {...}, ...}, entry2: {...}, ...}

# Create subplots
fig, axes = plt.subplots(5, 1, figsize=(7, 7))

# Iterate over dictionary entries
for i, (entry, cell_data) in enumerate(data.items()):
    # Iterate over cells in each entry
    for cell, cell_values in cell_data.items():
        # Assuming 'toutS' and 'xoutS' are keys in the cell_values dictionary
        axes[i].plot(cell_values['cycA_Cdk2_total']['toutS']/3600, cell_values['cycA_Cdk2_total']['xoutS'], label=f'Cell {cell}')
        # if i == 4:
            # axes[i].set_xlabel('Time (h)', fontsize=18)

    # axes[i].set_title(entry)
    axes[i].set_ylim(0, 50)
    axes[i].axvline(x=24, color='black', linestyle='--', label='24-hour point')  # Add vertical dashed line
    # axes[i].set_xlabel('Time')
    # axes[i].set_ylabel('CycA/Cdk2 (nM)')
    # axes[i].legend()

plt.tight_layout()
# plt.subplots_adjust(hspace=-0.00005)
plt.show()


# Create the proliferation bar plots
num_cells = 30
CellsInSphase = {}
labels = ['Sim', 'Exp']
exp_data = [36.4773, 3.899, 49.2645, 1.55545, 0.81705]
fig, axes = plt.subplots(5, 1, figsize=(1.5, 7))
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

    axes[i].bar(labels, CellsInSphase[condition]['Ratio'], label=condition, color = ['blue', 'grey'])
    axes[i].set_ylim(0, 60)
    # axes[i].set_ylabel('Ratio of cells in S phase')
    # axes[i].set_xlabel('Sim vs. Exp')
    axes[i].set_xticks([])
    # axes[i].legend()
plt.tight_layout()