# %%
import os 
import sys 
import jdata as jd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3F/results')

data = jd.load('Fig3F.json')
data

# %%
def cell_to_condition(results_dict):
    """Converts a cell replicate dictionary from cell:condition:results to condition:cell:results
    results_dict: dictionary of cell replicates
    """
    new_dict = {}
    for cell in results_dict:
        for condition in results_dict[cell]:
            if condition not in new_dict:
                new_dict[condition] = {}
            new_dict[condition][cell] = results_dict[cell][condition]
    return new_dict

def condition_to_cell(results_dict):
    """Converts a condition replicate dictionary from condition:cell:results to cell:condition:results
    results_dict: dictionary of condition replicat
    """
    new_dict = {}
    for condition in results_dict:
        for cell in results_dict[condition]:
            if cell not in new_dict:
                new_dict[cell] = {}
            new_dict[cell][condition] = results_dict[condition][cell]
    return new_dict

# %%
reorg = cell_to_condition(data)

#I need to extract the p53 data for each cell for this particular unit test...

p53_only = {}
for condition in reorg:
    p53_only[condition] = {}
    for cell in reorg[condition]:
        if cell not in p53_only:
            p53_only[condition][cell] = {}
        p53_only[condition][cell]['xoutS'] = reorg[condition][cell]['xoutS'][:, 2]
        p53_only[condition][cell]['toutS'] = reorg[condition][cell]['toutS']

# %%
# fig, axes = plt.subplots(4,1, figsize=(8, 4 * len(p53_only)))
fig, axes = plt.subplots(4,1)
for i, condition in enumerate(p53_only):
    ax = axes[i]
    for cell in p53_only[condition]:
        ax.plot(p53_only[condition][cell]['toutS']/3600, p53_only[condition][cell]['xoutS'], label=cell, linewidth=4)
        if condition == 'DamageDSB_1':
            ax.set_ylim(0, 20) 
        else:
            ax.set_ylim(0, 1500)
        # ax.set_title(condition)
        # ax.legend()
        ax.set_xlim(0, 30)
plt.subplots_adjust(hspace=0.5)
plt.show()

# %%


# %%



