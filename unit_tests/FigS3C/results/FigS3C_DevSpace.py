import sys
import os 
import importlib
import jdata as jd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src')

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/FigS3C/results')
data = jd.load('FigS3C.json')

conditional_strings = ['EGF', 'HRG', 'HGF', 'PDGF', 'FGF', 'IGF', 'INS']
fig, axes = plt.subplots(7, 2, figsize=(8, 12))  # Adjust figsize as needed
conditions =[]
for condition in data:
    conditions.append(condition)

for i in range(7):
    upper = (i+1)*7
    lower = upper - 7
    for j in range(lower, upper):
        axes[i, 0].plot(data[conditions[j]]['cell 0']['ppERK_total']['toutS'] / 60, data[conditions[j]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
        axes[i, 1].plot(data[conditions[j]]['cell 0']['ppAKT_total']['toutS'] / 60, data[conditions[j]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)

    axes[i, 1].set_yticklabels(axes[i, 1].get_yticks().astype(int), fontsize=16, weight='bold')
    axes[i, 1].set_xticklabels(axes[i, 1].get_xticks().astype(int), fontsize=16, weight='bold')
    axes[i, 0].set_yticklabels(axes[i, 0].get_yticks().astype(int), fontsize=16, weight='bold')
    axes[i, 0].set_xticklabels(axes[i, 0].get_xticks().astype(int), fontsize=16, weight='bold')
    axes[i, 1].set_ylabel(conditional_strings[i], fontsize=20, weight='bold', labelpad=-250, rotation=0)
axes[3, 1].set_ylabel(conditional_strings[3], fontsize=20, weight='bold', labelpad=-275, rotation=0)
axes[4, 1].set_ylabel(conditional_strings[4], fontsize=20, weight='bold', labelpad=-268, rotation=0)
axes[6, 1].set_ylabel(conditional_strings[6], fontsize=20, weight='bold', labelpad=-235, rotation=0)

axes[3, 1].set_yticklabels(np.arange(0.0,1.0, step=0.25), fontsize=16, weight='bold')
axes[4, 1].set_yticklabels(axes[4, 1].get_yticks(), fontsize=16, weight='bold')
legend_values = ['0.001 nM', '0.01 nM', '0.1 nM', '1 nM', '10 nM', '100 nM', '1000 nM']
legend_properties = {'weight':'bold'}
axes[0, 1].legend(labels = legend_values, loc='upper left', bbox_to_anchor=(1.5, 1), frameon=False, fontsize=20, prop=legend_properties)


# plt.tight_layout()
plt.subplots_adjust(hspace=0.5, wspace=0.5)
fig.savefig('FigS3C.png', dpi=300, bbox_inches='tight')