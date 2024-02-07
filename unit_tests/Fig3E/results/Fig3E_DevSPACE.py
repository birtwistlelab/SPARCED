import os 
import sys 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import jdata as jd

sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src/')

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3E/petab_files')
data = jd.load('../results/Fig3E.json')

conds = []
for condition in data.keys():
    conds.append(condition)

# Set up a 1x3 subplot grid
fig, axs = plt.subplots(3, 1, figsize=(6, 6), sharex=True)

# Plot the first set of data
axs[0].plot(data[conds[0]]['cell 0']['p53_active']['toutS']/3600, data[conds[0]]['cell 0']['p53_active']['xoutS'], label='Without Repair', color='blue', linewidth=4)
axs[0].plot(data[conds[1]]['cell 0']['p53_active']['toutS']/3600, data[conds[1]]['cell 0']['p53_active']['xoutS'], label='With Repair', color='red', linewidth=4)
axs[0].grid(True)
axs[0].set_ylim(0, 400)
axs[0].set_yticklabels(axs[0].get_yticks(), weight='bold', fontsize=14)
axs[0].set_xticks(np.arange(0, 24, step=6))
axs[0].legend(frameon=False)

# Plot the second set of data
axs[1].plot(data[conds[2]]['cell 0']['p53_active']['toutS']/3600, data[conds[2]]['cell 0']['p53_active']['xoutS'], label='Without Repair', color='blue', linewidth=4)
axs[1].plot(data[conds[3]]['cell 0']['p53_active']['toutS']/3600, data[conds[3]]['cell 0']['p53_active']['xoutS'], label='With Repair', color='red', linewidth=4)
axs[1].grid(True)
axs[1].set_ylim(0, 1000)
axs[1].set_yticklabels(axs[1].get_yticks(), weight='bold', fontsize=14)
axs[1].set_xticks(np.arange(0, 24, step=6))
# axs[1].legend()

# Plot the third set of data
axs[2].plot(data[conds[4]]['cell 0']['p53_active']['toutS']/3600, data[conds[4]]['cell 0']['p53_active']['xoutS'], label='Without Repair', color='blue', linewidth=4)
axs[2].plot(data[conds[5]]['cell 0']['p53_active']['toutS']/3600, data[conds[5]]['cell 0']['p53_active']['xoutS'], label='With Repair', color='red', linewidth=4)
axs[2].grid(True)
axs[2].set_ylim(0, 1200)
axs[2].set_yticklabels(axs[2].get_yticks(), weight='bold', fontsize=14)
axs[2].set_xticklabels(axs[2].get_xticks(), weight='bold', fontsize=14)


fig.savefig('../results/Fig3E.png')
