import sys
import os 
import importlib
import jdata as jd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src')

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3J/results')
data = jd.load('Fig3J.json')

observables = []
conditions = []
for condition in data:
    conditions.append(condition)
    for observable in data[condition]['cell 0']:
        observables.append(observable)
print(observables)
fig, ax = plt.subplots(2,2, figsize=(12.5,8))

#Cyclin D plot
ax[0,0].plot(data[conditions[0]]['cell 0'][observables[0]]['toutS']/3600, data[conditions[0]]['cell 0'][observables[0]]['xoutS'], linewidth=4, color='b')
ax[0,0].plot(data[conditions[1]]['cell 0'][observables[0]]['toutS']/3600, data[conditions[1]]['cell 0'][observables[0]]['xoutS'], linewidth=4, color='r')
ax[0,0].plot(data[conditions[2]]['cell 0'][observables[0]]['toutS']/3600, data[conditions[2]]['cell 0'][observables[0]]['xoutS'], linewidth=4, color='tab:orange')
ax[0,0].set_yticklabels(ax[0,0].get_yticks().astype(int), weight='bold', fontsize=16)
ax[0,0].set_xticklabels(ax[0,0].get_xticks().astype(int), weight='bold', fontsize=16)
# ax[0,0].set_title('Cyclin D', weight='bold', fontsize=16)
ax[0,0].set_ylabel('Cyclin D', weight='bold', fontsize=20)

#Cyclin E plot
ax[0,1].plot(data[conditions[0]]['cell 0'][observables[1]]['toutS']/3600, data[conditions[0]]['cell 0'][observables[1]]['xoutS'], linewidth=4, color='b')
ax[0,1].plot(data[conditions[1]]['cell 0'][observables[1]]['toutS']/3600, data[conditions[1]]['cell 0'][observables[1]]['xoutS'], linewidth=4, color='r')
ax[0,1].plot(data[conditions[2]]['cell 0'][observables[1]]['toutS']/3600, data[conditions[2]]['cell 0'][observables[1]]['xoutS'], linewidth=4, color='tab:orange')
ax[0,1].set_yticklabels(ax[0,1].get_yticks().astype(int), weight='bold', fontsize=16)
ax[0,1].set_xticklabels(ax[0,1].get_xticks().astype(int), weight='bold', fontsize=16)
# ax[0,1].set_title('Cyclin E', weight='bold', fontsize=16)
ax[0,1].set_ylabel('Cyclin E', weight='bold', fontsize=20)
ax[0,1].legend(['Basal Cyclin D mRNA x1', 'x10', 'x60'], bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=16)

#Cyclin A plot
ax[1,0].plot(data[conditions[0]]['cell 0'][observables[2]]['toutS']/3600, data[conditions[0]]['cell 0'][observables[2]]['xoutS'], linewidth=4, color='b')
ax[1,0].plot(data[conditions[1]]['cell 0'][observables[2]]['toutS']/3600, data[conditions[1]]['cell 0'][observables[2]]['xoutS'], linewidth=4, color='r')
ax[1,0].plot(data[conditions[2]]['cell 0'][observables[2]]['toutS']/3600, data[conditions[2]]['cell 0'][observables[2]]['xoutS'], linewidth=4, color='tab:orange')
ax[1,0].set_yticklabels(ax[1,0].get_yticks().astype(int), weight='bold', fontsize=16)
ax[1,0].set_xticklabels(ax[1,0].get_xticks().astype(int), weight='bold', fontsize=16)
ax[1,0].set_ylabel('Cyclin A', weight='bold', fontsize=20)
# ax[1,0].set_title('Cyclin A', weight='bold', fontsize=16)

#Cyclin B plot
ax[1,1].plot(data[conditions[0]]['cell 0'][observables[3]]['toutS']/3600, data[conditions[0]]['cell 0'][observables[3]]['xoutS'], linewidth=4, color='b')
ax[1,1].plot(data[conditions[1]]['cell 0'][observables[3]]['toutS']/3600, data[conditions[1]]['cell 0'][observables[3]]['xoutS'], linewidth=4, color='r')
ax[1,1].plot(data[conditions[2]]['cell 0'][observables[3]]['toutS']/3600, data[conditions[2]]['cell 0'][observables[3]]['xoutS'], linewidth=4, color='tab:orange')
ax[1,1].set_yticklabels(ax[1,1].get_yticks().astype(int), weight='bold', fontsize=16)
ax[1,1].set_xticklabels(ax[1,1].get_xticks().astype(int), weight='bold', fontsize=16)
# ax[1,1].set_title('Cyclin B', weight='bold', fontsize=16)
ax[1,1].set_ylabel('Cyclin B', weight='bold', fontsize=20)

fig.tight_layout()
fig.savefig('Fig3J.png', dpi=300)


# %%



