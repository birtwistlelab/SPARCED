import sys
import os 
import importlib
import jdata as jd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src')

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3I/results')
data = jd.load('Fig3I.json')

colors = ['k', 'b', 'r']
plt.figure(figsize=(7, 4))
for i, condition in enumerate(data):
    plt.plot(data[condition]['cell 0']['cPARP_total']['toutS']/3600, data[condition]['cell 0']['cPARP_total']['xoutS'], label=condition, linewidth=4, color=colors[i])
    # plt.xlabel('Time (hr)')
    # plt.ylabel('cPARP (nM)')
    plt.grid(True)
    plt.xlim(0, 50)
    plt.ylim(0, 400)
    plt.xticks(fontsize=16, weight='bold')
    plt.yticks(fontsize=16, weight='bold')
    plt.legend(['Low TRAIL dose','+ppERK / +ppAKT','+PUMA / +NOXA'], bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    plt.tight_layout()
    plt.savefig('Fig3I.png', dpi=300)



