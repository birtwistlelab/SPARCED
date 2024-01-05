import os
import sys
import jdata as jd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/FigS3C/scripts')

sys.path.insert(0, '/home/jonah/Desktop/SPARCED/unit_tests/src')

data = jd.load('../results/FigS3C.json')

conditional_strings = ['EGF_', 'H_', 'HGF_', 'P_', 'F_', 'I_', 'INS_']
fig, axes = plt.subplots(7, 2, figsize=(6, 8))  # Adjust figsize as needed

for i, string in enumerate(conditional_strings):
    for j, key in enumerate(data):
        if key.startswith(string):
            axes[i, 0].plot(data[key]['ppERK_total']['toutS'] / 60, data[key]['ppERK_total']['xoutS'], label=string)
            # axes[i, 0].set_title(string)

            axes[i, 1].plot(data[key]['ppAKT_total']['toutS'] / 60, data[key]['ppAKT_total']['xoutS'], label=string)
            # axes[i, 1].set_title(string)
    # axes[i, 0].set_ylabel('ppERK_total xoutS')
    # axes[i, 1].set_ylabel('ppAKT_total xoutS')

    # axes[i, 0].set_xlabel('Time (h)')
    # axes[i, 1].set_xlabel('Time (h)')
    legend_values = [0.001, 0.01, 0.1, 1, 10, 100, 1000]
    axes[0, 1].legend(labels = legend_values, loc='upper left', bbox_to_anchor=(1.5, 1), frameon=False)

plt.tight_layout()
plt.subplots_adjust(hspace=0.6)
plt.show()