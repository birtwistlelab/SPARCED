#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd
import numpy as np
import jdata as jd
import matplotlib.pyplot as plt

os.chdir('/home/jrhuggi/projects/foreman/SPARCED/unit_tests/Fig4A/results/')

data = jd.load('Fig4A.json')

sys.path.append('/home/jrhuggi/projects/foreman/SPARCED/unit_tests/src')

from observable_calc import CellDeathMetrics

death_results = CellDeathMetrics(data, 'cPARP_100plus')
obj = death_results.time_to_death()
dead_cells = {}

print(obj['E_0_INS_0_k316_3'][59:79])

dead_cells['E_0_INS_0_k316_3'] = 0
for cell_death in range(60, len(obj['E_0_INS_0_k316_3'])):
    if obj['E_0_INS_0_k316_3'][cell_death] is not None:
        dead_cells['E_0_INS_0_k316_3'] += 1 

dead_cells['E_20_INS_0p001_k316_3']  = 0
for cell in range(len(obj['E_20_INS_0p001_k316_3'])):
    if obj['E_20_INS_0p001_k316_3'][cell] is not None:
        dead_cells['E_20_INS_0p001_k316_3'] += 1
        
labels = ['No Growth Factors', 'Cells with Growth Factors']

dead_cell_values = np.array(list(dead_cells.values()))
# print(dead_cell_values)
dc_sem = np.sqrt((dead_cell_values*(100-dead_cell_values))/100)

fig = plt.figure(figsize=(3, 5))
plt.bar(dead_cells.keys(), dead_cells.values(), yerr=dc_sem, capsize=5, color='blue')
plt.ylim(0, 100)
plt.yticks(np.arange(0, 150, 50), fontsize=16, weight='bold')
plt.ylabel('Number of Dead Cells', fontsize=12, weight='bold')
plt.xticks(np.arange(2), labels = labels, rotation=45, fontsize=16, weight='bold')
plt.tight_layout()

plt.savefig('Fig4A_Barplot.png', dpi=300)