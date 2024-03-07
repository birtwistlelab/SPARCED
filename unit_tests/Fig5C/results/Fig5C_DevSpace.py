import os
import sys
import jdata as jd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig5C/results')


sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src/')
from observable_calc import *

data = jd.load('../results/Fig5C.json')

dead_cells = CellDeathMetrics(data, 'cPARP_total').death_ratio(percent=True)
dead_cell_values = np.array(list(dead_cells.values()))
conditions = []
for condition in data:
    conditions.append(condition)

labels=['BAD-dependent', 'BIM-dependent']
print(dead_cell_values)
num_cells = len(data[conditions[0]])
dc_sem = np.sqrt((dead_cell_values*(100-dead_cell_values))/num_cells)

fig = plt.figure(figsize=(3, 4))
plt.bar(dead_cells.keys(), dead_cells.values(), yerr=dc_sem, capsize=5, color='blue')
plt.ylim(0, 100)
plt.yticks(np.arange(0, 150, 50), fontsize=16, weight='bold')
plt.ylabel('% Death', fontsize=16, weight='bold')
plt.xticks(ticks=range(len(dead_cells)), labels=labels,rotation=45, fontsize=16, weight='bold')
# plt.title('E + I @ 48 hours (SPARCED)', y=1.05)
plt.tight_layout()

plt.savefig('Fig5C.png', dpi=300)