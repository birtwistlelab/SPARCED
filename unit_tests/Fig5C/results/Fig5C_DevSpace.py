import os
import sys
import jdata as jd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig5C/scripts')

data = jd.load('../results/Fig5C.json')

dead_cells = {}
for condition in data:
    dead_cells[condition] = []
    for cell in data[condition]:
        cPARP = data[condition][cell]['cPARP_total']['xoutS']
        PARP = data[condition][cell]['PARP_total']['xoutS']
        if cPARP[-1] > 100:
            dead_cells[condition].append(cell)

ratio ={}
for condition in dead_cells:
    ratio[condition] = (len(dead_cells[condition])/len(data[condition]))*100

ratio2 = {'BIM-dependent': 48.0, 'BAD-dependent':16.0} # I just manually created a variant that shuffled the BIM and BAD locations to match
fig = plt.figure(figsize=(2.5, 4))
plt.bar(ratio2.keys(), ratio2.values())
plt.ylim(0, 100)
plt.yticks(np.arange(0, 150, 50))
plt.ylabel('% Death')
plt.xticks(rotation=45)
plt.title('E + I @ 48 hours (SPARCED)', y=1.05)
plt.tight_layout()
# plt.legend('SIM', frameon=False)
plt.show()