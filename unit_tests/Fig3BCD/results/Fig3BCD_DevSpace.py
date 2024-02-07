import os
import sys
import jdata as jd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3BCD/results')

sys.path.insert(0, '/home/jonah/Desktop/SPARCED/unit_tests/src')

from petab_file_loader import PEtabFileLoader
from observable_calc import ObservableCalculator

# Plot experiment beside simulation
data = jd.load('Fig3BCD.json')
variable = ObservableCalculator.experimental_comparator('../petab_files/Fig3BCD.yml')

conditions = []

for condition in data:
    conditions.append(condition)

fig, axs = plt.subplots(9, 2, figsize=(11, 10))

#Experimental EGF for ppERK
axs[0,0].plot(variable[conditions[0]]['ppERK_total']['time']/3600, variable[conditions[0]]['ppERK_total']['measurement'], linewidth=4)
axs[0,0].plot(variable[conditions[1]]['ppERK_total']['time']/3600, variable[conditions[1]]['ppERK_total']['measurement'], linewidth=4)
axs[0,0].plot(variable[conditions[2]]['ppERK_total']['time']/3600, variable[conditions[2]]['ppERK_total']['measurement'], linewidth=4)
axs[0,0].plot(variable[conditions[3]]['ppERK_total']['time']/3600, variable[conditions[3]]['ppERK_total']['measurement'], linewidth=4)

#Experimental INS for ppERK
axs[1,0].plot(variable[conditions[4]]['ppERK_total']['time']/3600, variable[conditions[4]]['ppERK_total']['measurement'], linewidth=4)
axs[1,0].plot(variable[conditions[5]]['ppERK_total']['time']/3600, variable[conditions[5]]['ppERK_total']['measurement'], linewidth=4)
axs[1,0].plot(variable[conditions[6]]['ppERK_total']['time']/3600, variable[conditions[6]]['ppERK_total']['measurement'], linewidth=4)
axs[1,0].plot(variable[conditions[7]]['ppERK_total']['time']/3600, variable[conditions[7]]['ppERK_total']['measurement'], linewidth=4)

#Experimental EGF+INS for ppERK
axs[2,0].plot(variable[conditions[8]]['ppERK_total']['time']/3600, variable[conditions[8]]['ppERK_total']['measurement'],linewidth=4)
axs[2,0].plot(variable[conditions[9]]['ppERK_total']['time']/3600, variable[conditions[9]]['ppERK_total']['measurement'], linewidth=4)
axs[2,0].plot(variable[conditions[10]]['ppERK_total']['time']/3600, variable[conditions[10]]['ppERK_total']['measurement'], linewidth=4)
axs[2,0].plot(variable[conditions[11]]['ppERK_total']['time']/3600, variable[conditions[11]]['ppERK_total']['measurement'], linewidth=4)

#Experimental EGF for ppAKT
axs[3,0].plot(variable[conditions[0]]['ppAKT_total']['time']/3600, variable[conditions[0]]['ppAKT_total']['measurement'], linewidth=4)
axs[3,0].plot(variable[conditions[1]]['ppAKT_total']['time']/3600, variable[conditions[1]]['ppAKT_total']['measurement'], linewidth=4)
axs[3,0].plot(variable[conditions[2]]['ppAKT_total']['time']/3600, variable[conditions[2]]['ppAKT_total']['measurement'], linewidth=4)
axs[3,0].plot(variable[conditions[3]]['ppAKT_total']['time']/3600, variable[conditions[3]]['ppAKT_total']['measurement'], linewidth=4)

#Experimental INS for ppAKT
axs[4,0].plot(variable[conditions[4]]['ppAKT_total']['time']/3600, variable[conditions[4]]['ppAKT_total']['measurement'], linewidth=4)
axs[4,0].plot(variable[conditions[5]]['ppAKT_total']['time']/3600, variable[conditions[5]]['ppAKT_total']['measurement'], linewidth=4)
axs[4,0].plot(variable[conditions[6]]['ppAKT_total']['time']/3600, variable[conditions[6]]['ppAKT_total']['measurement'], linewidth=4)
axs[4,0].plot(variable[conditions[7]]['ppAKT_total']['time']/3600, variable[conditions[7]]['ppAKT_total']['measurement'], linewidth=4)

#Experimental EGF+INS for ppAKT
axs[5,0].plot(variable[conditions[8]]['ppAKT_total']['time']/3600, variable[conditions[8]]['ppAKT_total']['measurement'], linewidth=4)
axs[5,0].plot(variable[conditions[9]]['ppAKT_total']['time']/3600, variable[conditions[9]]['ppAKT_total']['measurement'], linewidth=4)
axs[5,0].plot(variable[conditions[10]]['ppAKT_total']['time']/3600, variable[conditions[10]]['ppAKT_total']['measurement'], linewidth=4)
axs[5,0].plot(variable[conditions[11]]['ppAKT_total']['time']/3600, variable[conditions[11]]['ppAKT_total']['measurement'], linewidth=4)

#Experimental EGF for pEIF4BP1
axs[6,0].plot(variable[conditions[0]]['pEIF4BP1_total']['time']/3600, variable[conditions[0]]['pEIF4BP1_total']['measurement'], linewidth=4)
axs[6,0].plot(variable[conditions[1]]['pEIF4BP1_total']['time']/3600, variable[conditions[1]]['pEIF4BP1_total']['measurement'], linewidth=4)
axs[6,0].plot(variable[conditions[2]]['pEIF4BP1_total']['time']/3600, variable[conditions[2]]['pEIF4BP1_total']['measurement'], linewidth=4)
axs[6,0].plot(variable[conditions[3]]['pEIF4BP1_total']['time']/3600, variable[conditions[3]]['pEIF4BP1_total']['measurement'], linewidth=4)

#Experimental INS for pEIF4BP1
axs[7,0].plot(variable[conditions[4]]['pEIF4BP1_total']['time']/3600, variable[conditions[4]]['pEIF4BP1_total']['measurement'], linewidth=4)
axs[7,0].plot(variable[conditions[5]]['pEIF4BP1_total']['time']/3600, variable[conditions[5]]['pEIF4BP1_total']['measurement'], linewidth=4)
axs[7,0].plot(variable[conditions[6]]['pEIF4BP1_total']['time']/3600, variable[conditions[6]]['pEIF4BP1_total']['measurement'], linewidth=4)
axs[7,0].plot(variable[conditions[7]]['pEIF4BP1_total']['time']/3600, variable[conditions[7]]['pEIF4BP1_total']['measurement'], linewidth=4)

#Experimental EGF+INS for pEIF4BP1
axs[8,0].plot(variable[conditions[8]]['pEIF4BP1_total']['time']/60, variable[conditions[8]]['pEIF4BP1_total']['measurement'], linewidth=4)
axs[8,0].plot(variable[conditions[9]]['pEIF4BP1_total']['time']/60, variable[conditions[9]]['pEIF4BP1_total']['measurement'], linewidth=4)
axs[8,0].plot(variable[conditions[10]]['pEIF4BP1_total']['time']/60, variable[conditions[10]]['pEIF4BP1_total']['measurement'], linewidth=4)
axs[8,0].plot(variable[conditions[11]]['pEIF4BP1_total']['time']/60, variable[conditions[11]]['pEIF4BP1_total']['measurement'], linewidth=4)

#Simulated EGF for ppERK
axs[0,1].plot(data[conditions[0]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[0]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
axs[0,1].plot(data[conditions[1]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[1]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
axs[0,1].plot(data[conditions[2]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[2]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
axs[0,1].plot(data[conditions[3]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[3]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)

#Simulated INS for ppERK
axs[1,1].plot(data[conditions[4]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[4]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
axs[1,1].plot(data[conditions[5]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[5]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
axs[1,1].plot(data[conditions[6]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[6]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
axs[1,1].plot(data[conditions[7]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[7]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)

#Simulated EGF+INS for ppERK
axs[2,1].plot(data[conditions[8]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[8]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
axs[2,1].plot(data[conditions[9]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[9]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
axs[2,1].plot(data[conditions[10]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[10]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)
axs[2,1].plot(data[conditions[11]]['cell 0']['ppERK_total']['toutS']/3600, data[conditions[11]]['cell 0']['ppERK_total']['xoutS'], linewidth=4)

#Simulated EGF for ppAKT
axs[3,1].plot(data[conditions[0]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[0]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)
axs[3,1].plot(data[conditions[1]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[1]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)
axs[3,1].plot(data[conditions[2]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[2]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)
axs[3,1].plot(data[conditions[3]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[3]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)

#Simulated INS for ppAKT
axs[4,1].plot(data[conditions[4]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[4]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)
axs[4,1].plot(data[conditions[5]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[5]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)
axs[4,1].plot(data[conditions[6]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[6]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)
axs[4,1].plot(data[conditions[7]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[7]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)

#Simulated EGF+INS for ppAKT
axs[5,1].plot(data[conditions[8]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[8]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)
axs[5,1].plot(data[conditions[9]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[9]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)
axs[5,1].plot(data[conditions[10]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[10]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)
axs[5,1].plot(data[conditions[11]]['cell 0']['ppAKT_total']['toutS']/3600, data[conditions[11]]['cell 0']['ppAKT_total']['xoutS'], linewidth=4)

#Simulated EGF for pEIF4BP1
axs[6,1].plot(data[conditions[0]]['cell 0']['pEIF4BP1_total']['toutS']/3600, data[conditions[0]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)
axs[6,1].plot(data[conditions[1]]['cell 0']['pEIF4BP1_total']['toutS']/3600, data[conditions[1]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)
axs[6,1].plot(data[conditions[2]]['cell 0']['pEIF4BP1_total']['toutS']/3600, data[conditions[2]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)
axs[6,1].plot(data[conditions[3]]['cell 0']['pEIF4BP1_total']['toutS']/3600, data[conditions[3]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)

#Simulated INS for pEIF4BP1
axs[7,1].plot(data[conditions[4]]['cell 0']['pEIF4BP1_total']['toutS']/3600, data[conditions[4]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)
axs[7,1].plot(data[conditions[5]]['cell 0']['pEIF4BP1_total']['toutS']/3600, data[conditions[5]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)
axs[7,1].plot(data[conditions[6]]['cell 0']['pEIF4BP1_total']['toutS']/3600, data[conditions[6]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)
axs[7,1].plot(data[conditions[7]]['cell 0']['pEIF4BP1_total']['toutS']/3600, data[conditions[7]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)

#Simulated EGF+INS for pEIF4BP1
axs[8,1].plot(data[conditions[8]]['cell 0']['pEIF4BP1_total']['toutS']/60, data[conditions[8]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)
axs[8,1].plot(data[conditions[9]]['cell 0']['pEIF4BP1_total']['toutS']/60, data[conditions[9]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)
axs[8,1].plot(data[conditions[10]]['cell 0']['pEIF4BP1_total']['toutS']/60, data[conditions[10]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)
axs[8,1].plot(data[conditions[11]]['cell 0']['pEIF4BP1_total']['toutS']/60, data[conditions[11]]['cell 0']['pEIF4BP1_total']['xoutS'], linewidth=4)


num_of_ax = [0,1,2,3,4,5,6,7,8]
for i in num_of_ax:
    axs[i, 0].set_yticklabels(axs[i, 0].get_yticks(), weight='bold', fontsize=14)
    axs[i, 1].set_yticklabels(axs[i, 1].get_yticks(), weight='bold', fontsize=14)

    axs[i, 0].set_xticks([])
    axs[i, 1].set_xticks([])    

axs[8,0].set_xticks([0, 60, 120, 240, 360])
axs[8,0].set_xticklabels([0, 60, 120, 240, 360], weight='bold', fontsize=14)

axs[8,1].set_xticks([0, 60, 120, 240, 360])
axs[8,1].set_xticklabels([0, 60, 120, 240, 360], weight='bold', fontsize=14)

axs[0,1].legend(['E = low','E = med','E = high','E = highest'], bbox_to_anchor=(1.46, 1.2), loc='upper right', frameon=False)
axs[1,1].legend(['I = low','I = med','I = high','I = highest'], bbox_to_anchor=(1.45, 1.2), loc='upper right', frameon=False)
axs[2,1].legend(['E-INS = low-low','E-INS = low-highest','E-INS = highest-low','E-INS = highest-highest'], bbox_to_anchor=(1.7, 1.2), loc='upper right', frameon=False)


fig.savefig('Fig3BCD.png', bbox_inches='tight', dpi=300)
