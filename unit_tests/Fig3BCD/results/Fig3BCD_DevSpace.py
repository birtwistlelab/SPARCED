import os 
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import jdata as jd

import re

sys.path.insert(0, '/home/jonah/Desktop/SPARCED/unit_tests/src')

from observable_calc import ObservableCalculator

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3BCD/results/')

data = jd.load('Fig3BCD.json')

variable = ObservableCalculator.experimental_comparator('../petab_files/Fig3BCD.yml')

# observables = ['ppERK_total', 'ppAKT_total', 'pEIF4BP1_total']
observables = set([observable for condition in data.values() for observable in condition['cell 0'].keys()])
for observable in observables:
    for condition in data:
        if re.search(r'^EGF_[0-9]', condition):
            plt.plot(data[condition]['cell 0'][observable]['toutS']/60,data[condition]['cell 0'][observable]['xoutS'], label=condition, linewidth=4)
            plt.grid(True)
            plt.xlim(0, 360)
            plt.xticks(np.arange(0,360,step=60))
            plt.xlabel('Time (min)')
            if observable == 'ppERK_total':
                plt.ylim(0, 120)
                plt.yticks(np.arange(0,120,step=20))

            if observable == 'ppAKT_total':
                plt.ylim(0, 40)
                plt.yticks(np.arange(0,120,step=20))

            if observable == 'pEIF4BP1_total':
                plt.ylim(0, 10)
                plt.yticks(np.arange(0,120,step=20))                
            plt.ylabel(f'{observable} (nM)')
            plt.legend()
    plt.show()

    for condition in data:
        if re.search(r'^INS_[0-9]', condition):
            plt.plot(data[condition]['cell 0'][observable]['toutS']/60,data[condition]['cell 0'][observable]['xoutS'], label=condition, linewidth=4)
            plt.grid(True)
            plt.xlim(0, 360)
            plt.xticks(np.arange(0,360,step=60))
            plt.xlabel('Time (min)')
            plt.ylabel(f'{observable} (nM)')
            if observable == 'ppERK_total':
                plt.ylim(0, 120)
                plt.yticks(np.arange(0,120,step=20))

            if observable == 'ppAKT_total':
                plt.ylim(0, 40)
                plt.yticks(np.arange(0,120,step=20))

            if observable == 'pEIF4BP1_total':
                plt.ylim(0, 10)
            plt.legend()
    plt.show()

    for condition in data:
        if re.search(r'^EGF_INS', condition):
            plt.plot(data[condition]['cell 0'][observable]['toutS']/60,data[condition]['cell 0'][observable]['xoutS'], label=condition, linewidth=4)
            plt.grid(True)
            plt.xlim(0, 360)
            plt.xticks(np.arange(0,360,step=60)) 
            plt.xlabel('Time (min)')
            plt.ylabel(f'{observable} (nM)')
            if observable == 'ppERK_total':
                plt.ylim(0, 120)
                plt.yticks(np.arange(0,120,step=20))

            if observable == 'ppAKT_total':
                plt.ylim(0, 40)
                plt.yticks(np.arange(0,120,step=20))

            if observable == 'pEIF4BP1_total':
                plt.ylim(0, 10)
            plt.legend()
    plt.show()


observables = set([observable for condition in variable.values() for observable in condition.keys()])
for observable in observables:
    for condition in variable:
        if re.search(r'^EGF_[0-9]', condition):
            plt.plot(variable[condition][observable]['time']/60,variable[condition][observable]['measurement'], label=condition, linewidth=4)
            plt.grid(True)
            plt.xlim(0, 360)
            plt.xticks(np.arange(0,360,step=60))
            plt.xlabel('Time (min)')
            plt.ylabel(f'{observable} (AU)')
            if observable == 'ppERK_total':
                plt.ylim(0, 5)
                plt.yticks(np.arange(0,5,step=5))

            if observable == 'ppAKT_total':
                plt.ylim(0, 50)
                plt.yticks(np.arange(0,50,step=50))

            if observable == 'pEIF4BP1_total':
                plt.ylim(0, 4)
            plt.legend()
    plt.show()

    for condition in variable:
        if re.search(r'^INS_[0-9]', condition):
            plt.plot(variable[condition][observable]['time']/60,variable[condition][observable]['measurement'], label=condition, linewidth=4)
            plt.grid(True)
            plt.xlim(0, 360)
            plt.xticks(np.arange(0,360,step=60))
            plt.xlabel('Time (min)')
            plt.ylabel(f'{observable} (AU)')
            if observable == 'ppERK_total':
                plt.ylim(0, 5)
                plt.yticks(np.arange(0,5,step=5))

            if observable == 'ppAKT_total':
                plt.ylim(0, 50)
                plt.yticks(np.arange(0,50,step=50))

            if observable == 'pEIF4BP1_total':
                plt.ylim(0, 4)
            plt.legend()
    plt.show()

    for condition in variable:
        if re.search(r'^EGF_INS', condition):
            plt.plot(variable[condition][observable]['time']/60,variable[condition][observable]['measurement'], label=condition, linewidth=4)
            plt.grid(True)
            plt.xlim(0, 360)
            plt.xticks(np.arange(0,360,step=60))
            plt.xlabel('Time (min)')
            plt.ylabel(f'{observable} (AU)')
            if observable == 'ppERK_total':
                plt.ylim(0, 5)
                plt.yticks(np.arange(0,5,step=5))

            if observable == 'ppAKT_total':
                plt.ylim(0, 50)
                plt.yticks(np.arange(0,50,step=50))

            if observable == 'pEIF4BP1_total':
                plt.ylim(0, 4)
            plt.legend()
    plt.show()