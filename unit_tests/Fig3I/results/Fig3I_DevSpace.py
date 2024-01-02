# %%
import sys
import os 
import importlib
import jdata as jd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src')
from petab_file_loader import PEtabFileLoader

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3I/results')
data = jd.load('Fig3I.json')

# %%
for condition in data:
    plt.plot(data[condition]['toutS']/3600, data[condition]['xoutS'][:, 105], label=condition)
    plt.legend()
    plt.xlabel('Time (hr)')
    plt.ylabel('cPARP (nM)')
    plt.xlim(0, 60)


# %%



