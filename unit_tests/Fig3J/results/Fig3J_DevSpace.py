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

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3J/scripts')
data = jd.load('../results/Fig3J.json')


yaml_file = '../petab_files/Fig3J.yml'
# Load the PEtab files
sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

# Load the SBML model
model_name = os.path.basename(sbml_file).split('.')[0]
sys.path.insert(0, os.path.join(os.getcwd(), model_name))

model_module = importlib.import_module(model_name)
model = model_module.getModel()

species_ids = list(model.getStateIds())


# %%
observables = ['Cd', 'Ce', 'Ca', 'Cb']
for observable in observables:
    for condition in data:
        print(data[condition]['xoutS'][:, species_ids.index(observable)])
        plt.plot(data[condition]['toutS'], data[condition]['xoutS'][:, species_ids.index(observable)], label=condition)
        plt.ylabel(observable)
        plt.xlabel('time')
        # plt.ylim(0, 1.1*np.max(data[condition]['xoutS']))
        plt.ticklabel_format(style='plain', axis='y')
        plt.legend()
    plt.show()

# %%
import matplotlib.pyplot as plt

observables = ['Cd', 'Ce', 'Ca', 'Cb']

fig, axes = plt.subplots(2, 2, figsize=(10, 8))

for i, observable in enumerate(observables):
    for condition in data:
        axes[i // 2, i % 2].plot(data[condition]['toutS'], data[condition]['xoutS'][:, species_ids.index(observable)], label=condition)

    axes[i // 2, i % 2].set_ylabel(observable)
    axes[i // 2, i % 2].set_xlabel('time')
    axes[i // 2, i % 2].ticklabel_format(style='plain', axis='y')
    axes[i // 2, i % 2].legend()

plt.tight_layout()
plt.show()


# %%



