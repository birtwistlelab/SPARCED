import os 
import sys 
import pandas as pd
import matplotlib.pyplot as plt
import jdata as jd

sys.path.append('/home/jonah/Desktop/SPARCED/unit_tests/src/')

os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3G/results/')

from petab_file_loader import PEtabFileLoader

import importlib
import numpy as np
os.chdir('/home/jonah/Desktop/SPARCED/unit_tests/Fig3G/scripts/')
yaml_file = '../petab_files/Fig3G.yml'

data = jd.load('Fig3G.json')

def plot_species(yaml_file, species, data):
    # Load the PEtab files
    sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

    # Load the SBML model
    model_name = os.path.basename(sbml_file).split('.')[0]
    sys.path.insert(0, os.path.join(os.getcwd(), model_name))

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()

    solver = model.getSolver()
    solver.setMaxSteps = 1e10

    species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt

    species_initializations = np.array(model_module.getModel().getInitialStates()) # Get the initial states from the model

    index = species_ids.index(species)
    for condition in data: 
        plt.plot(data[condition]['toutS']/3600, data[condition]['xoutS'][:, index], label=condition)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xlabel('Time (hours)')
        plt.ylabel(f'{species} (nM)')
        plt.savefig(f'{species}_Fig3G.png')

plot_species(yaml_file, 'cPARP', data)

step_size = len(data['TRAIL_19p25']['toutS']) // 10
scaled_array = [np.mean(data['TRAIL_19p25']['toutS'][i:i+step_size]) for i in range(0, len(data['TRAIL_19p25']['toutS']), step_size)]

scaled_array2 = []
for i in scaled_array:
    scaled_array2.append(i/3600)

traildoses = np.array([0.000385, 0.001925, 0.00385, 0.01925, 0.0385, 0.1925, 0.385, 1.9250, 3.85, 19.25, 38.5])
dosesngperml = traildoses*2.597402597402597e+01
plt.figure(figsize=(7, 4))
plt.plot(np.log10(dosesngperml),scaled_array2,marker='s',linewidth=4)
plt.xlabel('log10(TRAIL dose)')
plt.ylabel('time (hrs)', multialignment='center')
plt.grid(True)
plt.savefig('cPARP_Fig3G3.png')



def time_to_death(yaml_file, data):   
    dead_simulations = {}
    time_of_death = []
    for condition in data:
        dead_simulation = np.argwhere(data[condition]['xoutS'][-1, 103] < data[condition]['xoutS'][-1, 105])
        if len(dead_simulation) > 0:
            dead_simulations[condition] = data[condition]# I need to grab the timepoints where each simulation died
    if len(dead_simulations) < len(data):
        print('Not all simulations died')
        first_condition = list(data.keys())[0]
        [time_of_death.append(data[first_condition]['toutS'].max()/3600)for alive_remainder in range(len(data) - len(dead_simulations))]

    for dead_condition in dead_simulations:
        #Grab the first instance in xoutS where PARP was less than cPARP
        point_of_death = np.argwhere(dead_simulations[dead_condition]['xoutS'][:, 103] < dead_simulations[dead_condition]['xoutS'][:, 105])[0]
        # I need to grab the timepoints where each simulation died

        time_of_death.append(int(dead_simulations[dead_condition]['toutS'][point_of_death]/3600))

    return time_of_death

time_of_death = time_to_death(yaml_file, data)

# print(time_of_death)
traildoses = np.array([0.000385, 0.001925, 0.00385, 0.01925, 0.0385, 0.1925, 0.385, 1.9250, 3.85, 19.25, 38.5])
dosesngperml = traildoses*2.597402597402597e+01
plt.figure(figsize=(7, 4))
plt.plot(np.log10(dosesngperml),time_of_death,marker='s',linewidth=4, color='red')
plt.xlabel('log10(TRAIL dose)')
plt.ylabel('Time To Death (hrs)', multialignment='center')
plt.title('Fig3G3 Unit Test')
plt.ylim(0, 160)
plt.grid(True)
plt.savefig('cPARP_Fig3G3.png')