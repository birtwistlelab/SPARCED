import os
import sys
import importlib
from typing import Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import jdata as jd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Get the directory path
wd = os.path.dirname(os.path.abspath(__file__))
# Create the formatted path
sparced_root = '/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep).index('SPARCED')+1])
sys.path.append(os.path.join(sparced_root, 'unit_tests/src'))
from petab_file_loader import PEtabFileLoader
from sparced_erm import SPARCED_ERM
from observable_calc import ObservableCalculator

class PlottingUnitTests:
    # Loop through the data and store unique conditions
    def plot_data(yaml_file, json_dict, observable: Optional[str] = None):
        # Load PEtab files
        sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

        perturbants = list(conditions_df.columns[2:])
        unique_conditions = conditions_df.drop_duplicates(subset=perturbants)

        # Get the number of unique conditions
        num_conditions = len(unique_conditions)

        # Create a colormap
        color_map = cm.get_cmap('gist_ncar', num_conditions)

        # Load the SBML model
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(os.getcwd(), model_name))
        
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()
        species_ids = list(model.getStateIds())

        if observable is not None:
            index = species_ids.index(perturbants[observable])
        else:
            index = species_ids.index(perturbants[0])
        
        # Plot the data and add a custom legend with unique conditions and colors
        for i, condition in enumerate(unique_conditions['conditionId']):
            perturbation_data = measurement_df[measurement_df['simulationConditionId'] == condition]
            color = color_map(i)
            plt.plot(perturbation_data['time'], perturbation_data['measurement'], label=condition, color=color)

        plt.legend()
        plt.xlabel('Time')
        plt.ylabel('Measurement')
        plt.title('experimental data')
        plt.show()  # Show the plot

        # for i, condition in enumerate(unique_conditions['conditionId']):
        data = jd.load(json_dict)
        for i, condition in enumerate(unique_conditions['conditionId']):
            color = color_map(i)
            plt.plot(data[condition]['toutS'], data['EGF_INS_010_1721p0']['xoutS'][:, index], label=condition, color = color)
        plt.legend()
        plt.xlabel('Time')
        plt.ylabel('Measurement')
        plt.title('simulation data')
        plt.show()