import os
import sys
import shutil
import importlib
import pandas as pd
import numpy as np
import amici
import jdata as jd
import argparse
from typing import Optional
from petab_file_loader import PEtabFileLoader
import libsbml
import amici.plotting
import matplotlib.pyplot as plt
from datetime import datetime
import scipy.stats

# Get the directory path
wd = os.path.dirname(os.path.abspath(__file__))

# Create the formatted path
sparced_root = '/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep).index('SPARCED')+1])
sys.path.append(os.path.join(sparced_root, 'bin'))
from modules.RunSPARCED import RunSPARCED


class SPARCED_ERM:
    def __init__(self, yaml_file: str):
        self.yaml_file = yaml_file
        self.results_dict = self.sparced_erm()

    def sparced_erm(yaml_file: str, flagD: Optional[int] = None):
            """Simulate the experimental replicate model."""
            
            current_directory = os.getcwd()

            # Load the PEtab files
            sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

            # Load the SBML model
            model_name = os.path.basename(sbml_file).split('.')[0]
            sys.path.insert(0, os.path.join(current_directory, model_name))

            # Set model mathematical representation
            if  flagD == None:
                flagD = 1

            flagD = int(flagD)


            model_module = importlib.import_module(model_name)
            model = model_module.getModel()

            solver = model.getSolver()
            solver.setMaxSteps = 1e10


            # Create dynamic unit tests based on PEtab files
            results_dict = {}

            perturbants = list(conditions_df.columns[2:])
            unique_conditions = conditions_df.drop_duplicates(subset=perturbants)

            # Timepoints are set by the number of unique timepoints and maximum timepoint in the measurement table
            simulation_time = measurement_df['time'].max()
            # Set the number of records as the number of unique timepoints
            model.setTimepoints(np.linspace(0, simulation_time, len(measurement_df['time'].unique())))

            species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt
            species_initializations = np.array(model_module.getModel().getInitialStates()) # Get the initial states from the model

            for index, condition in unique_conditions.iterrows(): # Iterate through the unique conditions
                species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0 # Set any initializations less than 1e-6 to 0.0
                for species in perturbants:
                    # Set the initial concentrations for the perturbants in the conditions
                    species_initializations[species_ids.index(species)] = condition[species] 


                # Run the simulation
                print(f"Running simulation for condition {condition['conditionId']}")
                xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,simulation_time,species_initializations,[],sbml_file,model)

                
                iteration_name = condition['conditionId']

                #Build the results dictionary
                results_dict[iteration_name] = {}
                results_dict[iteration_name]['xoutS'] = xoutS_all
                results_dict[iteration_name]['toutS'] = tout_all

            return results_dict
