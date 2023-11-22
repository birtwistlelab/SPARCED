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

# Ensure the SPARCED root and bin directories are in the system path
sparced_root = '/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep).index('SPARCED')+1])
sys.path.append(os.path.join(sparced_root, 'bin'))
from modules.RunSPARCED import RunSPARCED


class SPARCED_ERM:
    def __init__(self, yaml_file: str):
        self.yaml_file = yaml_file
        self.results_dict = self.sparced_erm()

    def preincubate(yaml_file, flagP: Optional[int] = None):
        """Simulate the preincubation step."""
        if flagP != None:
            th = int(flagP)
            ts = 30
            # Load the PEtab files
            sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

            # Load the SBML model
            model_name = os.path.basename(sbml_file).split('.')[0]
            sys.path.insert(0, os.path.join(os.getcwd(), model_name))

            # Here, we import the model's internal packages as a python module
            model_module = importlib.import_module(model_name)
            model = model_module.getModel()

            solver = model.getSolver()
            solver.setMaxSteps = 1e10

            model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

            # Serum starve the cell by setting the initial conditions of species to 0.0
            species_initializations = np.array(model_module.getModel().getInitialStates())
            species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0 

            # Run SPARCED for preincubation time (th) with stimulus concentrations set to 0
            xoutS_all, xoutG_all, tout_all = RunSPARCED(0, th,species_initializations,[],sbml_file,model)

            # Store the preincubation results in a dictionary
            return xoutS_all[-1]


    def delayed_stimulus(xoutS_all, yaml_file: str, flagD: Optional[int] = None, flagP: Optional[int] = None):
        # Load the PEtab files
        sbml_file, _, _, _, _ = PEtabFileLoader.load_petab_files(yaml_file)

        # Load the SBML model
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(os.getcwd(), model_name))
        
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt


        #Find out if SPARCED died during the first round of stimulus
        death_point = np.argwhere(xoutS_all[:, species_ids.index('cPARP')]>100.0)[0]
        if len(death_point)>0: #If cPARP reached a critical concentration; do nothing, the cell likely died
            pass
        
        else: #
            #Take the last values from the first round of simulations and use them as the starting values for the next stimulus phase
            stimulus2_initial_values = xoutS_all[-1]
            
            species_initializations = np.array(model_module.getModel().getInitialStates()) # Get the initial states from the model

            simulation_time = measurement_df['time'].max()/3600 #Note; timepoints in measurements_df should be in seconds, as defined by the SBML file

            # Set the number of records as the number of unique timepoints
            model.setTimepoints(np.linspace(0, simulation_time, len(measurement_df['time'].unique())))


            # Create dynamic unit tests based on PEtab files
            results_dict = {}

            # Here, we pull our unique conditions from the conditions file
            perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, or compartment



    def sparced_erm(yaml_file: str, flagD: Optional[int] = None, flagP: Optional[int] = None):
            """Simulate the experimental replicate model."""

            # Load the PEtab files
            sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

            # Load the SBML model
            model_name = os.path.basename(sbml_file).split('.')[0]
            sys.path.insert(0, os.path.join(os.getcwd(), model_name))
            
            model_module = importlib.import_module(model_name)
            model = model_module.getModel()

            solver = model.getSolver()
            solver.setMaxSteps = 1e10

            # Timepoints are set by the number of unique timepoints and maximum timepoint in the measurement table
            simulation_time = measurement_df['time'].max()/3600

            # Set the number of records as the number of unique timepoints
            model.setTimepoints(np.linspace(0, simulation_time, len(measurement_df['time'].unique())))

            species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt

                        # Set model mathematical representation
            if  flagD == None:
                flagD = 1

            flagD = int(flagD)
            
                        # Load the preincubation step
            if flagP != None:    
                preinc_xoutS_all = SPARCED_ERM.preincubate(yaml_file,flagP)
                species_initializations = preinc_xoutS_all

            else:
                species_initializations = np.array(model_module.getModel().getInitialStates()) # Get the initial states from the model

            simulation_time = measurement_df['time'].max()/3600 #Note; timepoints in measurements_df should be in seconds, as defined by the SBML file

            # Set the number of records as the number of unique timepoints
            model.setTimepoints(np.linspace(0, simulation_time, len(measurement_df['time'].unique())))


            # Create dynamic unit tests based on PEtab files
            results_dict = {}

            # Here, we pull our unique conditions from the conditions file
            perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, or compartment

            unique_conditions = conditions_df.drop_duplicates(subset=perturbants)
            print(unique_conditions)

            for index, condition in unique_conditions.iterrows(): # Iterate through the unique conditions

                species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0 # Set any initializations less than 1e-6 to 0.0

                for entity in perturbants:  # Set the initial concentrations for the perturbants in the conditions table
                    try:
                        # If the entity is a species: change that species value
                        index = species_ids.index(entity)
                        species_initializations[index] = condition[entity]
                    except ValueError:
                        # If the entity is not found in species_ids, move on to the next task
                        pass

                    # If the entity is a parameter: change that parameter's value
                    if entity in parameters_df is not None:
                        model.setParameterById(entity, condition[entity])

                    # If the entity is a compartment: change that compartment's value
                    elif entity in open(sparced_root + '/input_files/Compartments.txt') is not None:
                        compartment = model.getCompartment(entity)
                        compartment.setSize(condition[entity])


                print(f"Running simulation for condition {condition['conditionId']}")
                xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,simulation_time,species_initializations,[],sbml_file,model)
                print("fininshed running simulation")

                iteration_name = condition['conditionId']

                #Build the results dictionary
                results_dict[iteration_name] = {}
                results_dict[iteration_name]['xoutS'] = xoutS_all
                results_dict[iteration_name]['toutS'] = tout_all

            return results_dict


    def stochastic_cell_replicates(yaml_file: str, flagD: Optional[int] = None, flagP: Optional[int] = None, \
                                   num_cells: Optional[int] = None):
        """"create a mock cell population that consists of single cells being ran in unison
        yaml_file: str - path to the YAML file
        flagD: int - 1 for deterministic, 0 for stochastic
        flagP: int - preincubation time in hours
        num_cells: int - number of cells to simulate
        """

        # Load the PEtab files
        sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

        # Load the SBML model
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(os.getcwd(), model_name))
        
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()
        
        mock_cell_population = {}
        for cell in num_cells:
            results_dict = SPARCED_ERM.sparced_erm(yaml_file, flagD, flagP)
            mock_cell_population['cell ' + cell] = results_dict