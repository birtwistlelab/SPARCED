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

    def heterogenization(yaml_file, heterogenize: int):
        """Simulate the preincubation step."""
        th = int(heterogenize)
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


    def sparced_erm(yaml_file: str, flagD: Optional[int] = None, heterogenize: Optional[int] = None, \
                    secondary_conditions: Optional[int] = None):
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

        species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt

                    # Set model mathematical representation
        if  flagD == None:
            flagD = 1

        flagD = int(flagD)
        
                    # Load the preincubation step
        # if heterogenize != None:    
        #     preinc_xoutS_all = SPARCED_ERM.preincubate(yaml_file,heterogenize)
        #     species_initializations = preinc_xoutS_all

        # else:
            # species_initializations = np.array(model_module.getModel().getInitialStates()) # Get the initial states from the model
        species_initializations = np.array(model_module.getModel().getInitialStates()) # Get the initial states from the model
        # Create dynamic unit tests based on PEtab files
        results_dict = {}

        # Here, we pull our unique conditions from the conditions file
        perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, or compartment

        unique_conditions = conditions_df.drop_duplicates(subset=perturbants)
        print(unique_conditions)


        for index, condition in unique_conditions.iterrows(): # Iterate through the unique conditions

            for entity in perturbants:

                if entity == 'heterogenize':
                    heterogenize = condition[entity]
                    preinc_xoutS_all = SPARCED_ERM.heterogenization(yaml_file,heterogenize)
                    species_initializations = preinc_xoutS_all
                    print(f"Cell heterogenized for {heterogenize} hours")

                #Defines the mathematical representation of the model
                if entity == 'flagD':
                    flagD = condition[entity]
                    print(f"Setting solver to {flagD}")
                
                if entity == 'num_cells':
                    num_cells = condition[entity]
                    print(f"Setting number of cells to {num_cells}")

                # If the entity is a compartment: change that compartment's value
                if entity in open('Compartments.txt'):
                    compartment = model.getCompartment(entity)
                    compartment.setSize(condition[entity]) 

                # If the entity is a parameter: change that parameter's value
                elif entity.strip() in parameters_df.iloc[:, 0].str.strip().tolist():
                    try:
                        model.setParameterById(entity, condition[entity])
                        parameter_value = model.getParameterById(entity)
                        print(f"Parameter {entity} set to {parameter_value}")
                    except RuntimeError:
                        model.setFixedParameterById(entity, condition[entity])
                        print(f"Setting fixed parameter {entity} to {parameter_value}")
                    # except RuntimeError:
                    #     print(f"Parameter {entity}  not found")
                    #     pass

                elif entity in species_ids:
                    # Set the primary concentrations for the perturbants in the conditions table
                    try:
                        # If the entity is a species: change that species value
                        index = species_ids.index(entity)
                        species_initializations[index] = condition[entity]
                        print(f"Setting species {entity} to {condition[entity]}")
                    except ValueError:
                        print(f"Species {entity} not found")
                        pass

                elif entity.lower().strip() in open('OmicsData.txt', 'r').read().lower().strip(): #This will check the unit test OmicsData,txt file
                    print('found in omics data')
                    # if "m_" in entity:
                    #     if condition[entity] > 0.0:
                    #         entity = entity.replace("m_", "")
                    try:
                        with open('OmicsData.txt', 'r') as omics_data_file:
                            omics_data = pd.read_csv(omics_data_file, sep = '\t', index_col=0)
                            # omics_data.loc[entity, 'Exp RNA'] = condition[entity]
                            # print(f"Setting mRNA {entity} copy number to {condition[entity]}")
                            omics_data.loc[entity, 'kTCleak'] = 0.0
                            omics_data.loc[entity, 'kTCmaxs'] = 0.0
                            omics_data.loc[entity, 'kTCd'] = 0.0
                            omics_data.to_csv('OmicsData.txt', sep = '\t')
                            print(f"Turning mRNA {omics_data.loc[entity]} synthesis and degredation off")

                    except ValueError:
                        # If the entity is not found in OmicsData: cancel the simulation
                        print(f"Entity {entity} not found!")
                        sys.exit(1)         


            if secondary_conditions is None:
                # Timepoints are set by the number of unique timepoints and maximum timepoint in the measurement table
                simulation_time = (measurement_df['time'].max())/3600

                # Set the number of records as the number of unique timepoints
                model.setTimepoints(np.linspace(0, 30, 2))

                iteration_name = condition['conditionId']
                results_dict[iteration_name] = {}
                
                if 'num_cells' in locals() and num_cells is not None:
                    for cell in range(num_cells):
                        print(f"Running cell {cell}")
                        xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,simulation_time,species_initializations,[],sbml_file,model)
                        print("finished running simulation")

                        iteration_name = condition['conditionId']

                        #Build the results dictionary
                        results_dict[iteration_name][f"cell {cell}"] = {}
                        results_dict[iteration_name][f"cell {cell}"]['xoutS'] = xoutS_all
                        results_dict[iteration_name][f"cell {cell}"]['xoutG'] = xoutG_all
                        results_dict[iteration_name][f"cell {cell}"]['toutS'] = tout_all   

                else:
                    print(f"Running simulation for condition {condition['conditionId']}")   
                    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,simulation_time,species_initializations,[],sbml_file,model)
                    print("finished running simulation")



                    #Build the results dictionary
                    results_dict[iteration_name] = {}
                    results_dict[iteration_name]['xoutS'] = xoutS_all
                    results_dict[iteration_name]['xoutG'] = xoutG_all
                    results_dict[iteration_name]['toutS'] = tout_all 


            else:
                #Here, we simulate the model with the first round of stimulus conditions. 
                first_stimulus_time = secondary_conditions

                # Set the number of records as the number of unique timepoints
                model.setTimepoints(np.linspace(0, 30, 2))

                print(f"Running simulation for condition {condition['conditionId']}")
                xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,first_stimulus_time,species_initializations,[],sbml_file,model)

                iteration_name = condition['conditionId']

                #Build the results dictionary
                results_dict[iteration_name] = {}
                results_dict[iteration_name]['xoutS'] = xoutS_all
                results_dict[iteration_name]['xoutG'] = xoutG_all
                results_dict[iteration_name]['toutS'] = tout_all

                print("Pausing simulation; setting next conditions")
                
                conditions2_df, _, _ = PEtabFileLoader.load_secondary_files(yaml_file)
                secondary_perturbants = list(conditions2_df.columns[2:]) # can be a species, parameter, or compartment

                #Find out if SPARCED died during the first round of stimulus
                death_point = np.argwhere(xoutS_all[:, species_ids.index('cPARP')]>100.0)
                if len(death_point)>0: #If cPARP reached a critical concentration; do nothing, the cell likely died
                    print('Apoptosis occured')
                
                else: 
                    #Take the last values from the first round of simulations and use them as the starting values for the next stimulus phase
                    stimulus2_initial_values = xoutS_all[-1]  
                
                    species_initializations = stimulus2_initial_values


                    unique_conditions_2 = conditions2_df.drop_duplicates(subset=secondary_perturbants)
                    print(unique_conditions_2)

                    for index, second_condition in unique_conditions_2.iterrows(): # Iterate through the unique conditions

                        # species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0 # Set any initializations less than 1e-6 to 0.0

                        for entity in secondary_perturbants:
                            # If the entity is a compartment: change that compartment's value
                            if entity in open(sparced_root + '/input_files/Compartments.txt'):
                                compartment = model.getCompartment(entity)
                                compartment.setSize(second_condition[entity]) 

                            # If the entity is a parameter: change that parameter's value
                            elif any(entity in parameters for parameters in open('ParamsAll.txt', 'r')):
                                parameter_value = second_condition[entity]
                                try:
                                    model.setParameterById(entity, second_condition[entity])
                                except RuntimeError:
                                    pass
                                try:
                                    model.setFixedParameterById(entity, second_condition[entity])
                                except RuntimeError:
                                    pass

                                print(f"Setting parameter {entity} to {parameter_value}")

                            elif entity in species_ids.index(entity):
                                # Set the secondary concentrations for the perturbants in the conditions table
                                try:
                                    # If the entity is a species: change that species value
                                    index = species_ids.index(entity)
                                    species_initializations[index] = second_condition[entity]
                                    print(f"Setting species {entity} to {condition[entity]}")
                                except RuntimeError:
                                    pass
                            
                            elif entity in open('OmicsData.txt'): #This will check the unit test OmicsData,txt file
                                try:
                                    with open('OmicsData.txt', 'r') as omics_data_file:
                                        omics_data = pd.read_csv(omics_data_file, sep = '\t', index_col=0)
                                        omics_data.loc[entity, 'GCN'] = second_condition[entity]

                                except RuntimeError:     
                                    # If the entity is not found in OmicsData: cancel the simulation
                                    print(f"Entity {entity} not found!")
                                    sys.exit(1)

                        # after the second stimulus is added, run the simulation for the remainder of the max timepoint
                        secondary_timeframe = (measurement_df['time'].max()/3600) - first_stimulus_time
                        # Set the number of records as the number of unique timepoints
                        model.setTimepoints(np.linspace(0, 30, 2))
                        print(f"Running simulation for condition {second_condition['conditionId']}")
                        xoutS_all2, xoutG_all2, tout_all2 = RunSPARCED(flagD,secondary_timeframe,species_initializations,[],sbml_file,model)
                        print("Finished running simulation")

                        tout_all2 = tout_all2 + (first_stimulus_time * 3600 +30) # add the first stimulus time to the second stimulus time

                        iteration_name = condition['conditionId']

                        #append secondary conditions to the results dictionary
                        results_dict[iteration_name]['xoutS'] = np.append(results_dict[iteration_name]['xoutS'], xoutS_all2, axis=0)
                        results_dict[iteration_name]['xoutG'] = np.append(results_dict[iteration_name]['xoutG'], xoutG_all2, axis=0)
                        results_dict[iteration_name]['toutS'] = np.append(results_dict[iteration_name]['toutS'], tout_all2, axis=0)

        return results_dict


    def stochastic_cell_replicates(yaml_file: str, flagD: Optional[int] = None, heterogenize: Optional[int] = None, \
                                   secondary_conditions: Optional[int] = None, num_cells: Optional[int] = None):
        """"create a mock cell population that consists of single cells being ran in unison
        yaml_file: str - path to the YAML file
        flagD: int - 1 for deterministic, 0 for stochastic
        flagP: int - preincubation time in hours
        secondary_conditions: int - timeframe to pause the simulation and add any secondary stimuli
        num_cells: int - number of cells to simulate
        """

        stochastic_cell_replicates = {}
        for cell in range(num_cells):
            print(f"Running cell {cell}")
            results_dict = SPARCED_ERM.sparced_erm(yaml_file, flagD, heterogenize, secondary_conditions)
            stochastic_cell_replicates['cell ' + str(cell)] = results_dict

        return stochastic_cell_replicates