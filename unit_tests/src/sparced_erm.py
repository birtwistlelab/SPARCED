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
        print(f'heterogenized for {tout_all[-1]/3600} hours')
        # Store the preincubation results in a dictionary
        return xoutS_all[-1]
    
    def set_perturbations(yaml_file, condition, model):
        """Set the perturbations for the simulation."""
        # Load the PEtab files
        sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

        condition = conditions_df[conditions_df['conditionId'] == condition[0]]

        species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt
        
        species_initializations = np.array(model.getInitialStates()) # Get the initial states from the model

        # Pull our unique conditions from the conditions file
        perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, or compartment

        if 'flagD' not in perturbants:
            flagD = 1

        if 'heterogenize' in perturbants: # Dealt with in the sparced_erm function section
            pass
            
        # set simulation perturbants for the pre-equilibration condition
        for entity in perturbants:
            
            #Defines the mathematical representation of the model
            if entity == 'flagD':
                flagD = int(condition[entity])
                if flagD == 0:
                    print(f"Setting gene sampling to stochastic")
                else:
                    print(f"Setting gene sampling to deterministic")
                
            if entity == 'num_cells': 
                pass # This is more importantly used in sparced_erm, as its better this sits outside the for loop

            # If the entity is a compartment: change that compartment's value
            if entity in open('Compartments.txt'):
                compartment = model.getCompartment(entity)
                compartment.setSize(condition[entity]) 

            # If the entity is a parameter: change that parameter's value
            elif entity.strip() in parameters_df.iloc[:, 0].str.strip().tolist():
                try:
                    model.setParameterById(entity, condition[entity].item())
                    parameter_value = model.getParameterById(entity)
                    print(f"Parameter {entity} set to {parameter_value}")
                except RuntimeError:
                    model.setFixedParameterById(entity, str(condition[entity]))
                    print(f"Setting fixed parameter {entity} to {parameter_value}")

            elif entity in species_ids:
                # Set the primary concentrations for the perturbants in the conditions table
                try:
                    # If the entity is a species: change that species value
                    index = species_ids.index(entity)
                    species_initializations[index] = condition[entity]
                    print(f"Setting species {entity} to {species_initializations[index]}")
                except ValueError:
                    print(f"Species {entity} not found")
                    pass

            elif entity.lower().strip() in open('OmicsData.txt', 'r').read().lower().strip(): #This will check the unit test OmicsData,txt file
                print('found in omics data')
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

        return model, species_initializations, flagD

    def sparced_erm(yaml_file: str):
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

        # Create dynamic unit tests based on PEtab files
        results_dict = {}

        # Pull our unique conditions from the conditions file
        perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, or compartment
        
        unique_conditions = conditions_df.drop_duplicates(subset=perturbants)
        print(unique_conditions)


        # Filter out the preequilibration conditions to isolate all unique experimental conditions
        filtered_conditions = [condition for index, condition in unique_conditions.iterrows() \
                            if 'preequilibrationConditionId' not in measurement_df.columns \
                                or condition['conditionId'] not in measurement_df['preequilibrationConditionId'].values]


        for condition in filtered_conditions: # Iterate through the unique conditions 


            iteration_name = condition['conditionId']

            results_dict[iteration_name] = {}


            # Set the number of cells to simulate, 1 by default
            num_cells = condition['num_cells'] if 'num_cells' in condition and condition['num_cells'] is not None else 1


            for cell in range(num_cells):
                
                print(f"Running cell {cell}")  


                if 'heterogenize' in condition and condition['heterogenize'] is not None:
                    heterogenize = condition['heterogenize']
                    preinc_xoutS_all = SPARCED_ERM.heterogenization(yaml_file,heterogenize)
                    model.setInitialStates(preinc_xoutS_all)


                if 'preequilibrationConditionId' in measurement_df.columns:

                    #Isolate the preequilibration condition if included in the measurement table
                    preequilibrate_condition = measurement_df.loc[measurement_df['simulationConditionId'] == condition['conditionId'], \
                                                                'preequilibrationConditionId'].dropna().unique()


                    # Set the preequilibration concentrations for the perturbants in the conditions table
                    preequilibrate_model, species_initializations, flagD = SPARCED_ERM.set_perturbations(
                                                                                                        yaml_file=yaml_file, 
                                                                                                        condition=preequilibrate_condition, 
                                                                                                        model=model
                                                                                                            )
                    
                    
                    # Timepoints are set by the number of unique timepoints and maximum timepoint in the measurement table
                    preequilibrate_simulation_time = measurement_df['time'][measurement_df['preequilibrationConditionId']\
                                                                            .isin(preequilibrate_condition)]\
                                                                                .dropna()\
                                                                                    .unique()\
                                                                                        /3600

                    # Set the number of records as the number of unique timepoints
                    preequilibrate_model.setTimepoints(np.linspace(0, 30, 2))

                    
                    xoutS_all, xoutG_all, tout_all = RunSPARCED(
                                                            flagD=flagD,
                                                                th=preequilibrate_simulation_time,
                                                                    spdata=species_initializations,
                                                                        genedata=[],
                                                                            sbml_file=sbml_file,
                                                                                model=preequilibrate_model
                                                                                )


                    #Build the results dictionary
                    results_dict[iteration_name][f"cell {cell}"] = {}
                    results_dict[iteration_name][f"cell {cell}"]['xoutS'] = xoutS_all
                    results_dict[iteration_name][f"cell {cell}"]['xoutG'] = xoutG_all
                    results_dict[iteration_name][f"cell {cell}"]['toutS'] = tout_all   

                    print(f"setting secondary conditions for {condition['conditionId']}")

                    #Find out if SPARCED died during the first round of stimulus
                    death_point = np.argwhere(xoutS_all[:, species_ids.index('cPARP')]>100.0)
                    if len(death_point)>0: #If cPARP reached a critical concentration; do nothing, the cell likely died
                        pass
                    else:
                        #Use final values of preequilibration as the starting values for the next stimulus phase
                        species_initializations = xoutS_all[-1]

                        # Set the next time frame to simulate
                        secondary_timeframe = (measurement_df['time'][measurement_df['simulationConditionId']\
                                                                        .isin(condition)]\
                                                                            .max()/3600)


                        # Set the secondary concentrations for the perturbants in the conditions table
                        secondary_model, species_initializations2, flagD = SPARCED_ERM.set_perturbations(
                                                                                                    yaml_file=yaml_file,
                                                                                                    condition=condition,
                                                                                                    model=preequilibrate_model
                                                                                                        )                 

                        # Set the number of records as the number of unique timepoints
                        secondary_model.setTimepoints(np.linspace(0, 30, 2))

                        xoutS_all2, xoutG_all2, tout_all2 = RunSPARCED(
                                                                flagD=flagD,
                                                                    th=secondary_timeframe,
                                                                        spdata=species_initializations2,
                                                                            genedata=[],
                                                                                sbml_file=sbml_file,
                                                                                    model=secondary_model
                                                                                    )
                        
                        print("finished running simulation")

                        tout_all2 = tout_all2 + (preequilibrate_simulation_time * 3600 +30)

                        #append secondary conditions to the results dictionary
                        results_dict[iteration_name][f"cell {cell}"]['xoutS'] = np.append(
                                                                                            arr=results_dict[iteration_name][f"cell {cell}"]['xoutS'], 
                                                                                            values=xoutS_all2, 
                                                                                            axis=0
                                                                                            )
                        
                        results_dict[iteration_name][f"cell {cell}"]['xoutG'] = np.append(
                                                                                            arr=results_dict[iteration_name][f"cell {cell}"]['xoutG'], 
                                                                                            values=xoutG_all2, 
                                                                                            axis=0
                                                                                            )
                        
                        results_dict[iteration_name][f"cell {cell}"]['toutS'] = np.append(
                                                                                            arr=results_dict[iteration_name][f"cell {cell}"]['toutS'], 
                                                                                            values=tout_all2, 
                                                                                            axis=0
                                                                                            )
                        
                else: # If there are no preequilibration conditions in the measurement table, run the simulation as normal

                    erm_model, species_initializations, flagD = SPARCED_ERM.set_perturbations(
                                                                                                yaml_file=yaml_file, 
                                                                                                condition=condition, 
                                                                                                model=model
                                                                                                )
                    
                            # Set the next time frame to simulate
                    simulation_timeframe = (
                                            measurement_df['time'][measurement_df['simulationConditionId']\
                                                                .isin(condition)]
                                                                .max()/3600
                                                                    )
                    

                    # Set the number of records as the number of unique timepoints
                    erm_model.setTimepoints(np.linspace(0, 30, 2))
                    
                    # RunSPARCED(flagD,th,spdata,genedata,sbml_file,model)
                    xoutS_all, xoutG_all, tout_all = RunSPARCED(
                                                            flagD=flagD,
                                                                th=simulation_timeframe,
                                                                    spdata=species_initializations,
                                                                        genedata=[],
                                                                            sbml_file=sbml_file,
                                                                                model=erm_model
                                                                                )


                    print("finished running simulation")

                    iteration_name = condition['conditionId']

                    #Build the results dictionary
                    results_dict[iteration_name][f"cell {cell}"] = {}
                    results_dict[iteration_name][f"cell {cell}"]['xoutS'] = xoutS_all
                    results_dict[iteration_name][f"cell {cell}"]['xoutG'] = xoutG_all
                    results_dict[iteration_name][f"cell {cell}"]['toutS'] = tout_all   
                    

        return results_dict
