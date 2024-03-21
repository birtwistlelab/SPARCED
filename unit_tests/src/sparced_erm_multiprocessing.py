#!/usr/bin/env python
# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This file conducts user defined, condition-specific simulations using SPARCED and returns the results in a dictionary.

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Import required libraries append to path necessary directories
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
import os
import sys
import math
from mpi4py import MPI
import libsbml
import pandas as pd
import numpy as np

# Get the directory path
wd = os.path.dirname(os.path.abspath(__file__))

# Ensure the SPARCED root and bin directories are in the system path
sparced_root = '/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep).index('SPARCED')+1])
sys.path.append(os.path.join(sparced_root, 'bin'))
from modules.RunSPARCED import RunSPARCED

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Define the SPARCED_ERM class, all codes are internal-use only
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class SPARCED_ERM:
    def __init__(self, yaml_file: str, model: str, conditions_df: pd.DataFrame, measurement_df: pd.DataFrame, parameters_df: pd.DataFrame, sbml_file: str):
        """This class is designed to simulate the experimental replicate model.
        input:
            yaml_file: str - path to the YAML file
            model: str - path to the SBML model
            conditions_df: pd.DataFrame - the conditions dataframe
            measurement_df: pd.DataFrame - the measurement dataframe
            parameters_df: pd.DataFrame - the parameters dataframe
            sbml_file: str - path to the SBML file
            """
        self.yaml_file = yaml_file
        self.model = model
        self.conditions_df = conditions_df
        self.measurement_df = measurement_df
        self.parameters_df = parameters_df
        self.sbml_file = sbml_file


    # def __call__(self, condition: pd.Series, cell: int):
    def __call__(self, condition: pd.Series):
        """
        Isolate the preequilibration condition if included in the measurement table.
        conditions_df: pd.DataFrame: The conditions dataframe
        measurement_df: pd.DataFrame: The measurement dataframe
        
        Returns:
        results_dict: dict: The results dictionary
        """
        sys.stdout = open(os.devnull, "w") # REMOVE


        # iteration_name = condition['conditionId']

        model_copy = self.model.clone() #ensures the model is reset to its original state inbetween conditions, avoiding carryover. 

        # print(f"Running cell {cell} in condition {iteration_name}")  

        # Set the primary concentrations for species in the model to heterogenized values
        if 'heterogenize' in condition and not math.isnan(condition['heterogenize']): # handles int and empty (NaN) cells in the num_cells column
            
            preinc_xoutS_all = SPARCED_ERM._heterogenization(self, condition['heterogenize'], model_copy)
            
            model_copy.setInitialStates(preinc_xoutS_all) # Set the initial states to the heterogenized values

        if 'preequilibrationConditionId' in self.measurement_df.columns:

            #Use final values of preequilibration as the starting values for the next stimulus phase
            model_copy = SPARCED_ERM._preequilibration_condition(self, 
                                                                    condition=condition, 
                                                                    model=model_copy)
            
            print(f"setting secondary conditions for {condition['conditionId']}")

        species_ids = list(self.model.getStateIds()) # assign species IDs to a list

        #Find out if starting cPARP is too high to start the simulation
        cPARP_value = model_copy.getInitialStates()[species_ids.index('cPARP')]
        parp_value = model_copy.getInitialStates()[species_ids.index('PARP')]

        if parp_value < cPARP_value:
            pass # need to find a new way to move to next loop after this. 

        erm_model, species_initializations, flagD = SPARCED_ERM._set_perturbations(self,
                                                                                    condition=condition, 
                                                                                    model=model_copy
                                                                                    )
        
                # Set the next time frame to simulate
        simulation_timeframe = (
                                self.measurement_df['time'][self.measurement_df['simulationConditionId']\
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
                                                                sbml_file=self.sbml_file,
                                                                    model=erm_model
                                                                    )


        sys.stdout = sys.__stdout__ # Re-enables print statements
        return xoutS_all, xoutG_all, tout_all
    


    def _heterogenization(self, heterogenize: int, model):
        """Simulate the preincubation step.
        yaml_file: str: The path to the PEtab YAML file
        heterogenize: int: The time to heterogenize for
        
        Returns:
        xoutS_all[-1]: np.array: The last values of the simulation, these are the heterogenized values that inform the conditions start values"""
        th = int(heterogenize)/3600 # convert to hours
        ts = 30

        model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

        # Serum starve the cell by setting the initial conditions of species to 0.0
        species_initializations = np.array(model.getInitialStates())
        growth_factors = ['E', 'H', 'HGF', 'P', 'F', 'I', 'INS']
        species_ids = list(model.getStateIds())
        for species in growth_factors:
            index = species_ids.index(species)
            species_initializations[index] = 0



        # Run SPARCED for preincubation time (th) with stimulus concentrations set to 0
        xoutS_all, _, tout_all = RunSPARCED(0, th,species_initializations,[],self.sbml_file,model)
        print(f'heterogenized for {tout_all[-1]/3600} hours')
        # return only the last values of the simulation, these are the heterogenized values that inform the conditions start values
        return xoutS_all[-1]



    def _preequilibration_condition(self, condition, model):
        """Isolate the preequilibration condition if included in the measurement table.
        condition: pd.Series: The condition to set the perturbations for
        model_copy: libsbml.Model: The model to set the perturbations for
        measurement_df: pd.DataFrame: The measurement table
        sbml_file: str: The path to the SBML file
        species_ids: list: The species IDs
        
        Returns:
        preequilibrate_model: libsbml.Model: The model with the perturbations set"""

        # Isolate the preequilibration condition if included in the measurement table
        preequilibrate_condition = self.measurement_df.loc[self.measurement_df['simulationConditionId'] == condition['conditionId'], \
                                                           'preequilibrationConditionId'].dropna().unique()

        # account for no preequilibration condition being found 
        if len(preequilibrate_condition) == 0:
            return model

        # Set the preequilibration concentrations for the perturbants in the conditions table
        # sys.stdout = open(os.devnull, "w") # Suppresses the _set_perturbation print statements, just for this function
        preequilibrate_model, species_initializations, flagD = SPARCED_ERM._set_perturbations(self,
                                                                                            condition=preequilibrate_condition, 
                                                                                            model=model
                                                                                                )
        # sys.stdout = sys.__stdout__ # Re-enables print statements
        
        # Timepoints are set by the number of unique timepoints and maximum timepoint in the measurement table
        preequilibrate_simulation_time = self.measurement_df['time'][self.measurement_df['preequilibrationConditionId']\
                                                                .isin(preequilibrate_condition)]\
                                                                    .dropna()\
                                                                        .unique()\
                                                                            /3600

        # Set the number of records as the number of unique timepoints
        preequilibrate_model.setTimepoints(np.linspace(0, 30, 2))

        # simulation event
        xoutS_all, _, _ = RunSPARCED(
                                    flagD=flagD,
                                        th=preequilibrate_simulation_time,
                                            spdata=species_initializations,
                                                genedata=[],
                                                    sbml_file=self.sbml_file,
                                                        model=preequilibrate_model
                                                        )


        # Set the initial model states to be the last values of the preequilibrationCondition simulation
        preequilibrate_model.setInitialStates(xoutS_all[-1])

        return preequilibrate_model


    def _set_perturbations(self, condition, model):
        """Set the perturbations for the simulation.
        yaml_file: str: The path to the PEtab YAML file
        condition: pd.Series: The condition to set the perturbations for
        model: libsbml.Model: The model to set the perturbations for
        
        Returns:
        model: libsbml.Model: The model with the perturbations set
        species_initializations: np.array: The initial species concentrations
        flagD: int: The flagD value"""

        condition = self.conditions_df[self.conditions_df['conditionId'] == condition[0]] # grabs the first condition

        species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt
        
        species_initializations = np.array(model.getInitialStates()) # Get the initial states from the model

        # Pull our unique conditions from the conditions file
        perturbants = list(self.conditions_df.columns[2:]) # can be a species, parameter, or compartment

        sbml_reader = libsbml.SBMLReader()
        sbml_doc = sbml_reader.readSBML(self.sbml_file)
        sbml_model = sbml_doc.getModel()

        if 'flagD' not in perturbants:
            flagD = 1 # Default to deterministic gene sampling

        if 'heterogenize' in perturbants: # Ignore's this to be dealt with separately
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
                print('made it to line 246')
                pass # This is used in sparced_erm, its better this sits outside the for loop

            # If the entity is a compartment: change that compartment's value
            compartments = [i.getId() for i in sbml_model.getListOfCompartments()]
            if entity in compartments:
                compartment = sbml_doc.getModel().getCompartment(entity)
                compartment.setVolume(condition[entity].iloc[0]) 
                print(f"Compartment {entity} set to {compartment.getVolume()}")
                
            # If the entity is a parameter: change that parameter's value
                print('made it to line 258')
            elif entity in self.model.getParameterIds():
                try:
                    model.setParameterById(entity, condition[entity].item())
                    parameter_value = model.getParameterById(entity)
                    print(f"Parameter {entity} set to {parameter_value}")
                except RuntimeError:
                    model.setFixedParameterById(entity, str(condition[entity]))
                    print(f"Setting fixed parameter {entity} to {parameter_value}")

            elif entity in species_ids:
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

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------