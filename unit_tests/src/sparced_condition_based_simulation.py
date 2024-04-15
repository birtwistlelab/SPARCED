#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------#
# This file conducts user defined, condition-specific simulations using 
# SPARCED and returns the results as nested NumPy arrays.
# ----------------------------------------------------------------------------#

# Import the necessary modules
import os
import sys
import math
import libsbml
import pandas as pd
import numpy as np

# Get the directory path
wd = os.path.dirname(os.path.abspath(__file__))

# Ensure the SPARCED root and bin directories are in the system path
sparced_root = '/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep).index('SPARCED')+1])
sys.path.append(os.path.join(sparced_root, 'bin'))
from unit_test_modules import UnitTestModules as utm
from modules.RunSPARCED import RunSPARCED

class SPARCED_CBS:
    def __init__(self, yaml_file: str, model: str, conditions_df: pd.DataFrame,
                  measurement_df: pd.DataFrame, parameters_df: pd.DataFrame, 
                  sbml_file: str):
        """This class is designed to simulate the experimental replicate model.
        input:
            yaml_file: str - path to the YAML file
            model: str - path to the SBML model
            conditions_df: pd.DataFrame - the conditions dataframe
            measurement_df: pd.DataFrame - the measurement dataframe
            parameters_df: pd.DataFrame - the parameters dataframe
            sbml_file: str - path to the SBML file
            """
        self.model = model.clone()
        self.yaml_file = yaml_file
        self.conditions_df = conditions_df
        self.measurement_df = measurement_df
        self.parameters_df = parameters_df
        self.sbml_file = sbml_file

    def _run_condition_simulation(self, condition: pd.Series):
        """This function runs the simulation for a single condition.
        input:
            condition: pd.Series - the condition to simulate
        output:
            result: pd.DataFrame - the simulation results
        """


        # Look for heterogenize parameters in the condition
        if 'heterogenize' in condition and not math.isnan(condition['heterogenize']):
            
            self.model = self._heterogenize(condition)
            print('heterogenized')
        
        if 'preequilibrationConditionId' in condition and not math.isnan(
            condition['preequilibrationConditionId']):

            self.model = self._preequilibrate(condition)
            print('preequilibrated')

        species_ids = list(self.model.getStateIds())

        # Find out if starting cPARP is too high to start the simulation
        cPARP_value = self.model.getInitialStates()[species_ids.index('cPARP')]
        parp_value = self.model.getInitialStates()[species_ids.index('PARP')]

        if parp_value < cPARP_value:
            print(cPARP_value)
            pass # need to find a new way to move to next loop after this. 

        # Set the perturbations for the simulation

        self.model = self._set_perturbations(condition)

        # Set the timepoints for the simulation
        simulation_timeframe = (
                                self.measurement_df['time']
                                [self.measurement_df['simulationConditionId']\
                                                    .isin(condition)]
                                                    .max()
                                                        )
        
        self.model.setTimepoints(np.linspace(0, 30))

        species_initializations = np.array(self.model.getInitialStates())

        #Find gene sampling method, flagD
        perturbations = list(self.conditions_df.columns[2:]) 
        if 'flagD' in perturbations:
            flagD = condition['flagD']
        else:
            flagD = 1

        # Run the simulation
        if flagD == 0:
            xoutS_all, xoutG_all, tout_all = RunSPARCED(
                                        flagD=flagD,
                                            th=(simulation_timeframe/3600),
                                                spdata=[],
                                                    genedata=[],
                                                        sbml_file=self.sbml_file,
                                                            model=self.model
                                                            )
            
        else:
            #Deterministic simulation (integrated-SBML)
            import amici

            # Set the number of records as the number of unique timepoints
            self.model.setTimepoints(np.linspace(0, simulation_timeframe, 1000))

            solver = self.model.getSolver()
            solver.setMaxSteps = 1e10

            rdata_o4a = amici.runAmiciSimulation(self.model,solver)

            xoutS_all = rdata_o4a['x']
            tout_all = rdata_o4a['t']
            xoutG_all = []

        print(xoutS_all[-1, 105])
        return xoutS_all, tout_all, xoutG_all



    def _preequilibrate(self, condition: pd.Series):      
        """This function assigns a set of conditions that replicate 
        prior experimental conditions before the primary stimulus of 
        interest.
        
        input:
            condition: pd.Series - the condition to simulate
        
        output:
            preequilibrated_model: pd.DataFrame - the preequilibrated model
        """
        # Isolate the preequilibration condition if included in the measurement
        # table
        preequilibrate_condition = (self.measurement_df.loc[
            self.measurement_df['simulationConditionId'] 
            == condition['conditionId'], 'preequilibrationConditionId']
            .dropna()
            .unique())

        # account for no preequilibration condition being found 
        if len(preequilibrate_condition) == 0:
            return self.model
            
        # set perturbations for the simulation
        self.model = self._set_perturbations(condition)

        #Find gene sampling method, flagD
        flagD = (self.conditions_df.loc[
            self.conditions_df['conditionId'] == preequilibrate_condition[0], 
            'flagD']
            .values[0])
        
        # Find the time frame for the preequilibration simulation
        simulation_timeframe = (self.measurement_df['time']
                                [self.measurement_df['preequilibrationConditionId']
                                                    .isin(condition)]
                                                    .max()
                                                        )
        
        species_initializations = np.array(self.model.getInitialStates())


        self.model.setTimepoints(np.linspace(0,30))

        # Run the simulation
        xoutS_all, _, _ = RunSPARCED(flagD, 
                                     simulation_timeframe,
                                     [],
                                     [],
                                     self.sbml_file,
                                     self.model)
        
        # Return the final values
        self.model.setInitialStates(xoutS_all[-1])

        return self.model


    def _set_perturbations(self, condition: pd.Series):
        """This function sets the perturbations for the simulation.
        input:
            condition: pd.Series - the condition to simulate
        output:
            model: libsbml.Model - the updated SBML model
        """
        # Get the perturbations for the condition
        perturbations = list(self.conditions_df.columns[2:]) 

        for perturbant in perturbations:
            try:
                self.model = utm._set_species_value(self.model, perturbant, 
                                                     condition[perturbant])
            except:
                pass

            try:
                self.model = utm._set_parameter_value(self.model, perturbant, 
                                                     condition[perturbant])
            except:
                pass

            try:
                self.model = utm._set_compartmental_volume(self.model, perturbant, 
                                                     condition[perturbant])
            except:
                pass

        return self.model

        
    def _heterogenize(self, condition: pd.Series):
        """This function runs the 'runSPARCED function and returns the final
        values, thus creating the simulated appearance of asynchrony among 
        replicates.
        input:
            condition: pd.Series - the condition to simulate
        output:
            heterogenized_initial_values: pd.DataFrame - the heterogenized
            initial values
        """
        heterogenize = condition['heterogenize']

        simulation_time = int(heterogenize)/3600

        self.model.setTimepoints(np.linspace(0,30))

        growth_factors = ['E', 'H', 'HGF', 'P', 'F', 'I', 'INS']
        
        
        for species in growth_factors:
            self.model = utm._set_species_value(self.model, species, 0) 

        xoutS_all, _, _ = RunSPARCED(flagD=1, 
                                     th=simulation_time,
                                     spdata = [],
                                     genedata = [], 
                                     sbml_file = self.sbml_file, 
                                     model = self.model)


        self.model.setInitialStates(xoutS_all[-1])
    
        return self.model
    
