import os
import re
import sys
import importlib
import numpy as np
import pandas as pd
from typing import Optional
from petab_file_loader import PEtabFileLoader

class ObservableCalculator:
    def __init__(self, yaml_file:str, results_dict: dict, observable_df: pd.DataFrame, measurement_df: pd.DataFrame, model: str):
        """This class is designed to calculate observable values from simulation results.
        input:
            yaml_file: str - path to the YAML file
            results_dict: dict - dictionary containing the simulation results
            observable_df: str - path to the observable dataframe
            model: str - path to the SBML model
        """
        self.yaml_file = yaml_file
        self.results_dict = results_dict
        self.observable_df = observable_df
        self.measurement_df = measurement_df
        self.model = model


    def __call__(self):
        """isolate only the observables of interest from the simulation data. Primary function is to cut down on data.

        output: dictionary containing the observables of interest"""


        species_ids = list(self.model.getStateIds()) # assign species IDs to a list

        observable_dict = {} # Instantiate the observable dictionary, dictionary structure will be condition -> cell -> observable -> xoutS, toutS, xoutG
        for condition in self.results_dict:
            observable_dict[condition] = {} # Instatiate the condition dictionary
            for cell in self.results_dict[condition]:
                observable_dict[condition][cell] = {} # Instantiate the cell dictionary
                for _, observable in self.observable_df.iterrows():
                    obs_formula = str(observable['observableFormula'])
                    try:
                        # Search the obs formula for species names
                        species = re.findall(r'\b\w+(?:\.\w+)?\*\w+(?:\.\w+)?\b', obs_formula)

                        for species in species:
                            # Split the species into its name and compartment
                            species_name, species_compartment = species.split('*')
            #                 # Construct the regex pattern to match the species name exactly
                            pattern = r'\b{}\b'.format(re.escape(species_name))
            #                 # Replace only the exact matches of the species name in the formula
                            obs_formula = re.sub(pattern, f'self.results_dict[condition]["{cell}"]["xoutS"][:, species_ids.index("{species_name}")]', obs_formula)

                            #I need to also isolate the parallel gene data saved in xoutG ### JRH keep that in mind
                    except:
                        print(f'Error in observable formula: {obs_formula}')
                        break
                    obs = eval(obs_formula)

                    observable_name = observable['observableId']
                    observable_dict[condition][cell][observable_name] = {}
                    observable_dict[condition][cell][observable_name]['xoutS'] = obs
                    observable_dict[condition][cell][observable_name]['toutS'] = self.results_dict[condition][cell]['toutS']
                    observable_dict[condition][cell][observable_name]['xoutG'] = self.results_dict[condition][cell]['xoutG']

        return observable_dict

    def _sum_unique_dict_entries(self):
        """Sum the unique entries in the results dictionary."""
        unique_entries = {}
        for key, value in self.results_dict.items():
            unique_entries[f'{key}'] = []
            for sub_key, sub_value in value.items():
                unique_entries[f'{key}'].append(1)
            unique_entries[f'{key}'] =  sum(unique_entries[f'{key}'])
        return unique_entries
    
    def _add_experimental_data(self, observable_dict: dict):
        """
        Returns a dictionary of experimental data for each observable and condition,
        matching simulation results dictionary format.

        Parameters:
        - yaml_file (str): Path to the yaml file.

        Returns:
        - result_dict (dict): Dictionary containing experimental data.
        """

        result_dict = observable_dict

        # Group by observableId and simulationConditionId
        if 'preequilibrationConditionId' in self.measurement_df.columns:
            # grouped_data = self.measurement_df[self.measurement_df['preequilibrationConditionId'].isna()].groupby(['observableId', 'simulationConditionId'])
            grouped_data = self.measurement_df.dropna(subset=['preequilibrationConditionId']).groupby(['observableId', 'simulationConditionId'])

        else:
            grouped_data = self.measurement_df.groupby(['observableId', 'simulationConditionId'])
        # look for experimental data in the measurements file by exculding all NaN values in measurement_df['measurement']
        # if all values are NaN, then there is no experimental data to compare to
        if self.measurement_df['measurement'].isna().all():
            print('No experimental data to compare to')
            return observable_dict

        for (observable, condition), condition_data in grouped_data:
            # print(condition_data)
            # print(observable, condition)
            for i in range(0, self._sum_unique_dict_entries()[condition]):
                result_dict[condition][f'cell {i}'][observable]['experiment time'] = condition_data['time'].values
                result_dict[condition][f'cell {i}'][observable]['experiment measurement'] = condition_data['measurement'].values
                # print(f'Added experimental data for {observable} in condition {condition}')
        return result_dict


class CellDeathMetrics:
    def __init__(self, data, observable_name):
        """ This is extended functionality for the observable calculator class. 
        It is designed to calculate different death point metrics for each cell in the simulation results.

        data: dictionary containing the simulation results from ObservableCalculator
        observable_name: name of the observable to be used to determine time of death
        """
        self.data = data
        self.observable_name = observable_name

    def time_to_death(self):
        """"Returns the time for each simulated cell death for each condition in the results dictionary
        
        output: dictionary containing the times to death for each cell per condition"""

        time_of_death = {}
        for condition in self.data:
            time_of_death[condition] = []
            for cell in self.data[condition]:
                dead_simulation = np.array(self.data[condition][cell][self.observable_name]['toutS']\
                                           [self.data[condition][cell][self.observable_name]['xoutS'] > 100.0])
                if dead_simulation.size > 0:
                    dead_simulation_times = dead_simulation[0]/3600             
                    time_of_death[condition].extend(dead_simulation_times.flatten().tolist())
                else:
                    time_of_death[condition].append(None)
        return time_of_death
    
    def average_time_to_death(self):   
        """"Returns the time for the average simulated cell death for each condition in the results dictionary
        
        output: dictionary containing the average time to death for each condition"""

        time_of_death = self.time_to_death()
        for condition, cell_times in self.time_to_death().items():
            time_of_death[condition] = np.mean([time for time in time_of_death[condition] if time is not None])

        return time_of_death

    def death_ratio(self, percent:Optional[bool] = False):
        """Returns the ratio of dead cells, should be proceeded by collect_the_dead function
        time_to_death: dictionary containing the times to death for each cell per condition from the time_to_death function

        output: dictionary containing the ratio of dead cells for each condition"""

        dead_cells = {}
        for condition, cell_times in self.time_to_death().items():
            dead_cells[condition] = 0
            total_cells = len(self.data[condition])
            dead_cells[condition] += sum(1 for time in cell_times if time is not None)
            # dead_cells[condition] = [dead_cells[condition] * (1 / total_cells) if total_cells !=0 else None for condition in dead_cells]
            dead_cells[condition] = dead_cells[condition] / total_cells if total_cells != 0 else None
            if percent:
                dead_cells[condition] = dead_cells[condition] * 100

        return dead_cells
    
    def alive_ratio(self):
        """Returns the ratio of alive cells, should be proceeded by collect_the_dead function
        
        output: dictionary containing the ratio of alive cells for each condition"""

        death_ratio = self.death_ratio()
        alive_ratio = [(1 - x)*100 for x in death_ratio.values()]
        return alive_ratio