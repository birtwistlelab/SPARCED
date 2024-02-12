import os
import re
import sys
import importlib
import numpy as np
from typing import Optional
from petab_file_loader import PEtabFileLoader

class ObservableCalculator:
    def __init__(self, yaml_file:str, results_dict: str):
        """This class is designed to calculate observable values from simulation results."""
        self.yaml_file = yaml_file
        self.results_dict = results_dict


    def __call__(self):
        """isolate only the observables of interest from the simulation data. Primary function is to cut down on data.
        yaml_file: yaml file containing the PEtab files
        results_dict: dictionary containing the simulation results from SPARCED_ERM

        output: dictionary containing the observables of interest"""

        # Load the PEtab files
        sbml_file, _, _, _, observable_df = PEtabFileLoader(self.yaml_file).__call__()

        # Load the SBML model
        current_directory = os.getcwd()
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(current_directory, model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        species_ids = list(model.getStateIds())

        results_dict = self.results_dict

        observable_dict = {}
        for condition in results_dict:
            observable_dict[condition] = {}
            for cell in results_dict[condition]:
                observable_dict[condition][cell] = {}
                for _, observable in observable_df.iterrows():
                    obs_formula = str(observable['observableFormula'])
    
                    # Search the obs formula for species names
                    species = re.findall(r'\b\w+(?:\.\w+)?\*\w+(?:\.\w+)?\b', obs_formula)

                    for species in species:
                        # Split the species into its name and compartment
                        species_name, species_compartment = species.split('*')
        #                 # Construct the regex pattern to match the species name exactly
                        pattern = r'\b{}\b'.format(re.escape(species_name))
        #                 # Replace only the exact matches of the species name in the formula
                        obs_formula = re.sub(pattern, f'results_dict[condition][cell]["xoutS"][:, species_ids.index("{species_name}")]', obs_formula)

                    # print(obs_formula)

                    obs = eval(obs_formula)

                    observable_name = observable['observableId']
                    observable_dict[condition][cell][observable_name] = {}
                    observable_dict[condition][cell][observable_name]['xoutS'] = obs
                    observable_dict[condition][cell][observable_name]['toutS'] = results_dict[condition][cell]['toutS']
                    observable_dict[condition][cell][observable_name]['xoutG'] = results_dict[condition][cell]['xoutG']

        return observable_dict


    def experimental_comparator(self):
        """
        Returns a dictionary of experimental data for each observable and condition,
        matching simulation results dictionary format.

        Parameters:
        - yaml_file (str): Path to the yaml file.

        Returns:
        - result_dict (dict): Dictionary containing experimental data.
        """

        result_dict = {}

        # Load the PEtab files
        _, _, _, measurement_df, observable_df = PEtabFileLoader(self.yaml_file).__call__()

        # Group by observableId and simulationConditionId
        if 'preequilibrationConditionId' in measurement_df.columns:
            grouped_data = measurement_df[measurement_df['preequilibrationConditionId'].isna()].groupby(['observableId', 'simulationConditionId'])

        else:
            grouped_data = measurement_df.groupby(['observableId', 'simulationConditionId'])
            
        for (observable, condition), condition_data in grouped_data:
            if condition not in result_dict:
                result_dict[condition] = {}

            if observable not in result_dict[condition]:
                result_dict[condition][observable] = {}

            result_dict[condition][observable]['time'] = condition_data['time'].values
            result_dict[condition][observable]['measurement'] = condition_data['measurement'].values

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
                dead_simulation = np.argwhere(self.data[condition][cell][self.observable_name]['xoutS'] > 100.0)
                if len(dead_simulation) > 0:
                    time_of_death[condition].append(self.data[condition][cell][self.observable_name]['toutS'][dead_simulation]/3600)
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