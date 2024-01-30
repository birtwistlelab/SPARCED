import os
import sys
import importlib
import numpy as np
from petab_file_loader import PEtabFileLoader

class ObservableCalculator:
    def __init__(self, results_dict: str):
        self.results_dict = results_dict
        self.observable_dict = self.observable_calculator()

    def observable_isolator(yaml_file, results_dict):
        """isolate only the observables of interest from the simulation data. Primary function is to cut down on data.
        yaml_file: yaml file containing the PEtab files
        results_dict: dictionary containing the simulation results from SPARCED_ERM
        """
        # Load the PEtab files
        sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

        # Load the SBML model
        current_directory = os.getcwd()
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(current_directory, model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        species_ids = list(model.getStateIds())

        observable_dict = {}

        for condition in results_dict:
            observable_dict[condition] = {}
            for cell in results_dict[condition]:
                observable_dict[condition][cell] = {}
                for _, observable in observable_df.iterrows():
                    obs = [
                        sum(
                            np.array(
                                [
                                    results_dict[condition][cell]['xoutS'][:, species_ids.index(species_name)]
                                    * float(species_compartment)
                                    for species in observable['observableFormula'].split('+')
                                    for species_name, species_compartment in [species.split('*')]
                                ]
                            )
                        )
                    ]

                    observable_name = observable['observableId']
                    observable_dict[condition][cell][observable_name] = {}
                    observable_dict[condition][cell][observable_name]['xoutS'] = sum(obs)
                    observable_dict[condition][cell][observable_name]['toutS'] = results_dict[condition][cell]['toutS']
                    observable_dict[condition][cell][observable_name]['xoutG'] = results_dict[condition][cell]['xoutG']

        return observable_dict

    def species_summation(yaml_file, results_dict): 
        """Calculate observable values from simulation results.
        yaml_file: yaml file containing the PEtab files
        results_dict: dictionary containing the simulation results from SPARCED_ERM
        """
        # Load the PEtab files
        sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

        # Load the SBML model
        current_directory = os.getcwd()
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(current_directory, model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        # Get the condition names
        perturbants = list(conditions_df.columns[2:])

        unique_conditions = conditions_df.drop_duplicates(subset=perturbants)

        unique_timepoints = measurement_df['time'].unique()

        # The first key in the results dictionary is the condition name
        iteration_names = unique_conditions['conditionId'].tolist()

        species_ids = list(model.getStateIds())

        timepoints_of_interest = {}

        observable_dict = {}

        for _, observable in observable_df.iterrows():

            condition_dict = {}

            for condition in iteration_names:
                time_steps = 30 #Not a fan of hard coding, this is a "works for now" solution
                timepoints_of_interest[condition] = {}

                timepoints_of_interest[condition]['toutS'] = []
                timepoints_of_interest[condition]['xoutS'] = []

                for timepoint in unique_timepoints:
                    tp_by_ts = timepoint / time_steps
                    timepoints_of_interest[condition]['toutS'].append(timepoint / time_steps)
                    timepoints_of_interest[condition]['xoutS'].append(list(results_dict[condition]['xoutS'][int(tp_by_ts), :]))

                timepoints_of_interest[condition]['toutS'] = np.array(timepoints_of_interest[condition]['toutS'])
                timepoints_of_interest[condition]['xoutS'] = np.array(timepoints_of_interest[condition]['xoutS'])
                

                obs = [
                    sum(
                        np.array(
                            [
                                timepoints_of_interest[condition]['xoutS'][:, species_ids.index(species_name)]
                                * float(species_compartment)
                                for species in observable['observableFormula'].split('+')
                                for species_name, species_compartment in [species.split('*')]
                            ]
                        )
                    )
                ]
                condition_dict[condition] = {}
                condition_dict[condition]['xoutS'] = sum(obs)
                condition_dict[condition]['toutS'] = timepoints_of_interest[condition]['toutS']
            observable_dict[observable['observableId']] = condition_dict

        return observable_dict
    
    
    def time_to_death(data):   
        """"Returns the time for each simulated cell death for each condition in the results dictionary"""
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


    def death_rate(self, yaml_file, results_dict):
        """Calculate the death rate from simulation results."""
        # Load the PEtab files
        sbml_file, parameters_df, conditions_df, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

        # Load the SBML model
        current_directory = os.getcwd()
        model_name = os.path.basename(self.sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(current_directory, model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        # Get the condition names
        perturbants = list(conditions_df.columns[2:])

        unique_conditions = conditions_df.drop_duplicates(subset=perturbants)

        # The first key in the results dictionary is the condition name
        iteration_names = unique_conditions['conditionId'].tolist()

        species_ids = list(model.getStateIds())

        death_rate_dict = {}

        for condition in iteration_names:
            death_point = np.argwhere(results_dict[condition]['xoutS'][:, species_ids.index('cPARP')]>100.0)[0]
            if len(death_point)>0:
                time_of_death = results_dict[condition]['toutS'][death_point]

        return death_rate_dict
    


    def collect_the_dead(data):
        """Returns the time for each simulated cell death for each condition in the results dictionary"""
        unique_conditions = set([stim for cells in data for stim in data[cells]])
        perished_cells = {stim: [] for stim in unique_conditions}
        for cell in data:
            for stim in data[cell]:
                perished_cell = np.argwhere(data[cell][stim]['xoutS'][:, 105]>100.0)
                if len(perished_cell) > 0: 
                    death_point = perished_cell[0] #Instance where we rule apoptosis is irreversible
                    time_of_death = data[cell][stim]['toutS'][death_point] # Timepoint to match death point. 
                    perished_cells[stim].extend(time_of_death/3600)
                
        # Our returned array now has every cell listed, and the timepoints for when it died for each condition
        # Ex: cell 0: condition1: 48.5hours
        return perished_cells
    
    
    def death_ratios(data, population_size):
        """Returns the ratio of dead cells, should be proceeded by collect_the_dead function"""
        death_ratios = {}
        for condition in data:
            conditional_death_ratio = len(data[condition]) / population_size
            death_ratios[condition] = conditional_death_ratio
        # return death_ratios
        return death_ratios

    def experimental_comparator(yaml_file):
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
        _, _, _, measurement_df, observable_df = PEtabFileLoader.load_petab_files(yaml_file)

        # Group by observableId and simulationConditionId
        grouped_data = measurement_df[measurement_df['preequilibrationConditionId'].isna()].groupby(['observableId', 'simulationConditionId'])

        for (observable, condition), condition_data in grouped_data:
            if condition not in result_dict:
                result_dict[condition] = {}

            if observable not in result_dict[condition]:
                result_dict[condition][observable] = {}

            result_dict[condition][observable]['time'] = condition_data['time'].values
            result_dict[condition][observable]['measurement'] = condition_data['measurement'].values

        return result_dict