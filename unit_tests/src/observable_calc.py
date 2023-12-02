import os
import sys
import importlib
import numpy as np
from petab_file_loader import PEtabFileLoader

class ObservableCalculator:
    def __init__(self, results_dict: str):
        self.results_dict = results_dict
        self.observable_dict = self.observable_calculator()

    def species_summation(self, yaml_file, results_dict): 
        """Calculate observable values from simulation results."""
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


        observable_dict = {}

        for _, observable in observable_df.iterrows():
            condition_dict = {}
            for condition in iteration_names:
                obs = [
                    sum(
                        np.array(
                            [
                                results_dict[condition]['xoutS'][:, species_ids.index(species_name)]
                                * float(species_compartment)
                                for species in observable['observableFormula'].split('+')
                                for species_name, species_compartment in [species.split('*')]
                            ]
                        )
                    )
                ]
                condition_dict[condition] = {}
                condition_dict[condition]['xoutS'] = sum(obs)
                condition_dict[condition]['toutS'] = results_dict[condition]['toutS']
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
    


