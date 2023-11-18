import os
import sys
import importlib
import numpy as np
from petab_file_loader import PEtabFileLoader

class ObservableCalculator:
    def __init__(self, results_dict: str):
        self.results_dict = results_dict
        self.observable_dict = self.observable_calculator()

    def observable_calculator(self, yaml_file, results_dict):
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

