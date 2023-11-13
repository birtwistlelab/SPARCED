import os
import sys
import glob
import shutil
import importlib
import yaml
import pandas as pd
import numpy as np
import amici
import jdata as jd
import argparse
from typing import Optional

# copy the SBML model into the PEtab input files directory
shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(os.getcwd()), 'petab_files/SPARCED.xml'))


class SPARCEDUnitTester:
    def __init__(self, yaml_file: str):
        self.yaml_file = yaml_file
        self.sbml_file, self.parameter_df, self.conditions_df, self.measurement_df, self.observable_df = self.load_petab_files()

    def load_petab_files(self):
        """Load PETAB files from a YAML file."""
        yaml_directory = os.path.join(os.path.dirname(os.getcwd()), 'petab_files/')

        with open(self.yaml_file, 'r') as file:
            yaml_dict = yaml.safe_load(file)

        # Construct full paths to petab files based on the YAML file's directory
        sbml_file = os.path.join(yaml_directory, yaml_dict['problems'][0]['sbml_files'][0])
        parameter_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['parameter_file']), sep='\t')
        conditions_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['condition_files'][0]), sep='\t')
        measurement_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['measurement_files'][0]), sep='\t')
        observable_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['observable_files'][0]), sep='\t')

        return sbml_file, parameter_df, conditions_df, measurement_df, observable_df
    

    def sparced_erm(self):
        """Simulate the experimental replicate model."""
        current_directory = os.getcwd()
        model_name = os.path.basename(self.sbml_file).split('.')[0]

        sys.path.insert(0, os.path.join(current_directory, model_name))

        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        solver = model.getSolver()
        solver.setMaxSteps = 1e10

        results_dict = {}

        perturbants = list(self.conditions_df.columns[2:])
        unique_conditions = self.conditions_df.drop_duplicates(subset=perturbants)

        simulation_time = self.measurement_df['time'].max()
        model.setTimepoints(np.linspace(0, simulation_time, 1000))

        species_ids = list(model.getStateIds())
        species_initial_states = np.array(model_module.getModel().getInitialStates())

        for index, condition in unique_conditions.iterrows():
            for species in perturbants:
                species_initial_states[species_ids.index(species)] = condition[species]

            model.setInitialStates(species_initial_states)

            simulation = amici.runAmiciSimulation(model, solver)

            iteration_name = condition['conditionId']

            results_dict[iteration_name] = {}
            results_dict[iteration_name]['xoutS'] = simulation['x']
            results_dict[iteration_name]['toutS'] = simulation['t']

        return results_dict

    def observable_calculator(self, results_dict):
        """Calculate observable values from simulation results."""
        current_directory = os.getcwd()
        model_name = os.path.basename(self.sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(current_directory, model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        perturbants = list(self.conditions_df.columns[2:])
        unique_conditions = self.conditions_df.drop_duplicates(subset=perturbants)

        iteration_names = unique_conditions['conditionId'].tolist()

        species_ids = list(model.getStateIds())
        observable_dict = {}

        for _, observable in self.observable_df.iterrows():
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
    

    def unit_test(self, observable: Optional[str] = None):
        """Create a unit test for a given observable."""
        experimental_replicate_model = self.sparced_erm()
        observables_data = self.observable_calculator(experimental_replicate_model)

        yaml_name = os.path.basename(self.yaml_file).split('.')[0]

        results_directory = os.path.join(os.path.dirname(os.getcwd()), 'results')

        if not os.path.exists(results_directory):
            os.makedirs(results_directory)

        results_path = os.path.join(results_directory, f"{yaml_name}.json")

        if observable is not None:
            jd.save(observables_data, results_path)
        else:
            jd.save(experimental_replicate_model, results_path)

# Direct path to YAML files
yaml_files_path = os.path.join(os.path.dirname(os.getcwd()), 'petab_files/')

# Find all YAML files in the specified directory
yaml_files = glob.glob(os.path.join(yaml_files_path, '*.yml'))

# Create a unit test for each YAML file

sparc_unit_tester = SPARCEDUnitTester(yaml_files[0])
sparc_unit_tester.unit_test()


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="SPARCED Unit Tester")
#     parser.add_argument("--yaml_file", required=False, help="Path to the YAML file")

#     args = parser.parse_args()

#     sparc_unit_tester = SPARCEDUnitTester(args.yaml_file)
#     sparc_unit_tester.unit_test()
