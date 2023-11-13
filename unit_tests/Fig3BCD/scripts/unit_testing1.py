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
# Get the directory path
wd = os.path.dirname(os.path.abspath(__file__))

# Create the formatted path
sparced_root = '/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep).index('SPARCED')+1])
sys.path.append(os.path.join(sparced_root, 'bin'))
from modules.RunSPARCED import RunSPARCED

# copy the SBML model into the PEtab input files directory
shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(os.getcwd()), 'petab_files/SPARCED.xml'))


class SPARCEDUnitTester:
    def __init__(self, yaml_file: str):
        self.yaml_file = yaml_file
        self.sbml_file, self.parameter_df, self.conditions_df, self.measurement_df, self.observable_df = self.load_petab_files()

    def load_petab_files(self):
        """Load PETAB files from a YAML file."""
        yaml_directory = os.path.join(os.path.dirname(os.getcwd()), 'petab_files/')
        
        # copy the SBML model into the PEtab input files directory
        if not os.path.exists(os.path.join(yaml_directory, 'SPARCED.xml')):
            shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(os.getcwd()), 'petab_files/SPARCED.xml'))

        with open(self.yaml_file, 'r') as file:
            yaml_dict = yaml.safe_load(file)

        # Construct full paths to petab files based on the YAML file's directory
        sbml_file = os.path.join(yaml_directory, yaml_dict['problems'][0]['sbml_files'][0])

        parameter_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['parameter_file']), sep='\t')
        
        conditions_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['condition_files'][0]), sep='\t')
        
        measurement_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['measurement_files'][0]), sep='\t')
        
        observable_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['observable_files'][0]), sep='\t')

        return sbml_file, parameter_df, conditions_df, measurement_df, observable_df    
    

    def sparced_erm(self, flagD: Optional[int] = None):
            """Simulate the experimental replicate model."""
            
            current_directory = os.getcwd()

            # Load the PEtab files
            # self.sbml_file, self.parameters_df, self.conditions_df, self.measurement_df, self.observable_df = self.load_petab_files(self.yaml_file)

            # Load the SBML model
            model_name = os.path.basename(self.sbml_file).split('.')[0]
            sys.path.insert(0, os.path.join(current_directory, model_name))

            # Set model mathematical representation
            if  flagD == None:
                flagD = 1

            flagD = int(flagD)


            model_module = importlib.import_module(model_name)
            model = model_module.getModel()

            solver = model.getSolver()
            solver.setMaxSteps = 1e10


            # Create dynamic unit tests based on PEtab files
            results_dict = {}

            perturbants = list(self.conditions_df.columns[2:])
            unique_conditions = self.conditions_df.drop_duplicates(subset=perturbants)

            # Timepoints are set by the number of unique timepoints and maximum timepoint in the measurement table
            simulation_time = self.measurement_df['time'].max()
            # Set the number of records as the number of unique timepoints
            # model.setTimepoints(np.linspace(0, simulation_time, len(self.measurement_df['time'].unique())))
            model.setTimepoints(np.linspace(0, simulation_time/60, 1000))


            species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt
            species_initializations = np.array(model_module.getModel().getInitialStates()) # Get the initial states from the model

            for index, condition in unique_conditions.iterrows(): # Iterate through the unique conditions
                species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0 # Set any initializations less than 1e-6 to 0.0
                for species in perturbants:
                    # Set the initial concentrations for the pertu  rbants in the conditions
                    species_initializations[species_ids.index(species)] = condition[species] 


                # Run the simulation
                print(f"Running simulation for condition {condition['conditionId']}")
                xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,simulation_time,species_initializations,[],self.sbml_file,model)

                
                iteration_name = condition['conditionId']

                #Build the results dictionary
                results_dict[iteration_name] = {}
                results_dict[iteration_name]['xoutS'] = xoutS_all
                results_dict[iteration_name]['toutS'] = tout_all

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
