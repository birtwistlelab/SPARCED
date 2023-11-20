import os
import sys
import shutil
import yaml
import pandas as pd
import importlib

# Load the PEtab files
class PEtabFileLoader:
    def __init__(self, yaml_file: str):
        self.yaml_file = yaml_file
        self.sbml_file, self.parameter_df, self.conditions_df, self.measurement_df, self.observable_df = self.load_petab_files()

    # def find_file(file_name, search_path='.'):
    #     for root, dirs, files in os.walk(search_path):
    #         if file_name in files:
    #             return os.path.join(root, file_name)

    #     return None

    def load_petab_files(yaml_file: str):
        """Load PETAB files from a YAML file."""
        yaml_directory = os.path.join(os.path.dirname(os.getcwd()), 'petab_files/')
        
        # copy the SBML model into the PEtab input files directory
        if not os.path.exists(os.path.join(yaml_directory, 'SPARCED.xml')):
            shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(os.getcwd()), 'petab_files/SPARCED.xml'))

        with open(yaml_file, 'r') as file:
            yaml_dict = yaml.safe_load(file)

        # Construct full paths to petab files based on the YAML file's directory
        sbml_file = os.path.join(yaml_directory, yaml_dict['problems'][0]['sbml_files'][0])

        parameter_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['parameter_file']), sep='\t')
        
        conditions_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['condition_files'][0]), sep='\t')
        
        measurement_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['measurement_files'][0]), sep='\t')
        
        observable_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['observable_files'][0]), sep='\t')

        return sbml_file, parameter_df, conditions_df, measurement_df, observable_df
    

    def model_loader(yaml_file: str):
        """Calculate the death rate from simulation results."""
        # Load the PEtab files
        sbml_file, _,_,_,_ = load_petab_files(yaml_file)

        # Load the SBML model
        current_directory = os.getcwd()
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(current_directory, model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        return model