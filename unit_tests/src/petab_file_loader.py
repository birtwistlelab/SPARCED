import os
import sys
import shutil
from typing import Optional
import yaml
import pandas as pd
import importlib

# Load the PEtab files
class PEtabFileLoader:
    def __init__(self, yaml_file: str) -> str:
        self.yaml_file = yaml_file

    def __call__(self):
        """Load PETAB files from a YAML file.
        yaml_file: path to yaml file"""
        yaml_directory = os.path.join(os.path.dirname(self.yaml_file))
        
        # copy the SBML model into the PEtab input files directory
        if not os.path.exists(os.path.join(yaml_directory, 'SPARCED.xml')):
            shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(self.yaml_file), 'SPARCED.xml'))

        with open(self.yaml_file, 'r') as file:
            yaml_dict = yaml.safe_load(file)

        # Construct full paths to petab files based on the YAML file's directory
        sbml_file = os.path.join(yaml_directory, yaml_dict['problems'][0]['sbml_files'][0])

        parameter_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['parameter_file']), sep='\t')
        
        conditions_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['condition_files'][0]), sep='\t')
        
        measurement_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['measurement_files'][0]), sep='\t')
        
        observable_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['observable_files'][0]), sep='\t')

        if 'model_specification_files' in yaml_dict['problems'][0]:

            model_specs = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['model_specification_files'][0]), sep='\t')

            conditions_df = pd.merge(conditions_df, model_specs, on='conditionId')

        return sbml_file, parameter_df, conditions_df, measurement_df, observable_df