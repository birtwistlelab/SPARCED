import os
import sys
import shutil
from typing import Optional
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
        """Load PETAB files from a YAML file.
        yaml_file: path to yaml file"""
        yaml_directory = os.path.join(os.path.dirname(yaml_file))
        
        # copy the SBML model into the PEtab input files directory
        if not os.path.exists(os.path.join(yaml_directory, 'SPARCED.xml')):
            shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(yaml_file), 'SPARCED.xml'))

        with open(yaml_file, 'r') as file:
            yaml_dict = yaml.safe_load(file)

        # Construct full paths to petab files based on the YAML file's directory
        sbml_file = os.path.join(yaml_directory, yaml_dict['problems'][0]['sbml_files'][0])

        parameter_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['parameter_file']), sep='\t')
        
        conditions_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['condition_files'][0]), sep='\t')
        
        measurement_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['measurement_files'][0]), sep='\t')
        
        observable_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['observable_files'][0]), sep='\t')

        petab_file_names = ['condition_files', 'measurement_files', 'observable_files']

        return sbml_file, parameter_df, conditions_df, measurement_df, observable_df
    

    def load_secondary_files(yaml_file: str):
        """Load PETAB files from a YAML file.
        yaml_file: path to yaml file"""
        yaml_directory = os.path.join(os.path.dirname(yaml_file))
        
        # copy the SBML model into the PEtab input files directory
        if not os.path.exists(os.path.join(yaml_directory, 'SPARCED.xml')):
            shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(yaml_file), 'SPARCED.xml'))

        with open(yaml_file, 'r') as file:
            yaml_dict = yaml.safe_load(file)

        if len(yaml_dict['problems'][0]['condition_files']) > 1:
          conditions_df2 = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['condition_files'][1]), sep='\t')
        else: 
            conditions_df2 = None
        if len(yaml_dict['problems'][0]['measurement_files']) > 1:
            measurement_df2 = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['measurement_files'][1]), sep='\t')
        else:
            measurement_df2 = None

        if len(yaml_dict['problems'][0]['observable_files']) > 1:
            observable_df2 = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['observable_files'][1]), sep='\t')
        else:
            observable_df2 = None

        return conditions_df2, measurement_df2, observable_df2

    def model_loader(yaml_file: str):
        """preload model via normal instructions rather than continuously calling these lines in every file"""
        # Load the PEtab files
        sbml_file, _,_,_,_ = PEtabFileLoader.load_petab_files(yaml_file)

        # Load the SBML model
        current_directory = os.getcwd()
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(current_directory, model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        return model