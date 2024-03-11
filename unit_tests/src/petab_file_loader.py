import os
import shutil
import yaml
import pandas as pd

# Load the PEtab files
class PEtabFileLoader:
    """Load PEtab files from a given YAML file.
    input:
        yaml_file: str - path to the YAML file"""
    
    def __init__(self, yaml_file: str) -> str:
        self.yaml_file = yaml_file

    def __call__(self):
        """
        returns: sbml model object -- model object from the SBML file
        parameter dataframe -- currently this is a clone of the input_files/parameters.txt file
        conditions dataframe -- experimental conditions to be tested, based on PEtab formatting
        measurement dataframe -- experimental measurements, linked by observableId and conditionId, as well as simulation times
        observable dataframe -- details of the observables to be tested, including the observableId, observableName, and observableFormula
        model specifications dataframe -- details of the SPARCED model specifications such as gene sampling method,
                                            heterogeneous starting values, and the number of cells to simulate"""

        yaml_directory = os.path.join(os.path.dirname(self.yaml_file))
        
        # copy the SBML model into the PEtab input files directory
        if not os.path.exists(os.path.join(yaml_directory, 'SPARCED.xml')):
            shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(self.yaml_file), 'SPARCED.xml'))

        with open(self.yaml_file, 'r') as file:
            yaml_dict = yaml.safe_load(file)

        # Construct full paths to petab files based on the YAML file's directory
        self.sbml_file = os.path.join(yaml_directory, yaml_dict['problems'][0]['sbml_files'][0])
        self.parameter_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['parameter_file']), sep='\t')
        self.conditions_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['condition_files'][0]), sep='\t')
        self.measurement_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['measurement_files'][0]), sep='\t')
        self.observable_df = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['observable_files'][0]), sep='\t')

        # The model specification files detail 
        if 'model_specification_files' in yaml_dict['problems'][0]:            
            model_specs = pd.read_csv(os.path.join(yaml_directory, yaml_dict['problems'][0]['model_specification_files'][0]), sep='\t')
            self.conditions_df = pd.merge(self.conditions_df, model_specs, on='conditionId')
        
        return self