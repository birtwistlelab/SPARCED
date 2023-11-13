import os
import shutil
import yaml
import pandas as pd

# Load the PEtab files
class PEtabFileLoader:
    def __init__(self, yaml_file: str):
        self.yaml_file = yaml_file
        self.sbml_file, self.parameter_df, self.conditions_df, self.measurement_df, self.observable_df = self.load_petab_files()

    def load_petab_files( yaml_file: str):
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