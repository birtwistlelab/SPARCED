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
from petab_file_loader import PEtabFileLoader
from sparced_erm import SPARCED_ERM
from observable_calc import ObservableCalculator



# copy the SBML model into the PEtab input files directory
shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(os.getcwd()), 'petab_files/SPARCED.xml'))

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
    parser.add_argument('--deterministic', metavar='flagD', type=int, help='1 for deterministic run, 0 for stochastic', default=1)
    args = parser.parse_args()


class SPARCEDUnitTester:
    def __init__(self, yaml_file: str):
        self.yaml_file = yaml_file
    
    def unit_test(self, yaml_file, observable: Optional[str] = None):
        """Create a unit test for a given observable.
        yaml_file: str - path to the YAML file
        observable: str - name of the observable to test
        """
        flagD = args.deterministic
        # Here, we simulate the model 
        experimental_replicate_model = SPARCED_ERM.sparced_erm(yaml_file, flagD)

        if observable is not None:
            observables_data = ObservableCalculator.observable_calculator(yaml_file, experimental_replicate_model)

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
sparc_unit_tester = SPARCEDUnitTester(args.yaml_file)
sparc_unit_tester.unit_test()
