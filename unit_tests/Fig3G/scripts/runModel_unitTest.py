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
sys.path.append(os.path.join(sparced_root, 'unit_tests/src'))
from petab_file_loader import PEtabFileLoader
from sparced_erm import SPARCED_ERM
from observable_calc import ObservableCalculator



# copy the SBML model into the PEtab input files directory
shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(os.getcwd()), 'petab_files/SPARCED.xml'))

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
    parser.add_argument('--observable', required=False, type=int, help='the observable of interest from an experiment', default=1)
    parser.add_argument('--name', required=False, type=str, help='the name of the file to save the results', default=None)
    args = parser.parse_args()



def __call__(yaml_file: str, observable: Optional[int] = None, name: Optional[str] = None):
    """Create a unit test for a given observable.
    input:
        yaml_file: str - path to the YAML file
        observable: int - 1 for run with observable, 0 for run without observable
        name: str - name of the file to save the results
    
    output:
        returns the results of the SPARCED model unit test simulation
    """

    experimental_replicate_model = SPARCED_ERM(yaml_file).__call__()

    yaml_name = os.path.basename(yaml_file).split('.')[0]

    results_directory = os.path.join(os.path.dirname(os.getcwd()), 'results')

    if not os.path.exists(results_directory):
        os.makedirs(results_directory)

    results_path = os.path.join(results_directory, f"{yaml_name}.json")


    if observable == 0:
        if name is not None:
            jd.save(experimental_replicate_model, os.path.join(results_directory, f"{name}.json"))

        else:
            jd.save(experimental_replicate_model, results_path)
    else: 
        print("Calculating observable")
        observables_data = ObservableCalculator(yaml_file, experimental_replicate_model).__call__()
        # observables_data = observables_data()
        if name is not None:
            jd.save(observables_data, os.path.join(results_directory, f"{name}.json"))

        else:
            jd.save(observables_data, results_path)



# Direct path to YAML files
yaml_files_path = os.path.join(os.path.dirname(os.getcwd()), 'petab_files/')

# Find all YAML files in the specified directory
yaml_files = glob.glob(os.path.join(yaml_files_path, '*.yml'))

# Create a unit test for each YAML file
__call__(yaml_files[0], observable=args.observable, name=args.name)
