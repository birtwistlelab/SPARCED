#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------#
# This file orchestrates the Message Passing Interface (MPI) job querying 
# and data saving of unit tests defined in PEtab files. Input is user 
# defined, condition-specific simulations details, output returns the results
# in a dictionary.
# ----------------------------------------------------------------------------#
# Import required libraries and required internal functions, append to path 
# necessary directories, and instantiate the MPI communicator
# ----------------------------------------------------------------------------#
import os
import sys
import glob
import shutil
import pickle
import argparse
import importlib
from typing import Optional

# Get the directory path
wd = os.path.dirname(os.path.abspath(__file__))

# Create the formatted path to the SPARCED input files
sparced_root = ('/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep)
                                              .index('SPARCED')+1]))

parser = argparse.ArgumentParser(
    description='Provide arguments to build the SPARCED model')
parser.add_argument('--observable', 
                    required=False, 
                    type=int, 
                    help='only the observable in observables.tsv is \
                        calculated (1) or if the entire simulation is \
                            saved (0)', 
                    default=1)
parser.add_argument('--name', 
                    required=False, 
                    type=str, 
                    help='the name of the file to save the results', 
                    default=None)
parser.add_argument('--model', '-m', 
                    required = True, 
                    type=str, 
                    help='path to the model dirctory \
                    to be used for unit testing', 
                    default=(os.path.join(sparced_root, 
                                          'SPARCED/models/SPARCED-standard')))
parser.add_argument('--benchmark', '-b', 
                    required = True, 
                    type=str,
                    help='benchmark to evaluate the model against',
                    default='stochastic_expression')
parser.add_argument('-J', '--inputsim',      default=f"{sparced_root}/SPARCED/data/simulation/",
                    help="relative path to the simulation input directory")
parser.add_argument('-G', '--genereg',       default="GeneReg.txt",
                    help="name of the simulation gene regulation file")
parser.add_argument('-X', '--omics',        default="OmicsData.txt",
                    help="name of the simulation omics data file")


args = parser.parse_args()
f_genereg = args.inputsim + args.genereg
f_omics = args.inputsim + args.omics

# Append utilities and model directories to the path
benchmark_utils_dir = os.path.join(sparced_root, 'benchmarks/benchmark_utils')
sys.path.append(benchmark_utils_dir)
sys.path.append(args.model)
sys.path.append(os.path.join(args.model, 'amici_SPARCED/'))

# Import the required modules
from petab_file_loader import PEtabFileLoader
from unit_test_modules import UnitTestModules as utm
from sparced_condition_based_simulation import SPARCED_CBS as cbs
from observable_calc import ObservableCalculator
from visualization_plotting import VisualizationPlotting
SPARCED = importlib.import_module('SPARCED.SPARCED')

class RunUnitTest:
    """Input the PEtab files and broadcast them to all processes. Then, load 
        the SBML model and create a list of unique conditions. Assign tasks 
        to ranks based on the number of jobs and the number of ranks. Send 
        the results to the root rank and save the results. If applicable, 
        iterate through the observable calculator and save any experimental 
        data with the observable-results.
    input:
        yaml_file: str - path to the YAML file
        observable: int - 1 for run with observable, 
                    0 for run without observable
        name: str - name of the file to save the results
    
    output:
        returns the results of the SPARCED model unit test simulation
    """


    def __init__(self,
                 model_path: str, 
                 benchmark: str,
                 observable: Optional[int] = None, 
                 name: Optional[str] = None):
        yaml_path = os.path.join(sparced_root, 
                                 f'benchmarks/{benchmark}/{benchmark}.yml')
        
        self.yaml_file = yaml_path
        self.model_path = model_path
        self.benchmark = benchmark
        self.observable = observable
        self.name = name


    def __call__(self):
        """Create a unit test for a given observable.
        input:
            yaml_file: str - path to the YAML file
            observable: int - 1 for run with observable, 0 for run without
              observable
            name: str - name of the file to save the results
        
        output:
            returns the results of the SPARCED model unit test simulation
        """
        # ----------------------Broadcast PEtab Files-----------------------#
        # Load the PEtab files  
        petab_files = PEtabFileLoader(self.yaml_file, self.model_path).__call__()
        petab_files_data = {
            'sbml_file': petab_files.sbml_file,
            'conditions_df': petab_files.conditions_df,
            'measurement_df': petab_files.measurement_df,
            'observable_df': petab_files.observable_df,
            'parameter_df': petab_files.parameter_df
        }

        sbml_file = petab_files_data['sbml_file']
        conditions_df = petab_files_data['conditions_df']
        measurement_df = petab_files_data['measurement_df']
        observable_df = petab_files_data['observable_df']
        parameters_df = petab_files_data['parameter_df']

        if 'visualization_df' in petab_files.__dict__:
            petab_files_data['visualization_df'] = petab_files.visualization_df
            visualization_df = petab_files_data['visualization_df']

        results = utm._results_dictionary(conditions_df, measurement_df)
            
        
        # Load the SBML model
        model = SPARCED.getModel()
        solver = model.getSolver()
        solver.setMaxSteps = 1e10

        #---------------------------Simulation--------------------------------#
        # list out the number of conditional simulations
        for index, condition in conditions_df.iterrows():
            for cell in range(condition['num_cells']):

                condition_id = condition['conditionId']
                print(f"Running condition {condition_id} for cell {cell}")
                xoutS, tout, xoutG = (
                                cbs(
                                    yaml_file=self.yaml_file, 
                                    model=model, 
                                    conditions_df=conditions_df, 
                                    measurement_df=measurement_df, 
                                    parameters_df=parameters_df, 
                                    sbml_file=sbml_file,
                                    f_genereg=f_genereg,
                                    f_omics=f_omics)
                                ._run_condition_simulation(condition)
                                )
                print('we proceed to the next condition')
        #------------------------Results Catalogue----------------------------#

                results[condition_id][f'cell {cell}']['xoutS'] = xoutS
                results[condition_id][f'cell {cell}']['toutS'] = tout
                results[condition_id][f'cell {cell}']['xoutG'] = xoutG




        #------------------------Observable Calculation-----------------------#

        # Create a results directory adjacent to scripts directory
        yaml_name = os.path.basename(self.yaml_file).split('.')[0]
        results_directory = os.path.join(
                            self.model_path, 
                            'results'
                            )
        
        if not os.path.exists(results_directory): 
            os.makedirs(results_directory)

        # Final output is saved in pickle format
        results_path = os.path.join(
                        results_directory, 
                        f"{yaml_name}.pkl"
                        )

        # This opts out of the observable calculation, saving all 
        # gene/time/species data
        if self.observable == 0: 
            if self.name is not None:
                results_path = os.path.join(
                                results_directory, 
                                f"{self.name}.pkl"
                                )
                with open(results_path, 'wb') as f:
                    pickle.dump(results, f)                 
            else:
                with open(results_path, 'wb') as f:
                    pickle.dump(results, f)

        else: 
            print("Calculating observable")
            # Instantiate the observable calculator
            observable_calc = ObservableCalculator(
                                    yaml_file=self.yaml_file, 
                                    results_dict=results, 
                                    observable_df=observable_df,
                                    measurement_df=measurement_df,
                                    model=model
                                    )
            # Return the observable data                       
            observables_data = observable_calc.__call__()
            # Add the experimental data to the observable dictionary
            observables_data = (observable_calc
                                ._add_experimental_data(
                observable_dict=observables_data)
                            )
            if self.name is not None:
                results_path = (os.path.join(
                                results_directory, 
                                f"{self.name}.pkl")
                                )
                with open(results_path, 'wb') as f:
                    pickle.dump(observables_data, f)

                if 'visualization_df' in petab_files_data:
                    print('Generating Unit Test Plot')
                    fig = VisualizationPlotting(
                        yaml_file=self.yaml_file, 
                        results_dict=observables_data, 
                        visualization_df=visualization_df, 
                        observable_df=observable_df, 
                        measurement_df=measurement_df
                        ).dynamic_plot()
                
                    # save the figure
                    fig.savefig(
                        os.path.join(
                            results_directory, 
                            f"{self.name}.png"
                            )
                        )
                        
            else:
                # jd.save(observables_data, results_path)
                with open(results_path, 'wb') as f:
                    pickle.dump(observables_data, f)

                if 'visualization_df' in petab_files_data:
                    print('Generating Unit Test Plot')
                    fig = VisualizationPlotting(
                        yaml_file=self.yaml_file, 
                        results_dict=observables_data, 
                        visualization_df=visualization_df, 
                        observable_df=observable_df, 
                        measurement_df=measurement_df
                        ).dynamic_plot()
                
                    # save the figure
                    fig.savefig(
                        os.path.join(
                            results_directory, 
                            f"{yaml_name}.png"
                            )
                        )
# ----------------------------------------------------------------------------#

# Create a unit test for each YAML file
RunUnitTest( model_path=args.model, benchmark = args.benchmark,
            observable=args.observable, name=args.name).__call__()

