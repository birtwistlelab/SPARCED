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
from mpi4py import MPI
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


args = parser.parse_args()

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
# import SPARCED.SPARCED as SPARCED
SPARCED = importlib.import_module('SPARCED.SPARCED')

communicator = MPI.COMM_WORLD # Create the MPI communicator
rank = communicator.Get_rank() # The rank of the current process
size = communicator.Get_size() # Total number of processes assigned

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
        if rank == 0: # Load the PEtab files
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
              
        else:
            petab_files_data = None
            # Broadcast them to all ranks, avoiding Race Conditions
        petab_files_data = communicator.bcast(petab_files_data, root=0)

        if rank != 0:
            sbml_file = petab_files_data['sbml_file']
            conditions_df = petab_files_data['conditions_df']
            measurement_df = petab_files_data['measurement_df']
            observable_df = petab_files_data['observable_df']
            parameters_df = petab_files_data['parameter_df']
            if 'visualization_df' in petab_files_data:
                visualization_df = petab_files_data['visualization_df']

        communicator.Barrier()
        
        # Load the SBML model
        model = SPARCED.getModel()
        solver = model.getSolver()
        solver.setMaxSteps = 1e10

        #----------------------Job Asisgnment-------------------------------#
        list_of_jobs = utm._total_tasks(conditions_df, measurement_df)

        total_jobs = len(list_of_jobs)

        # Catalogue each rank's list of tasks, 
        if rank == 0:
  
            rank_jobs_directory = {}

            for i in range(size):

                start_cell, end_cell = utm._assign_tasks(i, total_jobs, size)
                
                rank_jobs_directory[i] = list_of_jobs[start_cell:end_cell]
                
        # In round iteration, assign each rank it's task for the round
        rounds_to_complete = utm._number_of_rounds(total_jobs, size)
        

        for round_i in range(rounds_to_complete):
        
            if rank == 0:
                print(f'rounds to complete: {rounds_to_complete}')
                print(f'Round {round_i}')
                # Assign each rank it's task for the round
                for i in range(size):

                    rank_jobs = rank_jobs_directory[i]
                    # it needs to handle if there are no more tasks
                    if round_i < len(rank_jobs):
                        rank_job_for_round = rank_jobs[round_i]
                    else:
                        rank_job_for_round = None
                    communicator.send(rank_job_for_round, dest=i, tag=round_i)

            # Receive the task for the round
            rank_task = communicator.recv(source=0, tag=round_i)

            communicator.Barrier()
        #---------------------------Simulation--------------------------------#
            if rank_task is None:
                continue

            condition, cell, condition_id = utm._condition_cell_id(rank_task, 
                                                            conditions_df, 
                                                            measurement_df)
            print(f"Rank {rank} is running {condition_id} for cell {cell}")
            
            print(model.getInitialStates())

            xoutS, tout, xoutG = (
                            cbs(
                                yaml_file=self.yaml_file, 
                                model=model, 
                                conditions_df=conditions_df, 
                                measurement_df=measurement_df, 
                                parameters_df=parameters_df, 
                                sbml_file=sbml_file)
                            ._run_condition_simulation(condition)
                            )
        
        #------------------------Results Catalogue----------------------------#
            rank_results = {
                'condition_name': condition_id,
                'cell': cell,
                'xoutS': xoutS,
                'toutS': tout,
            }
            if xoutG != []:
                rank_results['xoutG'] = xoutG


           
            rank_results['xoutS'], rank_results['toutS'] = utm._results_size_checker(
                                                                rank_results['xoutS'], 
                                                                rank_results['toutS']
                                                                                    )
            

            if rank == 0:
                results[condition_id][f'cell {cell}']['xoutS'] = rank_results['xoutS']
                results[condition_id][f'cell {cell}']['toutS'] = rank_results['toutS']
                if xoutG != []:
                    results[condition_id][f'cell {cell}']['xoutG'] = rank_results['xoutG']
                    print('rank 0 catalogued')

                tasks_this_round = utm._tasks_this_round(size, total_jobs, round_i) - 1
                print(f'tasks this round: {tasks_this_round + 1}')
                completed_tasks = 0
                while completed_tasks < tasks_this_round:
                    print('receiving')
                    rank_results = communicator.recv(source=MPI.ANY_SOURCE, tag = MPI.ANY_TAG)
                    print(f'received results')
                    condition_id = rank_results['condition_name']
                    cell = rank_results['cell']
                    results[condition_id][f'cell {cell}']['xoutS'] = rank_results['xoutS']
                    results[condition_id][f'cell {cell}']['toutS'] = rank_results['toutS']
                    if xoutG != []:
                        results[condition_id][f'cell {cell}']['xoutG'] = rank_results['xoutG']
                    completed_tasks += 1
                    print(f'completed tasks: {completed_tasks}')
                    if completed_tasks == tasks_this_round:
                        print('breaking')
                        break
            
            else:
                print(f'Rank {rank} sending')
                communicator.send(rank_results, dest=0, tag=round_i)
                print(f'Rank {rank} sent')


        #------------------------Observable Calculation-----------------------#
        if rank == 0:
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

