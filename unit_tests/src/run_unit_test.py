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
# sys.path.append(os.path.join(sparced_root, 'unit_tests/src'))
src_dir = os.path.join(sparced_root, 'unit_tests/src')
sys.path.append(src_dir)
sys.path.append(os.path.join(src_dir, 'SPARCED'))
from petab_file_loader import PEtabFileLoader
from sparced_erm_multiprocessing import SPARCED_ERM
from observable_calc import ObservableCalculator
import SPARCED

# copy the SBML model into the PEtab input files directory
shutil.copy(os.path.join(src_dir, 'SPARCEDo4a_v1.xml'), 
            os.path.join(os.path.dirname(os.getcwd()), 
                         'petab_files/SPARCED.xml'))


if __name__ == "__main__": 
    # These are our command line arguments
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
    
    args = parser.parse_args()
    communicator = MPI.COMM_WORLD # Create the MPI communicator
    rank = communicator.Get_rank() # The rank of the current process
    size = communicator.Get_size() # Total number of processes assigned

    shutil.copytree(sparced_root + '/input_files', os.path.join(os.getcwd()), 
                    dirs_exist_ok=True)


# ---------------------------------------------------------------------------#
class UnitTest:
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


    def __init__(self, yaml_file: str, observable: Optional[int] = None, 
                 name: Optional[str] = None):
        self.yaml_file = yaml_file
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
        # Broadcast all PEtab files from the root, prevents "Race Conditions".
        if rank == 0:
            petab_files = PEtabFileLoader(self.yaml_file).__call__()  
            petab_files_data = {
                'sbml_file': petab_files.sbml_file,
                'conditions_df': petab_files.conditions_df,
                'measurement_df': petab_files.measurement_df,
                'observable_df': petab_files.observable_df,
                'parameters_df': petab_files.parameter_df
            }
            sbml_file = petab_files_data['sbml_file']
            conditions_df = petab_files_data['conditions_df']
            measurement_df = petab_files_data['measurement_df']
            observable_df = petab_files_data['observable_df']
            parameters_df = petab_files_data['parameters_df']
        else:
            petab_files_data = None

        petab_files_data = communicator.bcast(petab_files_data, root=0)

        # Unpack the broadcasted data
        if rank != 0:
            sbml_file = petab_files_data['sbml_file']
            conditions_df = petab_files_data['conditions_df']
            measurement_df = petab_files_data['measurement_df']
            observable_df = petab_files_data['observable_df']
            parameters_df = petab_files_data['parameters_df']
    

        communicator.Barrier()

        # Load the SBML model
        model = SPARCED.getModel()
        solver = model.getSolver()
        solver.setMaxSteps = 1e10


        # Filter out the preequilibration conditions to isolate all unique
        # experimental conditions
        filtered_conditions = ([
            condition for index, condition in conditions_df.iterrows()
            if 'preequilibrationConditionId' not in measurement_df.columns
            or condition['conditionId'] not in
            measurement_df['preequilibrationConditionId'].values
            ])
        
# ---------------------------------------------------------------------------#
# Conditions are the unique combinations of perturbants in the conditions 
# file. Each condition is a unique simulation.                                 
# ---------------------------------------------------------------------------#

        # This tallies the number of tasks for each unit test
        list_of_jobs = []

        for condition in filtered_conditions: 
            if 'num_cells' not in condition:
                condition_cell = f"{condition['conditionId']}+0"
                list_of_jobs.append(condition_cell)
            else:
                for cell in range(condition['num_cells']):
                    condition_cell = f"{condition['conditionId']}+{cell}"
                    list_of_jobs.append(condition_cell)

        results_dict = {}

        # Assign each rank a set of tasks
        def assign_tasks(rank, total_jobs,size):
            jobs_per_rank = total_jobs // int(size)
            remainder = total_jobs % int(size)
            
            if rank < remainder:
                rank_i_jobs = jobs_per_rank + 1
                start_cell = rank * rank_i_jobs
            else:
                rank_i_jobs = jobs_per_rank
                start_cell = rank * jobs_per_rank + remainder
                
            return start_cell, start_cell + rank_i_jobs


    # The ranks will iterate through each task, then send results to the root
        start_cell, end_cell = assign_tasks(rank, len(list_of_jobs), size)
        print(f' Rank {rank} has {end_cell - start_cell} tasks')
        for i in range(start_cell, end_cell):
            condition_name = list_of_jobs[i].split('+')[0]
            results_dict[condition_name] = {}
            for i in range(start_cell, end_cell):
                cell = list_of_jobs[i].split('+')[1]
                results_dict[condition_name][f'cell {cell}'] = {}

        for i in range(start_cell, end_cell):
            print(f'Rank {rank} is working on task {i}') 
            current_task = list_of_jobs[i]
            condition_name = current_task.split('+')[0]
            cell = current_task.split('+')[1]

            #without the [0] it retruns list[Series] instead of Series
            condition = ([
                condition for condition in filtered_conditions if 
                condition['conditionId'] == condition_name][0]) 

            # Code can be found @ sparced_erm_multiprocessing.py, line 47
            xoutS_all, tout_all, xoutG_all = (
                                    SPARCED_ERM(
                                                yaml_file=self.yaml_file, 
                                                model=model, 
                                                conditions_df=conditions_df, 
                                                measurement_df=measurement_df, 
                                                parameters_df=parameters_df, 
                                                sbml_file=sbml_file)
                                    .__call__(condition)
                                    )
            
            
            results_dict[condition_name][f'cell {cell}']['xoutS'] = xoutS_all
            results_dict[condition_name][f'cell {cell}']['toutS'] = tout_all

            if xoutG_all != []:
                results_dict[condition_name][f'cell {cell}']['xoutG'] = xoutG_all
        
        if rank > 0:
            results_dict_i = results_dict
            communicator.send(results_dict_i, dest=0)
            print(f'Rank {rank} sent results to root')

        # if root; receive the results from all other ranks
        if rank == 0: 
            if size > 1:
                unused_ranks = (
                    size - len(list_of_jobs) if size > len(list_of_jobs) 
                    else 0
                    )
                
                for i in range(1, size - unused_ranks):
                    print('Rank 0 is receiving results from rank', i)
                    results_dict_i = communicator.recv(source=i)
                    print('Rank 0 received results from rank', i)

                    cell = list_of_jobs[i].split('_')[1]

                    for condition in results_dict_i:
                        if condition not in results_dict:
                            results_dict[condition] = {}
                        for cell in results_dict_i[condition]:
                            if cell not in results_dict[condition]:
                                        results_dict[condition][cell] = {}
                    
                    for condition in results_dict_i:
                        for cell in results_dict_i[condition]:
                            results_dict[condition][cell]['xoutS'] = (
                                results_dict_i[condition][cell]['xoutS']
                                )
                            results_dict[condition][cell]['toutS'] = (
                                results_dict_i[condition][cell]['toutS']
                                )
                            if 'xoutG' in results_dict_i[condition][cell]:
                                results_dict[condition][cell]['xoutG'] = (
                                    results_dict_i[condition][cell]['xoutG']
                                    )
                # Create a results directory adjacent to scripts directory
                yaml_name = os.path.basename(self.yaml_file).split('.')[0]
                results_directory = os.path.join(
                                    os.path.dirname(os.getcwd()), 
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
                            pickle.dump(results_dict, f)                 
                    else:
                        with open(results_path, 'wb') as f:
                            pickle.dump(results_dict, f)

                else: 
                    print("Calculating observable")
                    # Instantiate the observable calculator
                    observable_calc = ObservableCalculator(
                                            yaml_file=self.yaml_file, 
                                            results_dict=results_dict, 
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

                    else:
                        # jd.save(observables_data, results_path)
                        with open(results_path, 'wb') as f:
                            pickle.dump(observables_data, f)
# ----------------------------------------------------------------------------#      


# Direct path to YAML files
yaml_files_path = os.path.join(os.path.dirname(os.getcwd()), 'petab_files/')

# Find all YAML files in the specified directory
yaml_files = glob.glob(os.path.join(yaml_files_path, '*.yml'))

# Create a unit test for each YAML file
UnitTest(yaml_files[0], observable=args.observable, name=args.name).__call__()
