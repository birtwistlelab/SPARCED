import os
import sys
import math
import glob
import shutil
import pickle
import argparse
import importlib
import numpy as np
import pandas as pd
from mpi4py import MPI
from typing import Optional

# Get the directory path
wd = os.path.dirname(os.path.abspath(__file__))

# Create the formatted path to the SPARCED input files
sparced_root = '/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep).index('SPARCED')+1])
sys.path.append(os.path.join(sparced_root, 'unit_tests/src'))

from petab_file_loader import PEtabFileLoader
from sparced_erm_multiprocessing import SPARCED_ERM
from observable_calc import ObservableCalculator

# copy the SBML model into the PEtab input files directory
shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(os.getcwd()), 'petab_files/SPARCED.xml'))

if __name__ == "__main__": 
    # These are our command line arguments
    # -- observable dictates whether only the observable in observables.tsv is calculated (default, 1) or if the entire simulation is saved (0)
    # -- name is the name of the file to save the results, default is parent directory name
    parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
    parser.add_argument('--observable', required=False, type=int, help='the observable of interest from an experiment', default=1)
    parser.add_argument('--name', required=False, type=str, help='the name of the file to save the results', default=None)
    args = parser.parse_args()

    communicator = MPI.COMM_WORLD

    rank = communicator.Get_rank()

    size = communicator.Get_size()


class UnitTest:
    """Create a unit test for a given observable.
    input:
        yaml_file: str - path to the YAML file
        observable: int - 1 for run with observable, 0 for run without observable
        name: str - name of the file to save the results
    
    output:
        returns the results of the SPARCED model unit test simulation
    """
    def __init__(self, yaml_file: str, observable: Optional[int] = None, name: Optional[str] = None):
        self.yaml_file = yaml_file
        self.observable = observable
        self.name = name


    def __call__(self):
        """Create a unit test for a given observable.
        input:
            yaml_file: str - path to the YAML file
            observable: int - 1 for run with observable, 0 for run without observable
            name: str - name of the file to save the results
        
        output:
            returns the results of the SPARCED model unit test simulation
        """
        # Broadcast all PEtab files from the root process to all other processes
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
            print('Rank', rank, 'is receiving the broadcasted data')
            sbml_file = petab_files_data['sbml_file']
            conditions_df = petab_files_data['conditions_df']
            measurement_df = petab_files_data['measurement_df']
            observable_df = petab_files_data['observable_df']
            parameters_df = petab_files_data['parameters_df']
    

        communicator.Barrier()

        # Load the SBML model
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(os.getcwd(), model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()
        solver = model.getSolver()
        solver.setMaxSteps = 1e10

        # Pull our unique conditions from the conditions file
        perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, gene, compartment, or model-specific condition

        unique_conditions = conditions_df.drop_duplicates(subset=perturbants)
        # print(unique_conditions)
        
        # Filter out the preequilibration conditions to isolate all unique experimental conditions
        filtered_conditions = [condition for index, condition in unique_conditions.iterrows() \
                        if 'preequilibrationConditionId' not in measurement_df.columns \
                            or condition['conditionId'] not in measurement_df['preequilibrationConditionId'].values]

            # --------------------------------------------------------------------------------------------------------------------#
            # Conditions are the unique combinations of perturbants in the conditions file. Each condition is a unique simulation.#
            # --------------------------------------------------------------------------------------------------------------------#

        list_of_jobs = []

        for condition in filtered_conditions: 
            if 'num_cells' not in condition:
                condition_cell = f"{condition['conditionId']}+0"
                list_of_jobs.append(condition_cell)
            else:
                for cell in range(condition['num_cells']):
                    condition_cell = f"{condition['conditionId']}+{cell}"
                    list_of_jobs.append(condition_cell)
        print(list_of_jobs)
        # if rank > len(list_of_jobs):
        #     print(f'Rank {rank} has no tasks')
        #     return
        
        results_dict = {}

        # Organize ranks by the number of jobs
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


        # Assign tasks to ranks
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

            condition = [condition for condition in filtered_conditions if condition['conditionId'] == condition_name][0] #without the [0] it retruns list[Series] instead of Series

            xoutS_all, xoutG_all, tout_all = SPARCED_ERM(self.yaml_file, model, conditions_df, measurement_df, parameters_df, sbml_file)\
                .__call__(condition)
            
            
            print(f'making a dictionary for rank {rank}')
            results_dict[condition_name][f'cell {cell}']['xoutS'] = xoutS_all
            results_dict[condition_name][f'cell {cell}']['xoutG'] = xoutG_all
            results_dict[condition_name][f'cell {cell}']['toutS'] = tout_all
            # print(results_dict)
        if rank > 0:
            results_dict_i = results_dict
            print(f'Rank {rank} is sending results to root rank')
            communicator.send(results_dict_i, dest=0)
            print(f'Rank {rank} sent results to root rank')

            # print(f'finished simulation for cell {cell} in condition {iteration_name}')
        print(f'Rank {rank} finished its tasks')
        # communicator.Barrier()

        print('the code proceeded through the comm. barrier')
        if rank == 0: # if the root rank, receive the results from all other ranks
            print('on rank 0')
            # I need to make sure that if only one rank is used, th
            if size > 1:
                unused_ranks = size - len(list_of_jobs) if size > len(list_of_jobs) else 0
                print('size is greater than 1')
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
                            results_dict[condition][cell]['xoutS'] = results_dict_i[condition][cell]['xoutS']
                            results_dict[condition][cell]['xoutG'] = results_dict_i[condition][cell]['xoutG']
                            results_dict[condition][cell]['toutS'] = results_dict_i[condition][cell]['toutS']

            # Create a results directory adjacent to scripts directory
            yaml_name = os.path.basename(self.yaml_file).split('.')[0]
            results_directory = os.path.join(os.path.dirname(os.getcwd()), 'results')
            if not os.path.exists(results_directory): # ensures if you've already created the directory, it won't throw an error
                os.makedirs(results_directory)

            results_path = os.path.join(results_directory, f"{yaml_name}.pkl") # final output is saved in pickle format

            if self.observable == 0: # This opts out of the observable calculation, saving al gene/time/species data
                if self.name is not None:
                    results_path = os.path.join(results_directory, f"{self.name}.pkl")
                    with open(results_path, 'wb') as f:
                        pickle.dump(results_dict, f)
                    # jd.save(experimental_replicate_model, os.path.join(results_directory, f"{name}.json"))

                else:
                    # jd.save(experimental_replicate_model, results_path)
                    with open(results_path, 'wb') as f:
                        pickle.dump(results_dict, f)
            else: 
                print("Calculating observable")
                observables_data = ObservableCalculator(yaml_file=self.yaml_file, 
                                                        results_dict=results_dict, 
                                                        observable_df=observable_df,
                                                        measurement_df=measurement_df, 
                                                        model=model).__call__()
                

                observables_data = ObservableCalculator(yaml_file=self.yaml_file, 
                                                        results_dict=results_dict, 
                                                        observable_df=observable_df,
                                                        measurement_df=measurement_df, 
                                                        model=model)._add_experimental_data(observable_dict=observables_data)
                if self.name is not None:
                    results_path = os.path.join(results_directory, f"{self.name}.pkl")
                    # jd.save(observables_data, os.path.join(results_directory, f"{name}.json"))
                    with open(results_path, 'wb') as f:
                        pickle.dump(observables_data, f)

                else:
                    # jd.save(observables_data, results_path)
                    with open(results_path, 'wb') as f:
                        pickle.dump(observables_data, f)
            


# Direct path to YAML files
yaml_files_path = os.path.join(os.path.dirname(os.getcwd()), 'petab_files/')

# Find all YAML files in the specified directory
yaml_files = glob.glob(os.path.join(yaml_files_path, '*.yml'))

# Create a unit test for each YAML file
UnitTest(yaml_files[0], observable=args.observable, name=args.name).__call__()
