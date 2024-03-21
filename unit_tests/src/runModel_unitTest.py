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
        # Load the PEtab files
        petab_files = PEtabFileLoader(self.yaml_file).__call__()
        sbml_file = petab_files.sbml_file
        conditions_df = petab_files.conditions_df
        measurement_df = petab_files.measurement_df
        observable_df = petab_files.observable_df
        parameters_df = petab_files.parameter_df


        # Load the SBML model
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(os.getcwd(), model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()
        solver = model.getSolver()
        solver.setMaxSteps = 1e10
        
        
        # execute the unit test
        # experimental_replicate_model = SPARCED_ERM(self.yaml_file, model, conditions_df, measurement_df, parameters_df, sbml_file).__call__()

        results_dict = {} #instantiate the results dictionary

        # Pull our unique conditions from the conditions file
        perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, gene, compartment, or model-specific condition

        unique_conditions = conditions_df.drop_duplicates(subset=perturbants)
        print(unique_conditions)

        # Filter out the preequilibration conditions to isolate all unique experimental conditions
        filtered_conditions = [condition for index, condition in unique_conditions.iterrows() \
                            if 'preequilibrationConditionId' not in measurement_df.columns \
                                or condition['conditionId'] not in measurement_df['preequilibrationConditionId'].values]

        # --------------------------------------------------------------------------------------------------------------------#
        # Conditions are the unique combinations of perturbants in the conditions file. Each condition is a unique simulation.#
        # --------------------------------------------------------------------------------------------------------------------#
        for condition in filtered_conditions: 

            iteration_name = condition['conditionId']

            results_dict[iteration_name] = {} # each condition has its own results section

            # Set the number of cells to simulate, 1 by default
            num_cells = condition['num_cells'] if 'num_cells' in condition and condition['num_cells'] is not None else 1

            for cell in range(num_cells):
            

                #  cells should be ran in parallel, not conditions
                comm = MPI.COMM_WORLD

                rank = comm.Get_rank()

                size = comm.Get_size()

                xoutS_all, xoutG_all, tout_all = SPARCED_ERM(self.yaml_file, model, conditions_df, measurement_df, parameters_df, sbml_file)\
                    ._process_cell_condition(self, condition, cell, results_dict)
               

                all_results = comm.gather(xoutS_all, root=0)
                print('finished simulation')

                        #Build the results dictionary
                results_dict[iteration_name][f"cell {cell}"] = {}
                results_dict[iteration_name][f"cell {cell}"]['xoutS'] = xoutS_all
                results_dict[iteration_name][f"cell {cell}"]['xoutG'] = xoutG_all
                results_dict[iteration_name][f"cell {cell}"]['toutS'] = tout_all







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
                    pickle.dump(experimental_replicate_model, f)
                # jd.save(experimental_replicate_model, os.path.join(results_directory, f"{name}.json"))

            else:
                # jd.save(experimental_replicate_model, results_path)
                with open(results_path, 'wb') as f:
                    pickle.dump(experimental_replicate_model, f)
        else: 
            print("Calculating observable")
            observables_data = ObservableCalculator(yaml_file=self.yaml_file, 
                                                    results_dict=experimental_replicate_model, 
                                                    observable_df=observable_df, 
                                                    model=model).__call__()

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
