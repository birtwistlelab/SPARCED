import os
import sys
import glob
import shutil
import importlib
import libsbml
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
sys.path.append(os.path.join(sparced_root, 'bin'))
from modules.RunSPARCED import RunSPARCED

# copy the SBML model into the PEtab input files directory
shutil.copy(os.path.join(os.getcwd(), 'SPARCED.xml'), os.path.join(os.path.dirname(os.getcwd()), 'petab_files/SPARCED.xml'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SPARCED Unit Tester")
    parser.add_argument("--yaml_file", required=False, help="Path to the YAML file")
    parser.add_argument("--deterministic", required=False, help="Deterministic (1) or stochastic simulation (0)",default=1)
    parser.add_argument("--preincubation", required=False, help="Preincubation time in hours",default=24)

    args = parser.parse_args()



def load_petab_files(yaml_file: str):
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



def preincubate(yaml_file, flagP: Optional[int] = None):
    """Simulate the preincubation step."""
    if flagP != None:
        th = flagP
        # Load the PEtab files
        sbml_file, parameters_df, conditions_df, measurement_df, observable_df = load_petab_files(yaml_file)

        # Load the SBML model
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(os.getcwd(), model_name))

        # Here, we import the model's internal packages as a python module
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()

        solver = model.getSolver()
        solver.setMaxSteps = 1e10

        # Serum starve the cell by setting the initial conditions of species to 0.0
        species_initializations = np.array(model_module.getModel().getInitialStates())
        species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0 

        # Run SPARCED for preincubation time (th) with stimulus concentrations set to 0
        xoutS_all, xoutG_all, tout_all = RunSPARCED(0, th,species_initializations,[],sbml_file,model)

        # Store the preincubation results in a dictionary
        return xoutS_all[-1]


def sparced_erm(yaml_file: str, flagD: Optional[int] = None, flagP: Optional[int] = None):
            """Simulate the experimental replicate model."""

            # Load the PEtab files
            sbml_file, parameters_df, conditions_df, measurement_df, observable_df = load_petab_files(yaml_file)

            # Load the SBML model
            model_name = os.path.basename(sbml_file).split('.')[0]
            sys.path.insert(0, os.path.join(os.getcwd(), model_name))

            # Set model mathematical representation
            if  flagD == None:
                flagD = 1

            flagD = int(flagD)

            # Load the preincubation step
            if flagP != None:    
                preinc_xoutS_all = preincubate(yaml_file,flagP)
                species_initializations = preinc_xoutS_all
            ###JRHUGGI: This is a good question for Arnab, how I continue with a preincubation step.
            

            model_module = importlib.import_module(model_name)
            model = model_module.getModel()

            solver = model.getSolver()
            solver.setMaxSteps = 1e10

            # Timepoints are set by the number of unique timepoints and maximum timepoint in the measurement table
            simulation_time = measurement_df['time'].max()/3600

            # Set the number of records as the number of unique timepoints
            model.setTimepoints(np.linspace(0, simulation_time, len(measurement_df['time'].unique())))

            species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt
            species_initializations = np.array(model_module.getModel().getInitialStates()) # Get the initial states from the model


###JRHUGGI: This is where the preincubation step needs to be added to the model
            simulation_time = measurement_df['time'].max()/3600 #Note; timepoints in measurements_df should be in seconds, as defined by the SBML file

            # Set the number of records as the number of unique timepoints
            model.setTimepoints(np.linspace(0, simulation_time, len(measurement_df['time'].unique())))


            # Create dynamic unit tests based on PEtab files
            results_dict = {}

            # Here, we pull our unique conditions from the conditions file
            perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, or compartment

            unique_conditions = conditions_df.drop_duplicates(subset=perturbants)
            print(unique_conditions)

            species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt
            species_initializations = np.array(model_module.getModel().getInitialStates()) # Get the initial states from the model

            for index, condition in unique_conditions.iterrows(): # Iterate through the unique conditions

                species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0 # Set any initializations less than 1e-6 to 0.0

                for entity in perturbants:  # Set the initial concentrations for the perturbants in the conditions table
                    try:
                        # If the entity is a species: change that species value
                        index = species_ids.index(entity)
                        species_initializations[index] = condition[entity]
                    except ValueError:
                        # If the entity is not found in species_ids, move on to the next task
                        pass

                    # If the entity is a parameter: change that parameter's value
                    if entity in parameters_df is not None:
                        model.setParameterById(entity, condition[entity])

                    # If the entity is a compartment: change that compartment's value
                    elif entity in open(sparced_root + '/input_files/Compartments.txt') is not None:
                        compartment = model.getCompartment(entity)
                        compartment.setSize(condition[entity])



                print(f"Running simulation for condition {condition['conditionId']}")
                xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,simulation_time,species_initializations,[],sbml_file,model)
                print("fininshed running simulation")

                iteration_name = condition['conditionId']

                #Build the results dictionary
                results_dict[iteration_name] = {}
                results_dict[iteration_name]['xoutS'] = xoutS_all
                results_dict[iteration_name]['toutS'] = tout_all

            return results_dict


def observable_calculator(yaml_file, results_dict):
    """Calculate observable values from simulation results."""
    current_directory = os.getcwd()

    # Load the PEtab files
    sbml_file, parameters_df, conditions_df, measurement_df, observable_df = load_petab_files(yaml_file)

    model_name = os.path.basename(sbml_file).split('.')[0]
    sys.path.insert(0, os.path.join(current_directory, model_name))
    model_module = importlib.import_module(model_name)
    model = model_module.getModel()

    perturbants = list(conditions_df.columns[2:])
    unique_conditions = conditions_df.drop_duplicates(subset=perturbants)

    iteration_names = unique_conditions['conditionId'].tolist()

    species_ids = list(model.getStateIds())
    observable_dict = {}

    for _, observable in observable_df.iterrows():
        condition_dict = {}
        for condition in iteration_names:
            obs = [
                sum(
                    np.array(
                        [
                            results_dict[condition]['xoutS'][:, species_ids.index(species_name)]
                            * float(species_compartment)
                            for species in observable['observableFormula'].split('+')
                            for species_name, species_compartment in [species.split('*')]
                        ]
                    )
                )
            ]
            condition_dict[condition] = {}
            condition_dict[condition]['xoutS'] = sum(obs)
            condition_dict[condition]['toutS'] = results_dict[condition]['toutS']
        observable_dict[observable['observableId']] = condition_dict

    return observable_dict


def unit_test(yaml_file: str, flagP: Optional[str] = None, observable: Optional[str] = None):
    """Create a unit test for a given observable.
    yaml_file: str - path to the YAML file
    observable: str - name of the observable to test
    """
    if flagP != None:
        #visualing the preincubation step
        preinc = preincubate(yaml_file, flagP=args.preincubation)

    # Here, we simulate the model 
    experimental_replicate_model = sparced_erm(yaml_file)

    if observable is not None:
        observables_data = observable_calculator(experimental_replicate_model)

    yaml_name = os.path.basename(yaml_file).split('.')[0]

    results_directory = os.path.join(os.path.dirname(os.getcwd()), 'results')

    if not os.path.exists(results_directory):
        os.makedirs(results_directory)

    results_path = os.path.join(results_directory, f"{yaml_name}.json")

    if observable is not None:
        jd.save(observables_data, results_path)
    else:
        jd.save(experimental_replicate_model, results_path)

    #save the preinc data
    jd.save(preinc, os.path.join(results_directory, f"{yaml_name}_preinc.json"))
# Direct path to YAML files
yaml_files_path = os.path.join(os.path.dirname(os.getcwd()), 'petab_files/')

# Find all YAML files in the specified directory
yaml_files = glob.glob(os.path.join(yaml_files_path, '*.yml'))

# Create a unit test for each YAML file
unit_test(yaml_files[0], flagP=args.preincubation)