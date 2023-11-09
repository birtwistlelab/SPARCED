import sys
import os
import importlib
import shutil
import jdata as jd
import amici
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.legend_handler import Line2D


def load_petab_files(yaml_file):
    """Load petab files from a yaml file.
        yaml: path to yaml file"""
    # Load yaml file
    with open(yaml_file, 'r') as file:
        yaml_dict = yaml.safe_load(file)
    sbml_file = yaml_dict['problems'][0]['sbml_files'][0]
    parameter_df = pd.read_csv(yaml_dict['parameter_file'], sep='\t')
    condition_df = pd.read_csv(yaml_dict['problems'][0]['condition_files'][0], sep='\t')
    measurement_df = pd.read_csv(yaml_dict['problems'][0]['measurement_files'][0], sep='\t')
    observable_df = pd.read_csv(yaml_dict['problems'][0]['observable_files'][0], sep='\t')
    if 'visualization_file' in yaml_dict['problems'][0]:
        visualization_df = pd.read_csv(yaml_dict['problems'][0]['visualization_files'][0], sep='\t') 

    return sbml_file, parameter_df, condition_df, measurement_df, observable_df
    if 'visualization_df' in  locals():
        return visualization_df
    


def sparced_ERM(yaml_file):
    """
    Stands for SPARCED (E)xperimental (R)eplicate (M)odel; a function for replicating experimental data based on the PEtab format.
    model_file: amici.Model instance
    observable: (string) biological observable of interest
    yaml_data: Dictionary containing data from yaml file

    """
    currentDirectory = os.getcwd() #Gets current directory
    sbml_file, parameter_df, conditions_df, measurement_df, observable_file = load_petab_files(yaml_file) #load petab input files
    model_name = os.path.basename(sbml_file).split('.')[0] #Gets the name of the sbml file
    sys.path.insert(0, os.path.join(currentDirectory,model_name)) #inserts model folder created by amici (post model creation) into the path
    model_module = importlib.import_module(model_name) #manually assigns SPARCED into a module to be imported
    model = model_module.getModel() #Retrieves model from using libSBML.getModel() function

    solver = model.getSolver()###UNDERSTAND BETTER###
    solver.setMaxSteps = 1e10###UNDERSTAND BETTER### 
        
    results_dict = {} #create a dictionary to store simulation results
    perturbants = list(conditions_df.columns[2:]) #Columns 3 and onwards are perturbants corresponding to species within species.txt of the SPARCED model
    unique_conditions = conditions_df.drop_duplicates(subset=perturbants) #Drops duplicates, leaving only unique conditions

    simulationTime = measurement_df['time'].max() ###unit of time should be defined by the sbml in PETab's measurements file.
    model.setTimepoints(np.linspace(0,simulationTime,1000)) ###need to set interval based on experimental data; could be number of unique measurement time points? For example: uWEstern has 5
    
    speciesIds = list(model.getStateIds()) #contains a list of all species names within SPARCED 
    speciesInitialStates = np.array(model_module.getModel().getInitialStates()) #Creates an array of initial species states

    for index, condition in unique_conditions.iterrows(): #iterates through each unique condition
            for species in perturbants: #iterates through each perturbant
                    speciesInitialStates[speciesIds.index(species)] = condition[species] #sets the initial state of each species to the value of the perturbant

            model.setInitialStates(speciesInitialStates) #sets the initial states for simulation
            
            simulation = amici.runAmiciSimulation(model,solver) 

            TimePoints = simulation['t']

            iterationName = condition['conditionId']

            results_dict[iterationName] = simulation['x'] #stores simulation results in a dictionary

    return results_dict, TimePoints



def observableCalculator(yaml_file, results_dict):
    """
    Calculates observable values from simulation results.
    yaml_file: path to yaml file
    results_dict: dictionary of simulation results
    """
    currentDirectory = os.getcwd()
    sbml_file, parameter_df, conditions_df, measurement_df, observable_df = load_petab_files(yaml_file)
    model_name = os.path.basename(sbml_file).split('.')[0]
    sys.path.insert(0, os.path.join(currentDirectory,model_name))
    model_module = importlib.import_module(model_name)
    model = model_module.getModel()

    perturbants = list(conditions_df.columns[2:])
    unique_conditions = conditions_df.drop_duplicates(subset=perturbants)
    
    iterationName = unique_conditions['conditionId'].tolist()  # Convert to list for iteration

    speciesIds = list(model.getStateIds())
    
    observable_dict = {}


    for index, observable in observable_df.iterrows():
        condition_dict = {}
        for condition in iterationName:
            obs = [
                sum(
                    np.array(
                        [
                            results_dict[condition][:, speciesIds.index(species_name)]
                            * float(species_compartment)
                            for species in observable['observableFormula'].split('+')
                            for species_name, species_compartment in [species.split('*')]
                        ]
                    )
                )
            ]
            condition_dict[condition] = sum(obs)
        observable_dict[observable['observableId']] = condition_dict

    return observable_dict        


def unitTest(yaml_file):
    """
    Creates a unit test for a given observable.
    yaml_file: path to yaml file
    observable: (string) observableID of interest from observable table, 'observableId' column, defined in yaml file
    """
    #Import PEtab files
    sbml_file, parameter_df, conditions_df, measurement_df, observable_df = load_petab_files(yaml_file)

    #Simulate experiment replicate model (ERM)
    ERM_Data, ERM_Time = sparced_ERM(yaml_file)
 
    #Calculate observable values from ERM
    observables_data = observableCalculator(yaml_file, ERM_Data)
    
    #Create a dictionary of observables
    yamlName = yaml_file.split('.')[0]
    if os.path.exists(os.getcwd() + '/' + yamlName):
        shutil.rmtree(os.getcwd() + '/' + yamlName)
        os.makedirs(os.getcwd() + '/' + yamlName)
        os.chdir(os.getcwd() + '/' + yamlName)
        jd.save(observables_data, yamlName + '.json')
        shutil.copyfile(os.path.dirname(os.getcwd()) +'/'+ yaml_file, os.getcwd()+ '/' + yaml_file)
    else:
        os.makedirs(os.getcwd() + '/' + yamlName)
        os.chdir(os.getcwd() + '/' + yamlName)
        jd.save(observables_data, yamlName + '.json')
        shutil.copyfile(os.path.dirname(os.getcwd()) +'/'+ yaml_file, os.getcwd()+ '/' + yaml_file)