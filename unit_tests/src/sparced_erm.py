import os
import sys
import math
import importlib
import pandas as pd
import numpy as np
from petab_file_loader import PEtabFileLoader

# Get the directory path
wd = os.path.dirname(os.path.abspath(__file__))

# Ensure the SPARCED root and bin directories are in the system path
sparced_root = '/'.join(wd.split(os.path.sep)[:wd.split(os.path.sep).index('SPARCED')+1])
sys.path.append(os.path.join(sparced_root, 'bin'))
from modules.RunSPARCED import RunSPARCED


class SPARCED_ERM:
    def __init__(self, yaml_file: str):
        self.yaml_file = yaml_file

    def __call__(self):
        """Simulate the experimental replicate model."""

        # Load the PEtab files
        # sbml_file, _, conditions_df, measurement_df, _ = PEtabFileLoader(self.yaml_file).__call__()
        petab_files = PEtabFileLoader(self.yaml_file).__call__()
        sbml_file = petab_files.sbml_file
        conditions_df = petab_files.conditions_df
        measurement_df = petab_files.measurement_df


        # Load the SBML model
        model_name = os.path.basename(sbml_file).split('.')[0]
        sys.path.insert(0, os.path.join(os.getcwd(), model_name))
        model_module = importlib.import_module(model_name)
        model = model_module.getModel()
        solver = model.getSolver()
        solver.setMaxSteps = 1e10


        species_ids = list(model.getStateIds()) # assign species IDs to a list

        results_dict = {} #instantiate the results dictionary

        # Pull our unique conditions from the conditions file
        perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, gene, compartment, or model-specific condition

        unique_conditions = conditions_df.drop_duplicates(subset=perturbants)
        print(unique_conditions)

        # Filter out the preequilibration conditions to isolate all unique experimental conditions
        filtered_conditions = [condition for index, condition in unique_conditions.iterrows() \
                            if 'preequilibrationConditionId' not in measurement_df.columns \
                                or condition['conditionId'] not in measurement_df['preequilibrationConditionId'].values]


        for condition in filtered_conditions:  

            iteration_name = condition['conditionId']

            results_dict[iteration_name] = {} # each condition has its own results section

            # Set the number of cells to simulate, 1 by default
            num_cells = condition['num_cells'] if 'num_cells' in condition and condition['num_cells'] is not None else 1

            for cell in range(num_cells):
                
                model_copy = model.clone() #ensures the model is reset to its original state inbetween conditions, avoiding carryover. 

                print(f"Running cell {cell}")  

                # Set the primary concentrations for species in the model to heterogenized values
                if 'heterogenize' in condition and not math.isnan(condition['heterogenize']): # handles int and empty (NaN) cells in the num_cells column
                    
                    preinc_xoutS_all = SPARCED_ERM._heterogenization(sbml_file, model_copy, condition['heterogenize'])
                    
                    model_copy.setInitialStates(preinc_xoutS_all) # Set the initial states to the heterogenized values

                if 'preequilibrationConditionId' in measurement_df.columns:

                    #Use final values of preequilibration as the starting values for the next stimulus phase
                    model_copy = SPARCED_ERM._preequilibration_condition(self, condition, model_copy, measurement_df, sbml_file)
                    print(f"setting secondary conditions for {condition['conditionId']}")

                #Find out if starting cPARP is too high to start the simulation
                cPARP_value = model_copy.getInitialStates()[species_ids.index('cPARP')]
                parp_value = model_copy.getInitialStates()[species_ids.index('PARP')]

                if parp_value < cPARP_value:
                    break

                erm_model, species_initializations, flagD = SPARCED_ERM._set_perturbations(
                                                                                            yaml_file=self.yaml_file, 
                                                                                            condition=condition, 
                                                                                            model=model_copy
                                                                                            )
                
                        # Set the next time frame to simulate
                simulation_timeframe = (
                                        measurement_df['time'][measurement_df['simulationConditionId']\
                                                            .isin(condition)]
                                                            .max()/3600
                                                                )

                # Set the number of records as the number of unique timepoints
                erm_model.setTimepoints(np.linspace(0, 30, 2))
                
                # RunSPARCED(flagD,th,spdata,genedata,sbml_file,model)
                xoutS_all, xoutG_all, tout_all = RunSPARCED(
                                                        flagD=flagD,
                                                            th=simulation_timeframe,
                                                                spdata=species_initializations,
                                                                    genedata=[],
                                                                        sbml_file=sbml_file,
                                                                            model=erm_model
                                                                            )


                print("finished running simulation")

                iteration_name = condition['conditionId']

                #Build the results dictionary
                results_dict[iteration_name][f"cell {cell}"] = {}
                results_dict[iteration_name][f"cell {cell}"]['xoutS'] = xoutS_all
                results_dict[iteration_name][f"cell {cell}"]['xoutG'] = xoutG_all
                results_dict[iteration_name][f"cell {cell}"]['toutS'] = tout_all   
                

        return results_dict



    def _heterogenization(sbml_file, model, heterogenize: int):
        """Simulate the preincubation step.
        yaml_file: str: The path to the PEtab YAML file
        heterogenize: int: The time to heterogenize for
        
        Returns:
        xoutS_all[-1]: np.array: The last values of the simulation, these are the heterogenized values that inform the conditions start values"""
        th = int(heterogenize)/3600 # convert to hours
        ts = 30

        model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

        # Serum starve the cell by setting the initial conditions of species to 0.0
        species_initializations = np.array(model.getInitialStates())

        # Run SPARCED for preincubation time (th) with stimulus concentrations set to 0
        xoutS_all, _, tout_all = RunSPARCED(0, th,species_initializations,[],sbml_file,model)
        print(f'heterogenized for {tout_all[-1]/3600} hours')
        # return only the last values of the simulation, these are the heterogenized values that inform the conditions start values
        return xoutS_all[-1]



    def _preequilibration_condition(self, condition, model_copy, measurement_df, sbml_file):
        """Isolate the preequilibration condition if included in the measurement table.
        condition: pd.Series: The condition to set the perturbations for
        model_copy: libsbml.Model: The model to set the perturbations for
        measurement_df: pd.DataFrame: The measurement table
        sbml_file: str: The path to the SBML file
        species_ids: list: The species IDs
        
        Returns:
        preequilibrate_model: libsbml.Model: The model with the perturbations set"""

        # Isolate the preequilibration condition if included in the measurement table
        preequilibrate_condition = measurement_df.loc[measurement_df['simulationConditionId'] == condition['conditionId'], \
                                                    'preequilibrationConditionId'].dropna().unique()


        # Set the preequilibration concentrations for the perturbants in the conditions table
        sys.stdout = open(os.devnull, "w") # Suppresses the _set_perturbation print statements, just for this function
        preequilibrate_model, species_initializations, flagD = SPARCED_ERM._set_perturbations(
                                                                                            yaml_file=self.yaml_file, 
                                                                                            condition=preequilibrate_condition, 
                                                                                            model=model_copy
                                                                                                )
        sys.stdout = sys.__stdout__ # Re-enables print statements
        
        # Timepoints are set by the number of unique timepoints and maximum timepoint in the measurement table
        preequilibrate_simulation_time = measurement_df['time'][measurement_df['preequilibrationConditionId']\
                                                                .isin(preequilibrate_condition)]\
                                                                    .dropna()\
                                                                        .unique()\
                                                                            /3600

        # Set the number of records as the number of unique timepoints
        preequilibrate_model.setTimepoints(np.linspace(0, 30, 2))

        # simulation event
        xoutS_all, _, _ = RunSPARCED(
                                    flagD=flagD,
                                        th=preequilibrate_simulation_time,
                                            spdata=species_initializations,
                                                genedata=[],
                                                    sbml_file=sbml_file,
                                                        model=preequilibrate_model
                                                        )


        # Set the initial model states to be the last values of the preequilibrationCondition simulation
        preequilibrate_model.setInitialStates(xoutS_all[-1])

        return preequilibrate_model


    def _set_perturbations(yaml_file, condition, model):
        """Set the perturbations for the simulation.
        yaml_file: str: The path to the PEtab YAML file
        condition: pd.Series: The condition to set the perturbations for
        model: libsbml.Model: The model to set the perturbations for
        
        Returns:
        model: libsbml.Model: The model with the perturbations set
        species_initializations: np.array: The initial species concentrations
        flagD: int: The flagD value"""
        # Load the PEtab files
        # _, parameters_df, conditions_df, _, _= PEtabFileLoader(yaml_file).__call__()
        petab_files = PEtabFileLoader(yaml_file).__call__()
        parameters_df = petab_files.parameter_df
        conditions_df = petab_files.conditions_df


        condition = conditions_df[conditions_df['conditionId'] == condition[0]] # grabs the first condition

        species_ids = list(model.getStateIds()) # Get the species IDs built in from Species.txt
        
        species_initializations = np.array(model.getInitialStates()) # Get the initial states from the model

        # Pull our unique conditions from the conditions file
        perturbants = list(conditions_df.columns[2:]) # can be a species, parameter, or compartment

        if 'flagD' not in perturbants:
            flagD = 1 # Default to deterministic gene sampling

        if 'heterogenize' in perturbants: # Ignore's this to be dealt with separately
            pass
            
        # set simulation perturbants for the pre-equilibration condition
        for entity in perturbants:
            
            #Defines the mathematical representation of the model
            if entity == 'flagD':
                flagD = int(condition[entity])
                if flagD == 0:
                    print(f"Setting gene sampling to stochastic")
                else:
                    print(f"Setting gene sampling to deterministic")
                
            if entity == 'num_cells': 
                pass # This is used in sparced_erm, its better this sits outside the for loop

            # If the entity is a compartment: change that compartment's value
            if entity in open('Compartments.txt'):
                compartment = model.getCompartment(entity)
                compartment.setSize(condition[entity]) 
                print(f"Compartment {entity} set to {compartment.getSize()}")
                
            # If the entity is a parameter: change that parameter's value
            elif entity.strip() in parameters_df.iloc[:, 0].str.strip().tolist():
                try:
                    model.setParameterById(entity, condition[entity].item())
                    parameter_value = model.getParameterById(entity)
                    print(f"Parameter {entity} set to {parameter_value}")
                except RuntimeError:
                    model.setFixedParameterById(entity, str(condition[entity]))
                    print(f"Setting fixed parameter {entity} to {parameter_value}")

            elif entity in species_ids:
                # Set the primary concentrations for the perturbants in the conditions table
                try:
                    # If the entity is a species: change that species value
                    index = species_ids.index(entity)
                    species_initializations[index] = condition[entity]
                    print(f"Setting species {entity} to {species_initializations[index]}")
                except ValueError:
                    print(f"Species {entity} not found")
                    pass

            elif entity.lower().strip() in open('OmicsData.txt', 'r').read().lower().strip(): #This will check the unit test OmicsData,txt file
                print('found in omics data')
                try:
                    with open('OmicsData.txt', 'r') as omics_data_file:
                        omics_data = pd.read_csv(omics_data_file, sep = '\t', index_col=0)
                        # omics_data.loc[entity, 'Exp RNA'] = condition[entity]
                        # print(f"Setting mRNA {entity} copy number to {condition[entity]}")
                        omics_data.loc[entity, 'kTCleak'] = 0.0
                        omics_data.loc[entity, 'kTCmaxs'] = 0.0
                        omics_data.loc[entity, 'kTCd'] = 0.0
                        omics_data.to_csv('OmicsData.txt', sep = '\t')
                        print(f"Turning mRNA {omics_data.loc[entity]} synthesis and degredation off")

                except ValueError:
                    # If the entity is not found in OmicsData: cancel the simulation
                    print(f"Entity {entity} not found!")
                    sys.exit(1)    

        return model, species_initializations, flagD
