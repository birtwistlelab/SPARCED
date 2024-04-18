#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------#
import libsbml
import numpy as np
import sys

class UnitTestModules:
    """A class for storing helper functions for the unit tests
    """
    @staticmethod
    def _results_size_checker(species_data, time_trajectories):
        """Check the size of the results dictionary
        input:
            species_data: np.array - the species data
            time_trajectories: np.array - the time trajectories
        
        output:
            returns the species and time trajectories, if larger than 2.4gb; 
            trims the results dictionary to every 1/4th of the data
        """
        
        size = sys.getsizeof(species_data)

        threshold_bytes = 2.4 * 1024**3

        while size > threshold_bytes:
            # cut every ith element until the size is less than the threshold
            species_data = species_data[::2]
            time_trajectories = time_trajectories[::2]

            size = sys.getsizeof(species_data)
        print(f"Size of results dictionary: {size} bytes")
        return species_data, time_trajectories


    @staticmethod
    def _tasks_this_round(size, total_jobs, round_number):
        """Calculate the number of tasks for the current round
        input:
            size: int - the total number of processes assigned
            total_jobs: int - the total number of tasks
        
        output:
            returns the number of tasks for the current round
        """
        number_of_rounds = UnitTestModules._number_of_rounds(total_jobs, size)
        # tasks_per_round = total_jobs // size
        tasks_per_round = size
        remainder = total_jobs % size

        if round_number+1 < number_of_rounds:
            tasks_this_round = tasks_per_round
        elif round_number+1 == number_of_rounds and remainder != 0:
            tasks_this_round = remainder
        elif round_number+1 == number_of_rounds and remainder == 0:
            tasks_this_round = tasks_per_round

        return tasks_this_round
    

    @staticmethod
    def _condition_cell_id(rank_task, conditions_df, measurement_df):
        """
        Extract the condition for the task from the filtered_conditions
        output:
            returns the condition for the task
        """
        filtered_conditions = UnitTestModules._filter_conditions(conditions_df, measurement_df)
        
        cell = rank_task.split('+')[1]

        condition_id = rank_task.split('+')[0]
        
        condition = ([
            condition for condition in filtered_conditions if 
            condition['conditionId'] == condition_id][0]) 
    


        return condition, cell, condition_id

    @staticmethod
    def _results_dictionary(conditions_df, measurement_df):
        """Create an empty dictionary for storing results
        input:
            filtered_conditions: pd.DataFrame - filtered conditions dataframe
        
        output:
            returns the empty results dictionary, ready to be filled
        """

        filtered_conditions = UnitTestModules._filter_conditions(conditions_df, measurement_df)

        results = {}

        for condition in filtered_conditions:
            condition_id = condition['conditionId']
            results[condition_id] = {}
            num_cells = condition['num_cells'] if 'num_cells' in condition else 1
            for cell in range(num_cells):
                results[condition_id][f'cell {cell}'] = {}
        return results


    @staticmethod
    def _number_of_rounds(total_jobs, size):
        """Calculate the number of rounds
        input:
            total_jobs: int - the total number of tasks
            size: int - the total number of processes assigned
        
        output:_number_of_rounds
            returns the number of rounds
        """
        rounds_to_complete = total_jobs // size
        remainder = total_jobs % size

        if remainder > 0:
            rounds_to_complete += 1
        return rounds_to_complete


    @staticmethod
    def _total_tasks(conditions_df, measurement_df):
        """Calculate the total number of tasks
        input:
            conditions_df: pd.DataFrame - conditions dataframe
        
        output:
            returns the total number of tasks
        """

        filtered_conditions = UnitTestModules._filter_conditions(conditions_df, measurement_df)

        list_of_jobs = []

        for condition in filtered_conditions: 
            if 'num_cells' not in condition:
                condition_cell = f"{condition['conditionId']}+0"
                list_of_jobs.append(condition_cell)
            else:
                for cell in range(condition['num_cells']):
                    condition_cell = f"{condition['conditionId']}+{cell}"
                    list_of_jobs.append(condition_cell)

        return list_of_jobs


    @staticmethod
    def _filter_conditions(conditions_df, measurement_df):
        """Filter the conditions dataframe to only include unique conditions
        input:
            conditions_df: pd.DataFrame - conditions dataframe
        
        output:
            returns the unique conditions dataframe
        """
        filtered_conditions = ([
            condition for index, condition in conditions_df.iterrows()
            if 'preequilibrationConditionId' not in measurement_df.columns
            or condition['conditionId'] not in
            measurement_df['preequilibrationConditionId'].values
            ])
        
        return filtered_conditions


    @staticmethod
    def _assign_tasks(rank, total_jobs,size):
        """Assign tasks to ranks based on the number of jobs and the number of 
            ranks
        input:
            rank: int - the rank of the current process
            total_jobs: int - the total number of tasks
            size: int - the total number of processes assigned
        output:
            returns the start and end cell for each rank
        """
        jobs_per_rank = total_jobs // int(size)
        remainder = total_jobs % int(size)
        
        if rank < remainder:
            rank_i_jobs = jobs_per_rank + 1
            start_cell = rank * rank_i_jobs
        else:
            rank_i_jobs = jobs_per_rank
            start_cell = rank * jobs_per_rank + remainder
            
        return start_cell, start_cell + rank_i_jobs


    @staticmethod # Not even sure if this one works
    def _set_compartmental_volume(model: libsbml.Model, compartment: str, 
                                  compartment_volume: int):
        """This function sets the volume of a compartment within the SBML model.
        input:
            model: libsbml.Model - the SBML model
            compartment: str - the compartment to set
            compartment_volume: int - the volume to set the compartment to
        output:
            model: libsbml.Model - the updated SBML model
        """
        # Get the list of compartments
        compartment_ids = list(model.getCompartmentIds())

        # assign the compartment volume
        model.setCompartmentVolumeById(compartment, compartment_volume)

        return model
    

    @staticmethod
    def _set_parameter_value(model: libsbml.Model, parameter: str, 
                             parameter_value: int):
        """This function sets the value of a parameter within the SBML model.
        input:
            model: libsbml.Model - the SBML model
            parameter: str - the parameter to set
            parameter_value: int - the value to set the parameter to
        output:
            model: libsbml.Model - the updated SBML model
        """
        # Get the list of parameters
        parameter_ids = list(model.getParameterIds())

        try:# assign the parameter value
            model.setParameterById(parameter, parameter_value)
        except RuntimeError:
            model.setFixedParameterById(parameter, parameter_value)
        
        print(f"Parameter {parameter} set to {model.getParameterById(parameter).getValue()}")

        return model


    # Set this to static method to avoid the need to pass self
    @staticmethod
    def _set_species_value(model: libsbml.Model, species: str, 
                            species_value: int):
        """Thiss function sets the initial value of a species or list of species
        within the sbml model.
        input:
            model: libsbml.Model - the SBML model
            species_value: int - the value to set the species to
            
            output:
                model: libsbml.Model - the updated SBML model"""
        
        # Get the list of species
        species_ids = list(model.getStateIds())
        
        # Get the initial values
        species_initializations = np.array(model.getInitialStates())
        
        # Set the initial values
        index = species_ids.index(species)
        species_initializations[index] = species_value
        model.setInitialStates(species_initializations)
        
        return model