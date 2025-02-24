#!/bin/bash python3
"""
Filename: RunSPARCED.py
Created: ?
Author(s): Medhi Bouhaddou, Cemal Erdem, Jonah R. Huggins, Aurore Amrit, and Marc Birtwistle

Description: Simulate an SBML-based model using the SPARCED algorithm.

"""
#<----------------------------package import----------------------------------->
from core.sparced_prep import SPARCEDPrep
from core.threshold_handler import ThresholdHandler
from core.bioexception_handler import ModelExceptions
import numpy as np

#<----------------------------function definition------------------------------>

class RunSPARCED:
    """
    Simulate an AMICI model using the SPARCED algorithm.
    """

    def __init__(self, model_handler, singe_model, duration=0, exchange=30):
        """
        Initialize the RunSPARCED class.
        """
        self.model_handler = model_handler
        self.singe_model = singe_model

        self.prep = SPARCEDPrep(self.model_handler, duration, exchange)

        def _run(self):
            """runs an instance of the SPARCED model."""
            species_list = self.model_handler.sbml_model.get_species_ids()

            mrna_indices = self._get_mrna_indices()

            # TODO: Ignore most of below, copilot made it.
            for step in range(self.prep.step_number):
                # Get the current time
                time = self.prep.time_array[step]

                # Update the model with the current species concentrations
                self.singe_model.update_model(species_conc, self.prep.exchange)

                # Simulate the model for the current time step
                self.singe_model.simulate(self.prep.exchange)

                # Get the updated species concentrations
                species_conc = self.singe_model.get_results_concentrations()

                # Store the updated species concentrations
                self.prep.species_results[step+1, :] = species_conc

                # Get the current species concentrations
                species_conc = self.prep.species_results[step, :]

                # Evaluate exceptions and thresholds
                if self.singe_model.model_exceptions.thresholds is not None:
                    threshold_handler = ThresholdHandler(self.singe_model.model, self.singe_model.model_exceptions)
                    threshold_handler.handle_thresholds()

    def _get_mrna_indices(self, species_list):
        """Gets the index for every mRNA species within the list of global species indices."""
        return [ind for ind, ele in enumerate(species_list) 
                       if ele[:2] == "m_"] # find the indeces for mRNA species
        
