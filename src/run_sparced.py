#!/bin/bash python3
"""
Filename: RunSPARCED.py
Created: ?
Author(s): Medhi Bouhaddou, Cemal Erdem, Jonah R. Huggins, Aurore Amrit, and Marc Birtwistle

Description: Simulate an SBML-based model using the SPARCED algorithm.

"""
#<----------------------------package import----------------------------------->
from core.sparced_prep import SPARCEDPrep
import numpy as np

#<----------------------------function definition------------------------------>

class RunSPARCED:
    """
    Simulate an AMICI model using the SPARCED algorithm.
    """

    def __init__(self, model_handler, singe_engine):
        """
        Initialize the RunSPARCED class.
        """
        self.model_handler = model_handler
        self.singe_engine = singe_engine

    def run(self, duration=100, exchange=30):
        """
        Run the SPARCED algorithm.
        """

        # Create the Preparation point
        prep = SPARCEDPrep(self.model_handler)

        ### Find better place for this
        step_number = prep.retrieve_steps_number(duration, exchange)

        # Time trajectories for simulation
        time = prep.time_trajectories(duration, exchange)

        


        # Simulate the model using the SINGE engine
        self.singe_engine.update_model(self.model_handler.model, exchange)
        self.singe_engine.simulate(duration)

        # Return the model
        return self.singe_engine

