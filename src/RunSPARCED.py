#!/bin/bash python3
"""
Filename: RunSPARCED.py
Created: ?
Author(s): Medhi Bouhaddou, Cemal Erdem, Jonah R. Huggins, Aurore Amrit, and Marc Birtwistle

Description: Simulate an SBML-based model using the SPARCED algorithm.

"""
#<----------------------------package import----------------------------------->



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

    def run(self, exchange=30, duration=100):
        """
        Run the SPARCED algorithm.
        """
        # Simulate the model using the SINGE engine
        self.singe_engine.update_model(self.model_handler.model, exchange)
        self.singe_engine.simulate(duration)

        # Return the model
        return self.singe_engine

