#!/bin/bash python3
# -*- coding: utf-8 -*-
"""
Filename: sparced_prep.py
Created: 2025-02-20
Author(s): Jonah R. Huggins
Description: SPARCED - (S)BML (P)roliferation (A)poptosis (R)eceptor Signaling
             (C)ell Cycle (E)xpression (D)NA Damage.
            
            This file prepares the core SPARCED algorithm with miscelaneous variables
            needed for simulation. 
"""

#<-----------------------------Import Packages------------------------------->
import os
from types import SimpleNamespace

import numpy as np
#<------------------------------Parent Class--------------------------------->
class SPARCEDPrep:
    """
    Prepares a stochastic gene expression model for simulation within the SPARCED algorithm.
    """
    
    HOURS_TO_SECONDS = 3600

    def __init__(self, model_handler, duration, exchange):
        self.model_handler = model_handler # If TAR names and volumes embedded in SINGE, remove
        self.prep = SimpleNamespace(duration, exchange)
        self._retrieve_steps_number()
        self._retrieve_time_array()
        self._makeEmptyResultsMatrix()
        del self.model_handler


    def _retrieve_steps_number(self):
        """
        Retrieve the number of steps to simulate.
        """
        self.prep.step_number = int(self.prep.duration * self.HOURS_TO_SECONDS / self.prep.exchange)

    def _retrieve_time_array(self):
        """
        Generate the time trajectories.
        """
        self.prep.time_array = np.arange(0, self.prep.duration * self.HOURS_TO_SECONDS + 1, self.prep.exchange)

    def _makeEmptyResultsMatrix(self):
        "Makes empty matrix for the results to be stored in"
        self.prep.species = np.zeros(shape=(self.prep.step_number+1,
                                            len(self.model_handler.sbml_model.getInitialConcentrations())))
        self.prep.species[0,:] = self.model_handler.sbml_model.getInitialConcentrations() # 24hr time point

    def _get_sparced_vals(self):
        """extract hard coded (non-extensible) attributes from SPARCED. 
        Cytoplasmic and nuclear volume handled by SBMLModelHandler. 
        FUTURE: Embed data in SINGE model"""
        self.prep.cytoplasm_volume = self.model_handler.sbml_model.getVolume("Cytoplasm")
        self.prep.nuclear_volume = self.model_handler.sbml_model.getVolume("Nucleus")