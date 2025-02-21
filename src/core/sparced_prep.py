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
from src.utils.model_handler import SBMLModelHandler

import numpy as np
#<------------------------------Parent Class--------------------------------->
class SPARCEDPrep:
    """
    Prepares a stochastic gene expression model for simulation within the SPARCED algorithm.
    """
    
    HOURS_TO_SECONDS = 3600

    def __init__(self, model_handler):
        self.model_handler = model_handler

    def retrieve_steps_number(self, duration: float, exchange: float):
        """
        Retrieve the number of steps to simulate.
        """
        return int(duration * self.HOURS_TO_SECONDS / exchange)

    def time_trajectories(self, duration: float, exchange: float):
        """
        Generate the time trajectories.
        """
        return np.arange(0, duration * self.HOURS_TO_SECONDS + 1, exchange)
