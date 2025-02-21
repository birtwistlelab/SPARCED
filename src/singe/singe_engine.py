#!/bin/bash python3
#-*- coding: utf-8 -*-
"""
Filename: singe_engine.py
Created: 2025-02-19
Author(s): Jonah R. Huggins

Description: SINGE - (S)tochastic (I)ntegrated (N)etwork for (G)ene (E)xpression
            singe is the stochastic gene expression simulation engine. It is a custom
            implementation of stochastic gene expression, written as the stochastic half
            of the SPARCED algorithm.

            This module contains the SINGE engine, which is responsible for returning an instance
            of the singe engine, which is used to update ODE models with stochastic gene expression.
            at a user-defined exchange rate.
"""
import os
import numpy as np

from singe_prep import SINGEPrep


#<------------------------------Parent Class--------------------------------->
class SINGEEngine:
    """
    SINGE - (S)tochastic (I)ntegrated (N)etwork for (G)ene (E)xpression
    singe is the stochastic gene expression simulation engine. It is a custom
    implementation of stochastic gene expression, written as the stochastic half
    of the SPARCED algorithm.
    """
    def __init__(self, perturbed_model, gene_regulation, omics_data, solver_flag, model_exceptions):
        """
        Initialize the SINGEEngine class.
        """
        self.solver_flag = solver_flag
        self.model_exceptions = model_exceptions

        self.prep = SINGEPrep(perturbed_model, gene_regulation, omics_data, self.solver_flag, self.model_exceptions)

    def update_model(self, model, exchange):
        """
        Update the model with stochastic gene expression.
        """
        # Prepare the model for stochastic gene expression
        self.prep.prepare_model(model, exchange)

    def simulate(self, duration):
        """
        Simulate the model with stochastic gene expression.
        """
        # Simulate the model with stochastic gene expression
        self.prep.simulate(duration)