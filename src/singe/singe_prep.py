#!/bin/bash python3
# -*- coding: utf-8 -*-
"""
Filename: singe_prep.py
Created: 2025-02-19
Author(s): Jonah R. Huggins
Description: SINGE - (S)tochastic (I)ntegrated (N)etwork for (G)ene (E)xpression
            singe is the stochastic gene expression simulation engine. It is a custom
            implementation of stochastic gene expression, written as the stochastic half
            of the SPARCED algorithm.
            
            This file prepares a stochastic gene expression model for simulation using 
            the SPARCED algorithm.
"""

#<-----------------------------Import Packages------------------------------->
import os
from src.utils.model_handler import SBMLModelHandler

#<------------------------------Parent Class--------------------------------->
class SINGEPrep:
    """
    Prepares a stochastic gene expression model for simulation within the SPARCED algorithm.
    """
    def __init__(self, ):

