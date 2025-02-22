#!/bin/bash python3 
"""
Filename: bioexception_handler.py
Created: 2025-02-19
Author(s): Jonah R. Huggins

Description: SINGE - (S)tochastic (I)ntegrated (N)etwork for (G)ene (E)xpression
                singe is the stochastic gene expression simulation engine. It is a custom
                implementation of stochastic gene expression, written as the stochastic half
                of the SPARCED algorithm.
                
                This file specifies a class of functions that take exceptions pertaining
                to model components and attributes and handles them in a way that is
                consistent with the SINGE algorithm.
"""

class ModelExceptions:
    """
    Handles user-defined entities within the stochastic gene expression module that
    require special handling.
    """
    def __init__(self, data, exceptions): 
        self.data = data # A particular dataset that is being handled !!! Maybe too vague...
        self.exceptions = exceptions

    def _set_mRNA_copy_numbers(self):
        """for handling instances of the mRNA Copy Number attribute in the SINGE model"""
        for attr in self.exceptions.omics:
            if attr.key() is 'Exp_RNA':
                for entity, value in self.exceptions.omics.attr:
                    self._setAttr(self.data.omics.mRCN, value)

    def _setAttr(self, attr, value):
        """
        Sets an attribute of a class to a given value.
        """
        setattr(self, attr, value)