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
#<-----------------------------Import Packages------------------------------->

from threshold_handler import ThresholdHandler

class ModelExceptions:
    """
    Handles user-defined entities within the stochastic gene expression module that
    require special handling.
    """
    def __init__(self, data, exceptions):
        self.data = data # A particular dataset that is being handled !!! Maybe too vague...
        self.exceptions = exceptions
        self.threshold_handler = ThresholdHandler(model, exceptions.thresholds)  # Use ThresholdHandler
        self._exception_handler()

    def _exception_handler(self):
        # Define a mapping of exception keys to their corresponding handler functions
        exception_handlers = {
            'Exp GCN': self._handle_exp_gcn,  
            'Exp RNA': self._set_mrna_copy_numbers,
            'kGin': self._handle_kgin,
            'kGac': self._handle_kgac,
            'kTCleak': self._handle_ktcleak,
            'kTCmaxs': self._handle_ktcmaxs,
            'kTCd': self._handle_ktcd,
            'Exp Protein': self._handle_exp_protein
        }

        # Iterate through attributes and call the corresponding handler if it exists
        for attr in self.exceptions.omics:
            handler = exception_handlers.get(attr.key())  # Get the handler function, if any
            if handler:
                handler()

    def _set_mrna_copy_numbers(self):
        """for handling instances of the mRNA Copy Number attribute in the SINGE model"""

        for entity, value in self.exceptions.omics.attr:
            index = self.data.omics.mRCN.to_list().index(entity)
            self._setattr(self.data.omics.mRCN[index], value)

    def _handle_exp_gcn(self):
        """for handling instances of the Gene Copy Number attribute in the SINGE model"""

        for entity, value in self.exceptions.omics.attr:
            index = self.data.omics.GCN.to_list().index(entity)
            self._setattr(self.data.omics.GCN[index], value)
    
    def _handle_kgin(self):
        """for handling instances of the Gene Inactivation Rate attribute in the SINGE model"""

        for entity, value in self.exceptions.genereg.attr:
            index = self.data.genereg.kGin.to_list().index(entity)
            self._setattr(self.data.genereg.kGin[index], value)

    def _handle_kgac(self):
        """for handling instances of the Gene Activation Rate attribute in the SINGE model"""

        for entity, value in self.exceptions.genereg.attr:
            index = self.data.genereg.kGac.to_list().index(entity)
            self._setattr(self.data.genereg.kGac[index], value)

    def _handle_ktcleak(self):
        """for handling instances of the Transcriptional Cleakage Rate attribute in the SINGE model"""

        for entity, value in self.exceptions.genereg.attr:
            index = self.data.genereg.kTCleak.to_list().index(entity)
            self._setattr(self.data.genereg.kTCleak[index], value)

    def _handle_ktcmaxs(self):
        """for handling instances of the Maximum Transcription Rate attribute in the SINGE model"""

        for entity, value in self.exceptions.genereg.attr:
            index = self.data.genereg.kTCmaxs.to_list().index(entity)
            self._setattr(self.data.genereg.kTCmaxs[index], value)

    def _handle_ktcd(self):
        """for handling instances of the Transcriptional Decay Rate attribute in the SINGE model"""

        for entity, value in self.exceptions.genereg.attr:
            index = self.data.genereg.kTCd.to_list().index(entity)
            self._setattr(self.data.genereg.kTCd[index], value)

    def _handle_exp_protein(self):
        """for handling instances of the Protein Copy Number attribute in the SINGE model"""

        for entity, value in self.exceptions.protein.attr:
            index = self.data.protein.pRCN.to_list().index(entity)
            self._setattr(self.data.protein.pRCN[index], value)

    def _setattr(self, attr, value):
        """
        Sets an attribute of a class to a given value.
        """
        setattr(self, attr, value)
