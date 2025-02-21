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
from types import SimpleNamespace
import pandas as pd
#<------------------------------Parent Class--------------------------------->
class SINGEPrep:
    """
    Prepares a stochastic gene expression model for simulation within the SPARCED algorithm.
    """
    def __init__(self, perturbed_model, solver, omics, genereg, exceptions = None):
        self.perturbed_model = perturbed_model
        self.solver = solver
        self.model_exceptions = exceptions
        self.prep(omics, genereg)

    def prep(self, omics, genereg):
        """
        Prepares the first round of simulation for the SINGEEngine
        """
        
        # Cytoplasmic and nuclear volume handled by SBMLModelHandler. Let it do it's job!
        cytoplasm_volume = self.perturbed_model.sbml_model.getVolume("Cytoplasm")
        nuclear_volume = self.perturbed_model.sbml_model.getVolume("Nucleus")

        omics = self._extract_omics_vals(omics)
        genereg = self._extract_genereg_vals(genereg)

        number_of_genes = int(len(omics.GCN))

        ## Start here on model exceptions routine. 

        return cytoplasm_volume, nuclear_volume, omics

    def _extract_omics_vals(self, omics):
        """
        retrieves the omics data values as np.float64 arrays from data_handler.py, 
        returns them for easier calculations.

        Returns:
        GCN (np.float64): Gene Copy Number in molecules per cell units
        mRCN (np.float64): mRNA Copy Number in molecules per cell units
        kGin (np.float64): Rate of Gene inactivation
        kGac (np.float64): Rate of Gene activation
        kTCleak (np.float64): Rate of Transcriptional leakage
        kTCmaxs (np.float64): Rate of Transcriptional maximal production
        kTCd (np.float64): Rate of Transcriptional degradation 
        """
        return SimpleNamespace(
            GCN = omics.getColumn('Exp GCN'),  #(G)ene (C)opy (N)umber in molecules per cell units
            mRCN = omics.getColumn('Exp RNA'), # (mR)NA (C)opy (N)umber in molecules per cell units
            kGin = omics.getColumn('kGin'), # Rate, (k), of (G)ene (in)activation
            kGac = omics.getColumn('kGac'), # Rate, (k), of (G)ene (ac)tivation
            kTCleak = omics.getColumn('kTCleak'), #Rate, (k), of (T)rans(C)riptional leakage
            kTCmaxs = omics.getColumn('kTCmaxs'), #Rate, (k), of (T)rans(C)riptional (m)aximal production
            kTCd = omics.getColumn('kTCd'), #Rate, (k), of (T)rans(C)riptional (d)egradation
        )
    
    def _extract_genereg_vals(self, genereg):
        """
        retrieves the gene regulation data values as np.float64 arrays from data_handler.py, 
        """
        #(T)rans(C)riptional (A)ctivators & (R)epressors
        TARs = genereg.values

         # genereg file format makes column names TAR-species
        number_of_TARs = len(genereg.columns)

        species_names = [name for name in self.perturbed_model.model.getStateIds()]

        species_indices = [species_names.index(name) for name in TARs.columns]

        return TARs, number_of_TARs, species_indices