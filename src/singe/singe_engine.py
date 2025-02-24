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
from types import SimpleNamespace
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
    def __init__(self, gene_regulation, omics_data, solver_flag, 
                 duration, exchange, model_exceptions):
        """
        Initialize the SINGEEngine class.
        """
        self.solver_flag = solver_flag
        self.model_exceptions = model_exceptions
        self.model = SimpleNamespace()

        self.prep = SINGEPrep(gene_regulation, 
                              omics_data, 
                              self.solver_flag,
                              duration,
                              exchange,
                              self.model_exceptions)

    def update_model(self, species_)
        """
        """

        active_genes, inactive_genes = self.convertmRNA2mpc()

        mrna_x = self._get_mrna_conc(species_concentrations, mrna_indicies)

    def simulate(self, species_concentrations):
        """
        Simulate the model with stochastic gene expression.
        """
        temp_params = SimpleNamespace()

        self.model.updated_gene_state_data = self._update_gene_data()

        self.model.updated_gene_vector = self._update_gene_vector()

        self.model.updated_mrnas = self._update_mrna_concentrations()
        
    def convertmRNA2mpc(self):
        """Converts nanomolar mRNA to mpc for tau-leap."""
        init_state = 0 
        active_genes = self.prep.gene_state_data[init_state:init_state+self.prep.number_of_genes]
        
        init_state = init_state+self.prep.number_of_genes

        inactive_genes = self.prep.gene_state_data[init_state:init_state+self.prep.number_of_genes]
        
        return active_genes, inactive_genes
    
    def _get_mrna_conc(self, mrna_indicies):
        """calculates the current mRNA concentration for the current iteration."""
        species_conc = self.prep.sbml_model.getInitialConcentrations()

        return np.divide(species_conc[mrna_indicies:], self.prep.mpc2nmcf_vn)
    
    def _update_mrna_concentrations(self, mrna_x, new_births, new_deaths):
        """Updates the mRNA concentrations for a particular iteration."""
        new_mrna_x = mrna_x + new_births - new_deaths

        new_mrna_x[new_mrna_x<0.0] = 0.0

        return new_mrna_x * self.prep.mpc2nmcf_vc
    
    def _make_tar_array(self):
        """Makes an array of values for the Transcriptional Activators and 
        Repressors."""
        species_conc = self.prep.sbml_model.getInitialConcentrations()

        spIDs = []

        for tar in range(self.prep.model.genereg.number_of_TARs):
            species_names = self.prep.model.sbml_model.get_species_ids()
            sps = species_names.index(self.prep.model.genereg.TARs[tar])

            spIDs.append(sps)

        return np.array(species_conc[spIDs])
    