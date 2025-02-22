#!/bin/bash python3
# -*- coding: utf-8 -*-
"""
Filename: prep.py
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
from types import SimpleNamespace

import numpy as np

from src.utils.model_handler import SINGEModelHandler
from solver_handler import SolverHandler
from bioexception_handler import ModelExceptions
#<------------------------------Parent Class--------------------------------->
class SINGEPrep:
    """
    Prepares a stochastic gene expression model for simulation within the SPARCED algorithm.
    """

    def __init__(self, solver, omics, genereg, duration, exchange, exceptions = None):
        self.prep = SINGEModelHandler((omics, genereg)) #Currently tuple so model handler doesn't break parent class.
        self.duration = duration
        self.exchange = exchange
        self.prep.Vn = 1.7500E-12 # Nuclear Volume
        self.prep.Vc = 5.2500E-12 # Cytoplasmic Volume
        self.solver = solver
        self.model_exceptions = exceptions

        solver_data = SimpleNamespace()

        solver_data.number_of_genes = int(len(self.prep.model.omics.GCN))

        solver_data.sum_of_genes = int(sum(self.prep.model.omics.GCN)) # sum of gene copy numbers

        self._makeGenePositionMatrix(solver_data.number_of_genes,
                                    solver_data.sum_of_genes)

        return_data = SolverHandler(self.solver, solver_data, self.prep) # Handles solver-flag related tasks

        self.prep.gene_state_vector, self.prep.gene_state_data = (return_data.gene_state_vector, 
                                                                return_data.gene_state_data)
        
        self._staticGeneActivationRate()
        self._staticGeneInactivationRate()
        self._makeTARTrajectories(solver_data.number_of_genes)
        self._makeEmptyResultsMatrix

        self.prep = ModelExceptions(self.prep, self.model_exceptions)

        ## Start here on model exceptions routine. 
        ### here, mRNA species for cellcycle are turned to 17
        # mExp_mpc[indsDm] = 17.0 # modify cell cycle gene mRNA numbers to 17
    
    def _makeGenePositionMatrix(self, number_of_genes, sum_of_genes):
        """
        builds a positional matrix of gene locations using the numbers of genes provided.
        """
        index = 0

        # Create a template matrix for storing gene positions
        gene_position_matrix = np.zeros((number_of_genes, sum_of_genes))

        for gene in range(number_of_genes):
            gene_position_matrix[gene, index:index+int(self.prep.model.omics.GCN[gene])] = 1.0

            index = index + int(self.prep.model.omics.GCN[gene])

        self.prep.gene_position_matrix = gene_position_matrix

    def _staticGeneActivationRate(self):
        """For now, the SPARCED model only uses one instance of the gene activation rate
        this method drops the activation rate vector (singe_model.omics.kGac) and returns
        only the first instance in the vector"""
        self.prep.kGac = self.prep.model.omics.kGac[0]
    
    def _staticGeneInactivationRate(self):
        """For now, the SPARCED model only uses one instance of the gene inactivation rate
        this method drops the inactivation rate vector (singe_model.omics.kGac) and returns
        only the first instance in the vector"""
        self.prep.kGin = self.prep.model.omics.kGin[0]

    def _makeTARTrajectoreis(self, number_of_genes):
        """Makes a [number_of_genes, number_of_TARs] shape array for ... [Talk to marc]"""
        self.prep.tcnas = np.ones((number_of_genes, self.prep.model.genereg.number_of_TARs))
        self.prep.tck50as = np.zeros((number_of_genes, self.prep.model.genereg.number_of_TARs))
        self.prep.tcrs = np.zeros((number_of_genes, self.prep.model.genereg.number_of_TARs))
        self.prep.tck50rs = np.zeros((number_of_genes, self.prep.model.genereg.number_of_TARs))

        for gene in range(number_of_genes):
            for TAR in range(self.prep.model.genereg.number_of_TARs)

            partial_ARs = self.prep.model.genereg.TARs[gene, TAR].find(';')

            if partial_ARs > 0:
                nH = np.float(self.prep.model.genereg.TARs[gene,TAR][0:partial_ARs])
                kH = np.float(self.prep.model.genereg.TARs[gene,TAR][partial_ARs+2::])
                if nH>0:
                    self.prep.tcnas[gene,TAR] = nH
                    self.prep.tck50as[gene,TAR] = kH
                else:
                    self.prep.tcnrs[gene,TAR] = abs(nH)
                    self.prep.tck50rs[gene,TAR] = kH

            self.prep.nanomoles = 1.0E9
            self.prep.AVAGADRO = 6.023E+23

            mpc2nmcf_Vn = self.prep.nanomoles/(self.prep.Vn*self.prep.AVAGADRO)

            self.prep.tck50as = self.prep.tck50as*(1/mpc2nmcf_Vn)
            self.prep.tck50rs = self.prep.tck50rs*(1/mpc2nmcf_Vn)
    
    def _makeEmptyResultsMatrix(self):
        "Makes empty matrix for the results to be stored in"
        self.prep.genes = np.zeros(shape=((int(self.duration*3600/self.exchange)+1, 
                                            len(self.prep.gene_state_data))))
        
        self.prep.genes[0, :] = self.prep.gene_state_data
