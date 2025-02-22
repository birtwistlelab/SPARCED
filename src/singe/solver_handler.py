#!/bin/bash python3 
"""
Filename: solver_handler.py
Created: 2025-02-21
Author(s): Jonah R. Huggins

Description: SINGE - (S)tochastic (I)ntegrated (N)etwork for (G)ene (E)xpression
                singe is the stochastic gene expression simulation engine. It is a custom
                implementation of stochastic gene expression, written as the stochastic half
                of the SPARCED algorithm.
                
                This file specifies a class of functions that manage the solver instead for 
                SPARCED, so that users can specify new instances of solvers. 
"""
from types import SimpleNamespace

import numpy as np

class SolverHandler:
    """Parent class for handling gene expression decision flag."""
    GENE_ON = 1.0

    def __init__(self, solver, input_data, singe_model):
        self.input_data = input_data
        self.singe_model = singe_model
        self.return_data = _get_gene_expression_solver(solver, input_data, self.singe_model)
        del self.input_data # input data already stored in higher level
        del self.singe_model #import info is returned in return_model

    def _makeGeneVector(self):
        """Must be implemented within child classes"""
        raise NotImplementedError("Subclasses must implement `makeGeneVector` method.")

    def _calcGeneStateData(self):
        """Must be implemented within child classes."""
        raise NotImplementedError("Subclasses must implement `calcGeneData` method.")
        

class DeterministicSolverHandler(SolverHandler):
    """Handles calculations for the SINGE engine if the hybrid flag is set to False"""
    
    def __init__(self, solver, input_data, singe_model):
        super().__init__(solver, input_data, singe_model)
        self._makeGeneVector()
        self._calcGeneStateData()

    def _makeGeneVector(self):
        """Makes a zeros-array of shape [number of genes,1] that will store the 
        states of each gene over time."""
        self.return_data.gene_state_vector = np.zeros(shape = (self.input_data.sum_of_genes, 1))

        on_gene_indeces = np.random.choice(self.input_data.sum_of_genes, 
                                           size=int(round(
                                               self.input_data.sum_of_genes 
                                               * self.singe_model.omics.kGac[0] #Should be dynamic
                                               / self.singe_model.omics.kGin[0] #Should be dynamic
                                               )
                                            ),
                                            replace = False
        )

        self.return_data.gene_state_vector[on_gene_indeces] = self.GENE_ON

    def _calcGeneStateData(self):
        """
        Calculates the initial active gene concentration (mpc) and 
        initial inactive gene concentration (mpc) for the Deterministic Setting
        """
        self.return_data.gene_state_data = []

        initial_active_genes = ((self.singe_model.omics.kGac * self.singe_model.omics.GCN)
                                   / (self.singe_model.omics.kGin + self.singe_model.omics.kGac))

        initial_inactive_genes = (self.singe_model.omics.GCN - initial_active_genes)

        self.return_data.gene_state_data = np.concatenate((initial_active_genes, 
                                                           initial_inactive_genes), axis=None)


class StochasticSolverHandler(SolverHandler):
    """Handles calculations for the SINGE engine if the hybrid flag is set to True"""

    def __init__(self, solver, input_data, singe_model):
        super().__init__(solver, input_data, singe_model)
        self._makeGeneVector()
        self._calcGeneStateData()

    def _makeGeneVector(self):
        """Makes a zeros-array of shape [number of genes,1]"""
        self.return_data.gene_state_vector = np.zeros(shape = (self.input_data.sum_of_genes, 1))


        on_gene_indeces = np.random.choice(self.input_data.sum_of_genes, 
                                           size=int(round(
                                               self.input_data.sum_of_genes 
                                               * self.singe_model.omics.kGac[0] #Should be dynamic
                                               / self.singe_model.omics.kGin[0] #Should be dynamic
                                               )
                                            ),
                                            replace = False
        )

        self.return_data.gene_state_vector[on_gene_indeces] = self.GENE_ON


    def _calcGeneStateData(self):
        """
        Calculates the initial active gene concentration (mpc) and 
        initial inactive gene concentration (mpc) for the Deterministic Setting
        """
        self.return_data.gene_state_data = []

        initial_active_genes = np.dot(self.input_data.gene_position_matrix,
                                    self.input_data.gene_state_vector).ravel()
        
        initial_inactive_genes = (self.singe_model.omics.GCN
                                    - self.input_data.initial_active_genes).ravel()

        self.return_data.gene_state_data = np.concatenate((initial_active_genes, 
                                                           initial_inactive_genes), 
                                                        axis = None
        )


@staticmethod
def _get_gene_expression_solver(solver, input_data, singe_model):
    """Factory function to select the correct solver type."""
    if solver == "False":
        return DeterministicSolverHandler(solver, input_data, singe_model)
    elif solver == "True":
        return StochasticSolverHandler(solver, input_data, singe_model)
    else:
        raise ValueError(f"Unknown solver type: {solver}")
