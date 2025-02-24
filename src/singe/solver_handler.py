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
import numpy as np

class SolverHandler:
    """Parent class for handling gene expression decision flag."""
    GENE_ON = 1.0

    def __init__(self, solver, input_params):
        self.input_params = input_params
        self.return_data = _get_gene_expression_solver(solver, self.input_params)
        del self.input_params #import info is returned in return_data

    def make_gene_vector(self):
        """Must be implemented within child classes"""
        raise NotImplementedError("Subclasses must implement `makeGeneVector` method.")

    def calc_gene_state_data(self):
        """Must be implemented within child classes."""
        raise NotImplementedError("Subclasses must implement `calcGeneData` method.")
        

class DeterministicSolverHandler(SolverHandler):
    """Handles calculations for the SINGE engine if the hybrid flag is set to False"""
    
    def __init__(self, solver, input_params):
        super().__init__(solver, input_params)

    def make_gene_vector(self):
        """Makes a zeros-array of shape [number of genes,1] that will store the 
        states of each gene over time."""
        gene_state_vector = np.zeros(shape = (self.input_params.sum_of_genes, 1))

        on_gene_indeces = np.random.choice(self.input_params.sum_of_genes, 
                                           size=int(round(
                                               self.input_params.sum_of_genes 
                                               * self.input_params.omics.kGac[0] #TODO: Should be dynamic
                                               / self.input_params.omics.kGin[0] #TODO: Should be dynamic
                                               )
                                            ),
                                            replace = False
        )

        gene_state_vector[on_gene_indeces] = self.GENE_ON

        return gene_state_vector

    def calc_gene_state_data(self):
        """
        Calculates the initial active gene concentration (mpc) and 
        initial inactive gene concentration (mpc) for the Deterministic Setting
        """
        gene_state_data = []

        initial_active_genes = ((self.input_params.omics.kGac * self.input_params.omics.GCN)
                                   / (self.input_params.omics.kGin + self.input_params.omics.kGac))

        initial_inactive_genes = (self.input_params.omics.GCN - initial_active_genes)

        gene_state_data = np.concatenate((initial_active_genes,
                                        initial_inactive_genes), axis=None)
        
        return gene_state_data


class StochasticSolverHandler(SolverHandler):
    """Handles calculations for the SINGE engine if the hybrid flag is set to True"""

    def __init__(self, solver, input_params):
        super().__init__(solver, input_params)

    def make_gene_vector(self):
        """Makes a zeros-array of shape [number of genes,1]"""
        gene_state_vector = np.zeros(shape = (self.input_params.sum_of_genes, 1))


        on_gene_indeces = np.random.choice(self.input_params.sum_of_genes, 
                                           size=int(round(
                                               self.input_params.sum_of_genes 
                                               * self.input_params.omics.kGac[0] #Should be dynamic
                                               / self.input_params.omics.kGin[0] #Should be dynamic
                                               )
                                            ),
                                            replace = False
        )

        gene_state_vector[on_gene_indeces] = self.GENE_ON

        return gene_state_vector

    def calc_gene_state_data(self):
        """
        Calculates the initial active gene concentration (mpc) and 
        initial inactive gene concentration (mpc) for the Deterministic Setting
        """
        gene_state_data = []

        initial_active_genes = np.dot(self.input_params.gene_position_matrix,
                                    self.input_params.gene_state_vector).ravel()
        
        initial_inactive_genes = (self.input_params.omics.GCN
                                    - self.input_params.initial_active_genes).ravel()

        gene_state_data = np.concatenate((initial_active_genes, 
                                        initial_inactive_genes), 
                                        axis = None
        )

        return gene_state_data


@staticmethod
def _get_gene_expression_solver(solver, input_params):
    """Factory function to select the correct solver type."""
    if solver == "False":
        return DeterministicSolverHandler(solver, input_params)
    elif solver == "True":
        return StochasticSolverHandler(solver, input_params)
    else:
        raise ValueError(f"Unknown solver type: {solver}")
