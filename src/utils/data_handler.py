#!/bin/bash python3
# -*- coding: utf-8 -*-

"""
Filename:
Author(s):
Created:
Description:
"""
#<-----------------------------Import Packages------------------------------->
import os
from types import SimpleNamespace

import pandas as pd
import numpy as np

#<-----------------------------Parent Class------------------------------->
class DataHandler:
    """
    Loads and performs CRUD operations on data.
    """
    def __init__(self, data_path: str):
        self.data_path = data_path
        self.data = None  # Placeholder for the actual data object

    def load_data(self):
        """Method to load data - must be implemented in child classes."""
        raise NotImplementedError("Subclasses must implement this method.")

    def modify_data(self, data):
        """Modify data - must be implemented in child classes."""
        raise NotImplementedError("Subclasses must implement this method.")

    def save_data(self):
        """Save data - must be implemented in child classes."""
        raise NotImplementedError("Subclasses must implement this method.")
    
#<-----------------------------Child Classes------------------------------->
class Perturbations(DataHandler):
    """
    Loads and performs CRUD operations on perturbations table. Perturbations table 
    consists of a single row of perturbed model entities, formated according to PEtab
    conditions file standards. 
    """
    def __init__(self, data_path: str):
        super().__init__(data_path)
        self.load_data()

    def load_data(self):
        """Load perturbations data."""
        self.data = pd.read_csv(self.data_path, sep='\t')

    def modify_data(self, data):
        """Modify perturbations data."""
        self.data = data

    def save_data(self):
        """Save perturbations data."""
        self.data.to_csv(self.data_path, sep='\t', index=False)

class Omics(DataHandler):
    """
    Loads and performs CRUD operations on omicsdata table. Omics data is a table of
    parameters pertaining to gene expression, mRNA copy number, protein abundance, 
    as well as degradation, synthesis, and translation rates.
    """
    def __init__(self, data_path: str):
        super().__init__(data_path)
        self.load_data()

    def load_data(self):
        """Load omics data."""
        omics_df = pd.read_csv(self.data_path, header = 0, index_col= 0, sep='\t')
        self.data = self._extract_omics_vals(omics_df)

    def modify_data(self, data):
        """Modify omics data."""
        self.data = data

    def getColumn(self, omics_df, cname):
        """
        Retrieves column as np.float64 vector array.
        """
        return np.array(omics_df[cname].values, dtype=np.float64)
    
    def _extract_omics_vals(self, omics_df):
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
            GCN = self.getColumn(omics_df, 'Exp GCN'),  #(G)ene (C)opy (N)umber in molecules per cell units
            mRCN = self.getColumn(omics_df,'Exp RNA'), # (mR)NA (C)opy (N)umber in molecules per cell units
            kGin = self.getColumn(omics_df,'kGin'), # Rate, (k), of (G)ene (in)activation
            kGac = self.getColumn(omics_df,'kGac'), # Rate, (k), of (G)ene (ac)tivation
            kTCleak = self.getColumn(omics_df,'kTCleak'), #Rate, (k), of (T)rans(C)riptional leakage
            kTCmaxs = self.getColumn(omics_df,'kTCmaxs'), #Rate, (k), of (T)rans(C)riptional (m)aximal production
            kTCd = self.getColumn(omics_df, 'kTCd'), #Rate, (k), of (T)rans(C)riptional (d)egradation
        )

    def save_data(self):
        """Save omics data."""
        self.data.to_csv(self.data_path, sep='\t', index=False)

class GeneRegulation(DataHandler):
    """
    Loads and performs CRUD operations on gene regulation data. Gene regulation data
    is a table of parameters informing the regulation of genes (rows) by individual 
    species (columns).
    """
    def __init__(self, data_path: str):
        super().__init__(data_path)
        self.load_data()

    def load_data(self):
        """Load gene regulation data."""
        genereg_df = pd.read_csv(self.data_path, header = 0, index_col= 0, sep='\t')
        self.data = self._extract_genereg_vals(genereg_df)

    def modify_data(self, data):
        """Modify gene regulation data."""
        self.data = data

    def getColumn(self, cname):
        """
        Retrieves column as np.float64 vector array.
        """
        return np.array(self.data[cname].values, dtype=np.float64)
    
    def _extract_genereg_vals(self, genereg_df):
        """
        retrieves the gene regulation data values as np.float64 arrays from data_handler.py, 
        """
        #(T)ranscriptional (A)ctivators & (R)epressors
        TARs = genereg_df.values

         # genereg file format makes column names TAR-species
        number_of_TARs = len(genereg_df.columns)

        return SimpleNamespace(TARs, number_of_TARs)

    def save_data(self):
        """Save gene regulation data."""
        self.data.to_csv(self.data_path, sep='\t', index=False)

class Results(DataHandler):
    """
    Performs Basic operations on the results trajectories provided by 
    the SPARCED algorithm.
    """
    def __init__(self, data_path: str):
        super().__init__(data_path)

    def load_data(self):
        """Load results data."""
        self.data = pd.read_csv(self.data_path, sep='\t')

    def modify_data(self, data):
        """Modify results data."""
        self.data = data

    def save_data(self, results: tuple):
        """Save results data."""

        sparced_results = self.combine_results(results)
        sparced_results.to_csv(self.data_path, sep='\t', index=False)

    def combine_results(self, results:tuple) -> pd.DataFrame:
        """
        Unpacks a tuple of results and combines them into a single dataframe.
        """
        sparced_results = pd.DataFrame()
        for data in results:
            if isinstance(data, pd.DataFrame):
                sparced_results = pd.concat([sparced_results, data], axis=1)
            elif isinstance(data, np.ndarray):
                sparced_results = pd.concat([sparced_results, pd.DataFrame(data)], axis=1)
            else:
                raise ValueError("Data must be a pandas dataframe or numpy array.")
        return sparced_results

# Storing here temporarily while SPARCED-specific code gets worked out. 
    # def combine_results(model, species: np.ndarray, 
    #                     genes: np.ndarray, time: np.ndarray) -> None:
    #     """
    #     Takes nested array results of genes, mRNA, proteins, and time values and 
    #     concatentates them into a single dataframe. 
    #     Parameters:
    #         - model: model object
    #         - xoutS_all (np.ndarray): array of species values
    #         - xoutG_all (np.ndarray): array of gene values
    #         - toutS_all (np.ndarray): array of time values
    #     Returns:
    #         - pd.DataFrame: dataframe of concatenated results
    #     """
    #     # Create a dataframe to store the results
    #     sparced_results = pd.DataFrame(data = time, columns = ['time'])

    #     # Get the gene and species names from the model
    #     species_names = model.getStateIds()
    #     gene_data = [x for n, x in enumerate(species_names) if 'm_' in x]

    #     # Add species to the results dictionary
    #     sparced_results = pd.concat([sparced_results, pd.DataFrame(data = species, columns = species_names)], axis=1)

    #     # Add genes to the results dictionary
    #     gene_data = gene_data[1:] # Skip header
    #     resa = [sub.replace('m_', 'ag_') for sub in gene_data]
    #     resi = [sub.replace('m_', 'ig_') for sub in gene_data]
    #     gene_data2 = np.concatenate((resa, resi), axis=None)

    #     sparced_results = pd.concat([sparced_results, pd.DataFrame(data = genes, columns = gene_data2)], axis=1)
        
    #     return sparced_results