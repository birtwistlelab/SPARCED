#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from utils.data_handling import load_input_data_file


def load_species_initial_conditions(f_species: str, ligands: np.ndarray,
                                    compound: str=None, dose: float=None) -> np.ndarray:
    """
    Create species initial conditions array
    """
    # Warning: first line is considered as a header and hence skipped
    # Species names should be on first column
    # Species initial conditions, should be on third column
    species = load_input_data_file(f_species)
    species_names = []
    species_initial_conditions = []
    for row in species[1:]:
        species_names = np.append(species_names, str(row[0]))
        species_initial_conditions = np.append(species_initial_conditions, float(row[2]))
    # Any concentration bellow 1e-6 is considered as zero (0)
    species_initial_conditions[np.argwhere(species_initial_conditions <= 1e-6)] = 0.0
    # Adjust ligands concentration
    for l in ligands:
        species_initial_conditions[np.argwhere(species_names == l[0])] = l[1]
    if compound is not None:
        species_initial_conditions[np.argwhere(species_names == compound)] = float(dose)
    return(species_initial_conditions)

