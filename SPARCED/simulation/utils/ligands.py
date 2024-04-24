#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from utils.arguments import parse_args
from utils.data_handling import load_input_data_file


def basic_ligands(c_egf: float=0.0, c_ins: float=0.0) -> np.ndarray:
    """Basic ligands setup

    Use passed arguments to set EGF and insulin values.
    Set all remaining ligands (HER, HGF, PDGF, FGF and IGF) to zero.

    Warning:
        SPARCED symbols for ligands are hard-coded.

    Arguments:
        c_egf: The extracellular concentration of EGF (nM).
        c_ins: The extracellular concentration of insulin (nM).

    Returns:
        A numpy array containing ligands symbols and concentrations.
    """
    return(np.array([['E',      float(c_egf)],                  # EGF
                     ['H',      float(0.0)],                    # HER
                     ['HGF',    float(0.0)],                    # HGF
                     ['P',      float(0.0)],                    # PDGF
                     ['F',      float(0.0)],                    # FGF
                     ['I',      float(0.0)],                    # IGF
                     ['INS',    float(c_ins)]], dtype=object))  # Insulin

def basic_ligands_from_arguments() -> np.ndarray:
    """
    Call basic ligands function while passing arguments values
    """
    args = parse_args()
    return(basic_ligands(args.egf, args.insulin))

def load_input_ligands(f_input: str) -> np.ndarray:
    """Load ligands concentrations from the given input dat a file

    Load an input ligands concentrations file structured as tab separated,
    and remove the first column (human-readable) to keep only the second
    (SPARCED symbol) and third (concentration in nM) columns.

    Arguments:
        f_input: The input ligands concentrations file.

    Returns:
        A numpy array containing ligands symbols and concentrations.
    """

    data = load_input_data_file(f_input)
    # Convert concentrations from strings to floats
    for row in data:
        row[2] = float(row[2])
    return(data[:,1:])

