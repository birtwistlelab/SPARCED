#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import IO

import numpy as np
import pandas as pd
import re


def antimony_write_compartments_IC(f_antimony: IO[str], compartments: np.ndarray) -> None:
    """ Write compartments initial conditions in the given Antimony file

    Note:
        The first row is considered as a header, and hence it is skipped.
        Names should be located on the first column of the array.
        Volumes should be located on the second column of the array.

    Arguments:
        f_antimony: The open Antimony file.
        compartments: The content of the input compartments file.

    Returns:
        Nothing.
    """

    f_antimony.write("# Compartments initialization:\n")
    for i, value in enumerate(compartments[1:]):
        f_antimony.write("{name} = {volume:.6e};\n{name} has volume;\n"
                         .format(name=value[0], volume=np.double(value[1])))
    f_antimony.write("\n")

def antimony_write_reactions_IC(f_antimony: IO[str], p_names: np.ndarray, p_vals: np.ndarray) -> None:
    """Write reactions parameters initial conditions in the given Antimony file

    Warning:
        TODO use only one array for the parameters instead of taking the risk to
        separate names from values.

    Arguments:
        f_antimony: The open Antimony file.
        p_names: The parameters names.
        p_values: The parameters values.

    Returns:
        Nothing.
    """

    f_antimony.write("# Parameters initialization:\n ")
    for i, val in enumerate(p_names):
        f_antimony.write("{name} = {value:.6e};\n"
                         .format(name=val, value=np.double(p_vals[i])))
    f_antimony.write("\n")

def antimony_write_species_IC(f_antimony: IO[str], species: np.ndarray) -> None:
    """Write species initial concentrations in the given Antimony file

    Note:
        The first row is considered as a header, and hence it is skipped.
        Species names should be located on the first column of the array.
        Species concentrations should be located on the third column of the array.

    Arguments:
        f_antimony: The open Antimony file.
        species: Content of the input species file.

    Returns:
        Nothing.
    """

    f_antimony.write("# Species initialization:\n")
    for i, value in enumerate(species[1:]):
        f_antimony.write("{name} = {concentration:.6e};\n"
                         .format(name=value[0], concentration=np.double(value[2])))
    f_antimony.write("\n")

