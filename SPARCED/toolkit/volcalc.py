#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compartment's Volume Related Functions"""

def adjust_ec_vol(ec, nb_cell):
    """Compute the extracellular volume adjusted to one cell

    Args:
        ec: the total extracellular volume
        nb_cell: the number of cells in the plate

    Returns:
        A float.
    """
    return(ec / nb_cell)
