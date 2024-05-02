#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Protein-Ligand Binding Related Functions"""

from numpy import sqrt

def bound_fraction(kd, pt, lt):
    """"Compute a protein's bound fraction

    Args:
        kd: the constant of dissociation as a numerical
        pt: the total concentration of protein as a numerical
        lt: the total concentration of ligand as a numerical

    Returns:
        A float.
    """
    a = pt
    b = (kd + lt + pt)
    c = lt
    return((b - sqrt(b*b - 4*a*c)) / (2*a))

def koff(kd, kon):
    """Compute the koff

    Args:
        kd: the constant of dissociation as a numerical
        kon: the kon as a numerical

    Returns:
        A float.
    """
    return(kd * kon)
