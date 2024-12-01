#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests of the protlig module"""

from SPARCED.toolkit import protlig as stk

def test_bound_fraction():
    # Tests from  https://jiang.bio.purdue.edu/binding
    assert (stk.bound_fraction(1.0, 0.1, 12.402) >= 0.92) and (stk.bound_fraction(1.0, 0.1, 12.402) <= 0.93)
    assert (stk.bound_fraction(1.0, 1.0, 1.685) >= 0.53) and (stk.bound_fraction(1.0, 1.0,  1.685) <= 0.54)

def test_koff():
    assert stk.koff(3, 4) == 12
