#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests of the volcalc module"""

from SPARCED.toolkit import volcalc as stk

def test_adjust_ec_vol():
    assert stk.adjust_ec_vol(10, 2) == 5
