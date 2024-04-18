#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


def load_input_data_file(f_input: str) -> np.ndarray:
    """Load the given input data file 

    Load an input data file structured as tab separated.

    Arguments:
        f_input: The input data file.

    Returns:
        A numpy array containing the data.
    """

    data = np.array([np.array(line.strip().split("\t"))
                    for line in open(f_input)], dtype="object")
    return(data)

