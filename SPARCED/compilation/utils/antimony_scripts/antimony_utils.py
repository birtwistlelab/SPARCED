#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import re


def extract_antimony_model_name(f_antimony: str) -> str:
    """Extract the model name from the given Antimony file path

    Note:
        The model name is considered to be preceded by the last "ant_" prefix
        found and followed by the first dot "." found in the file path string.

    Arguments:
        f_antimony: The name of the Antimony file.

    Returns:
        The extracted model name.
    """

    sub_file_name = f_antimony.split("ant_")
    try:
        assert len(sub_file_name) > 0
    except:
        print("ERROR: Antimony model name could not be extracted.\n Please make \
               sure you passed all the required arguments. Current value for \
               the file name is: {name}.\n".format(name=f_antimony))
        sys.exit(0)
    model_name = sub_file_name[len(sub_file_name)-1].split(".")[0]
    return(model_name)

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

