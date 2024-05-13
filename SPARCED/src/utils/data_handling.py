#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

import numpy as np
import pandas as pd
from pathlib import Path
import yaml

def append_subfolder(folder: str | os.PathLike, subfolder: str,
                     abort_on_error: bool=False) -> str | os.PathLike:          
    """Append a subfolder to a folder path

    Arguments:
        folder: The folder path.
        subfolder: The subfolder name.
        abort_on_error: Abort process when encountering an error.

    Returns:
        The subfolder's path.
    """

    folder = Path(folder)
    try:
        assert folder.exists()
    except:
        print("WARNING: Folder doesn't exist. This is never normal.\nFolder name:{name}."
              .format(name=folder))
        if abort_on_error:
            print("Aborting now.")
            sys.exit(0)

    subfolder_path = folder / subfolder
    try:
        assert subfolder_path.exists()
    except:
        print("WARNING: Subfolder doesn't exist yet. This is normal if you are creating a new subfolder.\nSubfolder name: {name}."
              .format(name=subfolder))
        if abort_on_error:
            print("Aborting now.")
            sys.exit(0)

    return(subfolder_path)

def convert_excel_to_tsv(f_excel: str) -> None:
    """Convert an Excel file to TSV (SPARCED's standard input format)           

    This function creates a new .txt file at the same location than the passed  
    Excel file.

    Warning:
        This is some old code written four years ago, it hasn't been tested since.

    Arguments:
        f_excel: The Excel sheet path.

    Returns:
        Nothing.
    """

    data = pd.read_excel(f_excel, header=0, index_col=0)
    data.to_csv((f_excel.split("."))[0] + ".txt", sep="\t") 

def load_input_data_config(data_path: str | os.PathLike, yaml_name: str) -> dict[str, str | os.PathLike]:
    """Load input data files paths configuration

    Note:
        File structure is assumed to be organized as follow:
        > model folder
        > data subfolder containing a YAML configuration file describing input
        data organization
        > model compilation and simulation sub-subfolders containing the input
        data files

    Arguments:
        data_path: The input data files folder path.
        yaml_name: The YAML configuration file name.

    Returns:
        A dictionnary containing all the input data file paths.
    """

    # Load data and YAML paths
    yaml_path = append_subfolder(data_path, yaml_name, True)
    # Read input data files structure in YAML configuration file
    with yaml_path.open() as f:
        input_files_configuration = yaml.safe_load(f)
    return(input_files_configuration)

def load_input_data_file(f_input: str | os.PathLike) -> np.ndarray:
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

