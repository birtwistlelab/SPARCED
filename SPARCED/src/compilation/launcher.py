#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

from compilation.create_model import create_model
from utils.arguments import parse_args
from utils.data_handling import *
from utils.sanitize import sanitize_model_name


def launch_model_creation() -> None:
    """Launch model creation function

    Small routine to process the parsed arguments and call create_model()

    Arguments:
        None.
    
    Returns:
        Nothing.
    """

    args = parse_args()
    # Model
    model_name = sanitize_model_name(args.name)
    model_path = append_subfolder(args.model, model_name, True)
    # Input data files
    input_folder = append_subfolder(model_path, args.input_data, True)
    input_files = load_compilation_input_files(input_folder, args.yaml)
    # Output parameters
    output_parameters_path = append_subfolder(model_path, args.output_parameters)
    # Runtime booleans
    is_SPARCEd = not args.wild  # if it's not wild then it's SPARCED
    verbose = args.verbose
    # Call create_model
    create_model(model_name, model_path, input_files, output_parameters_path,
                 verbose, is_SPARCED)

def load_compilation_input_files(data_folder: str | os.PathLike, yaml_name: str) -> dict[str, str | os.PathLike]:
    """Load compilation input data files paths dictionnary

    Arguments:
        model_path: The model folder path.
        data_folder: The input data files folder name.
        yaml_name: The YAML configuration file name.

    Returns:
        A dictionnary containing all the input data file paths.
    """

    input_files_configuration = load_input_data_config(data_folder, yaml_name)
    compilation_files = input_files_configuration["compilation"]
    # Load root of compilation input data files
    compilation_root = compilation_files.pop("root", None)
    compilation_data_path = append_subfolder(data_folder, compilation_root, True)
    # Reconstruct full path for each input file listed in the dictionnary
    for input_file in compilation_files.keys():
        if compilation_files[input_file] != None:
            compilation_files[input_file] = append_subfolder(compilation_data_path, compilation_files[input_file])
    return(compilation_files)


if __name__ == '__main__':
    launch_model_creation()

