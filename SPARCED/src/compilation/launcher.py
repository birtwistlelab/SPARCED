#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import yaml

from compilation.create_model import create_model
from utils.arguments import parse_args
from utils.paths_handling import append_subfolder


def launch_model_creation() -> None:
    """Launch model creation function

    Small routine to process the parsed arguments and call create_model()

    Arguments:
        None.
    
    Returns:
        Nothing.
    """

    args = parse_args()
    # Runtime booleans
    is_SPARCED = not args.wild  # if its not wild then it's SPARCED
    verbose = args.verbose
    # Model
    model_name = load_model_name(args.name)
    model_path = load_model_path(args.model, model_name)
    # Input data files
    input_files = load_compilation_input_files(model_path, args.input_data, args.yaml)
    # Call create_model
    create_model(model_name, model_path, input_files, args.output_parameters,
                 verbose, is_SPARCED)

def load_compilation_input_files(model_path: str | os.PathLike, data_folder: str,
                                 yaml_name: str) -> dict[str, str | os.PathLike]:
    """Load compilation input data files paths dictionnary

    Note:
        File structure is assumed to be organized as follow:
        > model folder
        > data subfolder with YAML configuration file describing input data
        organization
        > model compilation sub-subfolder containing the input data files
        required for model compilation

    Arguments:
        model_path: The model folder path.
        data_folder: The input data files folder name.
        yaml_name: The YAML configuration file name.

    Returns:
        A dictionnary containing all the input data file paths.
    """

    # Load data and YAML paths
    data_path, error = append_subfolder(model_path, data_folder)
    yaml_path, error = append_subfolder(data_path, yaml_name)
    try:
        assert not error
    except:
        print("ERROR: Data or YAML path doesn't exist. Aborting now.")
        sys.exit(0)
    # Read input data files structure in YAML configuration file
    with yaml_path.open() as f:
        input_files_structure = yaml.safe_load(f)
    compilation_files = input_files_structure["compilation"]
    # Load root of compilation input data files
    compilation_root = compilation_files.pop("root", None)
    compilation_data_path, error = append_subfolder(data_path, compilation_root)
    try:
        assert not error
    except:
        print("ERROR: Compilation data input files folder doesn't exist. Aborting now.")
        sys.exit(0)
    # Reconstruct full path for each input file listed in the dictionnary
    for input_file in compilation_files.keys():
        if compilation_files[input_file] != None:
            compilation_files[input_file], error = append_subfolder(compilation_data_path, compilation_files[input_file])
    return(compilation_files)

def load_model_name(name: str) -> str:
    """Ensure conformity of model name

    Note:
        As Antimony cannot handle the dash ('-') character in a model name, any
        occurence of this character is replace with an underscore ('_').

    Arguments:
        name: The model name.

    Returns:
        The conform model name.
    """

    # Ensure that a name is specified
    try:
        assert not name == None
    except:
        print("ERROR: Please specify a model name. Aborting now.")
        sys.exit(0)
    # Clean name
    model_name = name.replace('-', '_')
    return(model_name)

def load_model_path(models_folder: str, model_name: str) -> str | os.PathLike:
    """Ensure correctness of model path

    Arguments:
        models_folder: The path to the models folder.
        model_name: The conform model name.

    Returns:
        The correct model's resulting path.
    """
    
    model_path, error = append_subfolder(models_folder, model_name)
    try:
        assert not error
    except:
        print("{name}: Model path doesn't exist. Aborting now."
               .format(name=model_name))
        sys.exit(0)
    return(model_path)


if __name__ == '__main__':
    launch_model_creation()

