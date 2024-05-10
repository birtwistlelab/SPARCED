#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from pathlib import Path


def append_subfolder(folder: str | os.PathLike, subfolder: str) -> tuple[str | os.PathLike, bool]:
    """Append a subfolder to a folder path

    Arguments:
        folder: The folder path.
        subfolder: The subfolder name.

    Returns:
        A tuple composed of the subfolder's path and a boolean stating if an
        error occured.
    """

    error_occured = False
    folder = Path(folder)
    
    try:
        assert folder.exists()
    except:
        print("WARNING: Folder doesn't exist.\nFolder name:{name}."
              .format(name=folder))
        error_occured = True
    
    subfolder_path = folder / subfolder
    
    try:
        assert subfolder_path.exists()
    except:
        print("WARNING: Subfolder doesn't exist.\nSubfolder name: {name}."
              .format(name=subfolder))
        error_occured = True

    return(subfolder_path, error_occured)

