#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from pathlib import Path
import sys


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
        print("WARNING: Folder doesn't exist. This is never normal.\nFolder \
               name:{name}.".format(name=folder))
        if abort_on_error:
            print("Aborting now.")
            sys.exit(0)

    subfolder_path = folder / subfolder
    
    try:
        assert subfolder_path.exists()
    except:
        print("WARNING: Subfolder doesn't exist. This is normal if you are
               creating a new subfolder.\nSubfolder name: {name}."
               .format(name=subfolder))
        if abort_on_error:
            print("Aborting now.")
            sys.exit(0)

    return(subfolder_path)

