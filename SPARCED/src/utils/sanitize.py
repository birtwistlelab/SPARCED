#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


def sanitize_model_name(name: str) -> str:
    """Ensure conformity of model name

    Note:
        As Antimony cannot handle the dash ('-') character in a model name, any
        occurence of this character is replaced with an underscore ('_').

    Arguments:
        name: The model name.

    Returns:
        The corresponding conform model name.
    """

    # Ensure that a name is specified
    try:
        assert not name == None
    except:
        print("ERROR: Please specify a model name. Aborting now.")
        sys.exit(0)
    # Sanitize name
    model_name = name.replace('-', '_')
    return(model_name)

def sanitize_popsize(popsize: int, abort_on_error: bool=False) -> int:
    """Ensure validity of population size

    Arguments:
        popsize: The population size.

    Returns:
        The population size if valid.
    """

    try:
        assert int(popsize) > 0
    except:
        print("ERROR: Cell population size should be superior to zero (0). \
                Current population size is: {size}."
              .format(size = popsize))
        if abort_on_error:
            print("Aborting now.")
            sys.exit(0)
        else:
            print("Setting automatically population size to one (1).")
            popsize = 1
    return(int(popsize))

