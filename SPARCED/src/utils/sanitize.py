#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


def sanitize_model_name(name: str) -> str:
    """Ensure conformity of model name

    Note:
        As Antimony cannot handle the dash ('-') character in a model name, any
        occurence of this character is replaced with an underscore ('_').

    Arguments:
        name: THe model name.

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

