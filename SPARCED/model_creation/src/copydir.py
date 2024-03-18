#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Copy Directory"""

import sys
from distutils.dir_util import copy_tree

def copy_directory(src, dest):
    try:
        copy_tree(src, dest)
    except:
        print("Error: failed to copy directory")
