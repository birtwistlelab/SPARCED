#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This script is an example of how to build and run the SPARCED model.

# Model creation
from model_creation.utils import *
from model_creation.create_model import *


# First step is to build the model.
# You may also want to comment this section out and run a model that was
# previously built instead.
# Note that building the model may take several minutes.
launch_model_creation() # Process parsed arguments and launch model creation

