import sys
import os
sys.path.append(os.getcwd()[0:os.getcwd().rfind('/')]+'/bin')

import libsbml
import importlib
import amici
import amici.plotting
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import scipy.stats
import argparse

from modules.RunSPARCED import RunSPARCED

#%%

# SBML model we want to import
sbml_file = 'SPARCED_U87.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name
# used for filename to save
cellNumber = 0  

