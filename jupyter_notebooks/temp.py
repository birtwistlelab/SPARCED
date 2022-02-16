import pandas as pd
import numpy as np
# import re
# import libsbml
import os
import sys
import importlib
import amici

# from scipy.signal import find_peaks

# import matplotlib as mpl

# from mpi4py import MPI




# mpl.rcParams['figure.dpi'] = 300


wd = str(os.getcwd()).replace("jupyter_notebooks","")




sbml_file = "SPARCED_au565.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.join(wd,model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()