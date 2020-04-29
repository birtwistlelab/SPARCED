#!/usr/bin/env python3

import pandas as pd
import argparse

def changeSpeciesVals(valString):
    file_data = None #initialize outside of if statement
    file_data = pd.read_csv("Species.txt", sep="\t", header=0, index_col=0, encoding="latin-1")
    for item in valString.split(','):
        paramName, paramVal = tuple(item.split(":"))
        file_data.at[paramName, 'IC_Xinitialized'] =  paramVal
    file_data.to_csv("Species.txt", sep="\t")


parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--speciesVals', metavar='speciesVals', help='values to change in the species input data file')
args = parser.parse_args()

directive = None
with open(args.paramFile,"r") as f:
    directive = f.readline().strip()
    
changeSpeciesVals(directive)
