#!/usr/bin/env python3

import pandas as pd
import argparse
import numpy as np
import csv

import glob

def changeSpeciesVals(valString):
    if len(valString) == 0:
        return
    file_data = None #initialize outside of if statement
    file_data = pd.read_csv("Species.txt", sep="\t", header=0, index_col=0, encoding="latin-1")
    for item in valString.split(','):
        paramName, paramVal = tuple(item.split(":"))
        file_data.at[paramName, 'IC_Xinitialized'] =  paramVal
    file_data.to_csv("Species.txt", sep="\t")


def changeRatelawVals(valString):
    if len(valString) == 0:
        return
    file_data = None #initialize outside of if statement
    file_data = np.array([np.array(line.strip().split("\t")) for line in open('Ratelaws.txt')])
    for item in valString.split(','):
        paramName, paramOffset, paramVal = tuple(item.split(":"))
        for idx,line in enumerate(file_data):
            if file_data[idx][0] == paramName:
                file_data[idx][2+int(paramOffset)] = paramVal
    with open('Ratelaws.txt', 'w') as f:
        csv.writer(f, delimiter="\t", lineterminator="\n").writerows(file_data)


# parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
# parser.add_argument('--paramfile', metavar='paramfile', help='file contains values to change in the species input data file')
# args = parser.parse_args()

paramfile = glob.glob('sweep*.txt')[0]
    

speciesDirective = None
ratelawDirective = None
with open(paramfile,"r") as f:
    speciesDirective = f.readline().strip()
    ratelawDirective = f.readline().strip()


changeSpeciesVals(speciesDirective)
changeRatelawVals(ratelawDirective)
