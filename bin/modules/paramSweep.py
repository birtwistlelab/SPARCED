import pandas as pd

def paramSweep(sweepString):
    fname, rowName, colName, paramVal = tuple(sweepString.split(":"))
    file_data = pd.read_csv(fname+".txt", sep="\t", header=0, index_col=0, encoding="latin-1")
    file_data.at[rowName, colName] =  paramVal
    file_data.to_csv(fname+".txt", sep="\t")
