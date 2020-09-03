import pandas as pd
import sys

def excelToTSV(fname):
    file_data = pd.read_excel(fname, header=0, index_col=0)
    file_data.to_csv((fname.split("."))[0] + ".txt", sep="\t")

excelToTSV(sys.argv[1])
