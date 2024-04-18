#!/usr/bin/env python
#-*- coding: utf-8 -*-

import pandas as pd


def convert_excel_to_tsv(f_excel: str) -> None:
    """Convert an Excel file to TSV (SPARCED's standard input format)

    This function creates a new .txt file at the same location than the passed
    Excel file.

    Warning:
        This is some old code written four years ago, it hasn't been tested since.

    Arguments:
        f_excel: The Excel sheet path.

    Returns:
        Nothing.
    """

    data = pd.read_excel(f_excel, header=0, index_col=0)
    data.to_csv((f_excel.split("."))[0] + ".txt", sep="\t")

