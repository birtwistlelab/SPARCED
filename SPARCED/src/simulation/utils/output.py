#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd


def save_simulation_output(model, simulation_name: str, simulation_number: int,
                           output_sim: str, xoutS_all: np.ndarray,
                           xoutG_all: np.ndarray, tout_all: np.ndarray) -> None:
    """
    Save simulation results in three files: species, genes and time
    """

    # TODO: ensure that the last character of the output_sim parameter is a dash,
    # otherwise if verbose ask the user if he meant to add one
    # Species
    columnsS = [ele for ele in model.getStateIds()]
    condsSDF = pd.DataFrame(data = xoutS_all, columns = columnsS)
    condsSDF.to_csv(output_sim + simulation_name + '_S_' + str(simulation_number) + '.txt', sep="\t")
    condsSDF = None
    # Genes
    columnsG = [x for n, x in enumerate(columnsS) if 'm_' in x]
    columnsG = columnsG[1:] # Skip header
    resa = [sub.replace('m_', 'ag_') for sub in columnsG]
    resi = [sub.replace('m_', 'ig_') for sub in columnsG]
    columnsG2 = np.concatenate((resa, resi), axis=None)
    condsGDF = pd.DataFrame(data = xoutG_all, columns = columnsG2)
    condsGDF.to_csv(output_sim + simulation_name + '_G_' + str(simulation_number) + '.txt', sep="\t")
    condsGDF = None
    # Time
    np.savetxt(output_sim + simulation_name + '_T_' + str(simulation_number) + '.txt', tout_all, newline="\t", fmt="%s")

