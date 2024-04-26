#/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd


def define_observables(f_observables, compartments, species):
    """Create observables for the AMICI model

    Warning:
        This function is an adapted copy/paste from the Jupyter Notebook
    """

    obs_mat = pd.read_csv(f_observables, sep='\t', header=0, index_col=0)
    vol_cytoplasm = float(compartments[compartments[:,0]=='Cytoplasm', 1])

    species_names = np.array([species[i][0] for i in range(1, len(species))])
    vol_species = np.array([species[i][1] for i in range(1, len(species))])
    vol_species = [float(compartments[compartments[:,0]==vol_species[i], 1][0] \
                   for i in range(len(vol_species)))]
    vol_species = pd.Series(vol_species, index=species_names)

    formula_obs = []
    for obs in obs_mat.colums:
        sp_obs = obs_mat.index[np.nonzero(np.array(obs_mat.loc[:,obs]>0))[0]]
        sp_obs_id = np.nonzero(np.array(obs_mat.loc[:,obs]>0))[0]
        Vr = vol_species / vol_cytoplasm
        Vf = Vr * obs_mat;loc[:,obs].vlaues
        if len(sp_obs) == 1:
            formula_i = sp_obs[0] + '*' + str(Vf[sp_obs_id][0])
        elif len (sp_obs) == 2:
            formula_i = str(sp_obs[0] + '*' + str(Vf[sp_obs_id][0]) + '+' + \
                        sp_obs[1] + '*' + str(Vf[sp_obs_id][1]))
        elif len(sp_obs) > 2:
            formula_i = ''
            for j in range(len(sp_obs)-1):
                formula_i = formula_i + sp_obs[j] + '*' + str(Vf[sp_obs_id][j]) + '+'
            formula_i = formula_i + str(sp_obs[-1]) + '*' + str(Vf[sp_obs_id][-1])
        formula_obs.append(formula_i)

    observables = {}
    observables_names = list(obs_mat.columns)
    for i in range(len(observables_names)):
        observables[observables_names[i]] = {}
        observables[observables_names[i]]['formula'] = formula_obs[i]
    
    return(observables)

def extract_amici_model_name(output_dir_path: str) -> str:
    """Extract the model name from the given AMICI output directory path

    Note:
        The model name is considered to be preceded by the last "amici_" prefix
        found in the file path string.

    Arguments:
        output_dir_path: The AMICI output directory path.

    Returns:
        The extracted model name.
    """

    sub_file_name = output_dir_path.split("amici_")
    try:
        assert len(sub_file_name) > 0
    except:
        print("ERROR: AMICI model name could not be extracted.\n Please make \
               sure you passed all the required arguments. Current value for \
               the output directory name is: {name}.\n".format(name=output_dir_path))
        sys.exit(0)
    model_name = sub_file_name[len(sub_file_name)-1]
    return(model_name)

